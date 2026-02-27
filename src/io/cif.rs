use std::collections::HashMap;
use std::path::Path;

use crate::core::atom::AtomInfo;
use crate::core::element::element_by_symbol;
use crate::core::molecule::Molecule;
use crate::io::pdb::infer_bonds;

/// Parse an mmCIF file and return a Molecule.
pub fn load_cif(path: &Path) -> Result<Molecule, String> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| format!("Failed to read {}: {}", path.display(), e))?;
    parse_cif_string(&content, path)
}

/// Parse mmCIF content from a string.
pub fn parse_cif_string(content: &str, path: &Path) -> Result<Molecule, String> {
    let name = path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown")
        .to_string();

    let mut mol = Molecule::new(name);
    let mut tokens = Tokenizer::new(content);

    // Look for _atom_site loop
    while let Some(tok) = tokens.peek() {
        if tok == "loop_" {
            tokens.next(); // consume loop_
            if is_atom_site_loop(&mut tokens) {
                parse_atom_site_loop(&mut tokens, &mut mol)?;
            }
        } else {
            tokens.next();
        }
    }

    // Post-processing
    if mol.bonds.is_empty() {
        infer_bonds(&mut mol);
    }
    mol.build_residues();
    mol.apply_default_representation();

    log::info!(
        "Loaded mmCIF: {} atoms, {} bonds, {} residues",
        mol.atoms.len(),
        mol.bonds.len(),
        mol.residues.len()
    );

    Ok(mol)
}

fn is_atom_site_loop(tokens: &mut Tokenizer) -> bool {
    // Check if the next few keys start with _atom_site.
    // We peek without consuming.
    if let Some(key) = tokens.peek() {
        return key.starts_with("_atom_site.");
    }
    false
}

fn parse_atom_site_loop(tokens: &mut Tokenizer, mol: &mut Molecule) -> Result<(), String> {
    // 1. Read keys to build column mapping
    let mut columns: HashMap<String, usize> = HashMap::new();
    let mut col_idx = 0;

    while let Some(tok) = tokens.peek() {
        if tok.starts_with("_") {
            columns.insert(tok.to_string(), col_idx);
            tokens.next();
            col_idx += 1;
        } else {
            break; // End of header
        }
    }

    // Required columns
    let c_group = columns.get("_atom_site.group_PDB"); // ATOM/HETATM
    let c_symbol = columns.get("_atom_site.type_symbol");
    let c_atom = columns.get("_atom_site.label_atom_id"); // Name
    let c_resn = columns.get("_atom_site.label_comp_id");
    let c_chain = columns.get("_atom_site.auth_asym_id").or(columns.get("_atom_site.label_asym_id"));
    let c_resi = columns.get("_atom_site.auth_seq_id").or(columns.get("_atom_site.label_seq_id"));
    let c_x = columns.get("_atom_site.Cartn_x");
    let c_y = columns.get("_atom_site.Cartn_y");
    let c_z = columns.get("_atom_site.Cartn_z");
    let c_occ = columns.get("_atom_site.occupancy");
    let c_b = columns.get("_atom_site.B_iso_or_equiv");
    let c_id = columns.get("_atom_site.id"); // Serial

    if col_idx == 0 {
        return Err("No columns found in _atom_site loop".into());
    }

    let c_model = columns.get("_atom_site.pdbx_PDB_model_num"); // Model number

    if c_x.is_none() || c_y.is_none() || c_z.is_none() {
        return Err("Missing coordinate columns in _atom_site loop".into());
    }

    let mut row_count = 0;
    // 2. Read rows
    loop {
        // Check if we hit a new keyword or loop_ or end of file
        let peek = match tokens.peek() {
            Some(p) => p,
            None => break,
        };
        if peek.starts_with("_") || peek == "loop_" || peek == "data_" {
            break;
        }

        // Read one row (col_idx items)
        let mut row: Vec<&str> = Vec::with_capacity(col_idx);
        for _ in 0..col_idx {
            if let Some(val) = tokens.next() {
                row.push(val);
            } else {
                break;
            }
        }

        if row.len() != col_idx {
            break; // Incomplete row
        }

        row_count += 1;
        if row_count % 10000 == 0 {
             // println!("Parsed {} rows...", row_count);
        }

        // Parse row
        let group = c_group.map(|&i| row[i]).unwrap_or("ATOM");
        if group != "ATOM" && group != "HETATM" {
            continue;
        }
        
        // Only load the first model
        if let Some(model_idx) = c_model {
            let model_num: i32 = row[*model_idx].parse().unwrap_or(1);
            if model_num > 1 {
                continue;
            }
        }

        let is_hetatm = group == "HETATM";

        let x: f32 = c_x.and_then(|&i| row[i].parse().ok()).unwrap_or(0.0);
        let y: f32 = c_y.and_then(|&i| row[i].parse().ok()).unwrap_or(0.0);
        let z: f32 = c_z.and_then(|&i| row[i].parse().ok()).unwrap_or(0.0);

        let atom_name = c_atom.map(|&i| row[i]).unwrap_or("CA").to_string();
        let elem_sym = c_symbol.map(|&i| row[i]).unwrap_or("C").to_string();
        let resn = c_resn.map(|&i| row[i]).unwrap_or("UNK").to_string();
        let chain = c_chain.map(|&i| row[i].chars().next().unwrap_or(' ')).unwrap_or(' ');
        
        // resi might be '.'
        let resi_str = c_resi.map(|&i| row[i]).unwrap_or("0");
        let resi: i32 = if resi_str == "." || resi_str == "?" {
            0
        } else {
            resi_str.parse().unwrap_or(0)
        };

        let occupancy: f32 = c_occ.and_then(|&i| row[i].parse().ok()).unwrap_or(1.0);
        let b_factor: f32 = c_b.and_then(|&i| row[i].parse().ok()).unwrap_or(0.0);
        let serial: u32 = c_id.and_then(|&i| row[i].parse().ok()).unwrap_or(0);

        // Lookup element data
        let elem_data = element_by_symbol(&elem_sym);
        let atomic_number = elem_data
            .map(|e| {
                crate::core::element::ELEMENTS
                    .iter()
                    .position(|x| std::ptr::eq(x, e))
                    .unwrap_or(0) as u8
            })
            .unwrap_or(0);
        let vdw = elem_data.map(|e| e.vdw).unwrap_or(1.7);
        let color = elem_data.map(|e| e.color).unwrap_or([0.5, 0.5, 0.5]);

        let atom = AtomInfo {
            name: atom_name,
            element: atomic_number,
            elem_symbol: elem_sym,
            resn,
            resi,
            ins_code: '\0', // basic support
            chain,
            alt: ' ', // basic support
            b_factor,
            occupancy,
            vdw,
            color,
            is_hetatm,
            serial,
            ..Default::default()
        };

        mol.atoms.push(atom);
        mol.coords.push([x, y, z]);
    }

    Ok(())
}

// ── Simple Tokenizer ────────────────────────────────────────────────────────

struct Tokenizer<'a> {
    input: &'a str,
    cursor: usize,
    peeked: Option<&'a str>,
}

impl<'a> Tokenizer<'a> {
    fn new(input: &'a str) -> Self {
        let mut t = Self {
            input,
            cursor: 0,
            peeked: None,
        };
        t.advance();
        t
    }

    fn peek(&self) -> Option<&'a str> {
        self.peeked
    }

    fn next(&mut self) -> Option<&'a str> {
        let val = self.peeked;
        self.advance();
        val
    }

    fn advance(&mut self) {
        // Skip whitespace
        while self.cursor < self.input.len() {
            let c = self.input[self.cursor..].chars().next().unwrap();
            if c.is_whitespace() {
                self.cursor += c.len_utf8();
            } else if c == '#' {
                // Skip comment until newline
                while self.cursor < self.input.len() {
                    let nc = self.input[self.cursor..].chars().next().unwrap();
                    if nc == '\n' {
                        self.cursor += 1;
                        break;
                    }
                    self.cursor += nc.len_utf8();
                }
            } else {
                break;
            }
        }

        if self.cursor >= self.input.len() {
            self.peeked = None;
            return;
        }

        let start = self.cursor;
        let first = self.input[start..].chars().next().unwrap();

        if first == '\'' || first == '"' {
            // Quoted string
            let quote = first;
            self.cursor += 1; // skip open quote
            let mut end = self.cursor;
            while self.cursor < self.input.len() {
                let c = self.input[self.cursor..].chars().next().unwrap();
                if c == quote {
                    // Check if it's a closing quote (followed by whitespace or end)
                    let next_idx = self.cursor + 1;
                    let is_end = if next_idx >= self.input.len() {
                        true
                    } else {
                        let nc = self.input[next_idx..].chars().next().unwrap();
                        nc.is_whitespace()
                    };

                    if is_end {
                        end = self.cursor;
                        self.cursor += 1; // skip close quote
                        break;
                    }
                }
                self.cursor += c.len_utf8();
            }
            self.peeked = Some(&self.input[start + 1..end]);
        } else {
            // Unquoted string
            while self.cursor < self.input.len() {
                let c = self.input[self.cursor..].chars().next().unwrap();
                if c.is_whitespace() {
                    break;
                }
                self.cursor += c.len_utf8();
            }
            self.peeked = Some(&self.input[start..self.cursor]);
        }
    }
}
