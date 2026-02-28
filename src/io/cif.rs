use std::collections::HashMap;
use std::path::Path;

use crate::core::atom::AtomInfo;
use crate::core::element::{element_by_symbol, ELEMENTS};
use crate::core::molecule::Molecule;
use crate::io::pdb::infer_bonds;
use crate::io::LoadResult;

#[derive(Debug, Clone, PartialEq, Eq)]
struct AtomIdentity {
    is_hetatm: bool,
    chain: char,
    resn: String,
    resi: i32,
    ins_code: char,
    name: String,
    alt: char,
}

#[derive(Debug, Clone)]
struct ParsedAtom {
    identity: AtomIdentity,
    atom: AtomInfo,
    coord: [f32; 3],
}

/// Parse an mmCIF file and return a Molecule + loader metadata.
pub fn load_cif(path: &Path) -> Result<LoadResult, String> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| format!("Failed to read {}: {}", path.display(), e))?;
    parse_cif_string(&content, path)
}

/// Parse mmCIF content from a string.
pub fn parse_cif_string(content: &str, path: &Path) -> Result<LoadResult, String> {
    let name = path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown")
        .to_string();

    let mut mol = Molecule::new(name);
    let mut warnings = Vec::new();
    let mut tokens = Tokenizer::new(content);

    // State rows in encounter-order remapped model numbering (1..N).
    let mut state_rows: Vec<Vec<ParsedAtom>> = Vec::new();

    // Look for _atom_site loop(s)
    while let Some(tok) = tokens.peek() {
        if tok == "loop_" {
            tokens.next(); // consume loop_
            if is_atom_site_loop(&mut tokens) {
                parse_atom_site_loop(&mut tokens, &mut state_rows)?;
            }
        } else {
            tokens.next();
        }
    }

    let source_model_count = state_rows.len();

    if let Some(first_state) = state_rows.first() {
        mol.atoms = first_state.iter().map(|r| r.atom.clone()).collect();
        mol.coord_sets = vec![first_state.iter().map(|r| r.coord).collect()];

        let mut mismatch = false;
        for (state_idx, rows) in state_rows.iter().enumerate().skip(1) {
            if rows.len() != first_state.len() {
                warnings.push(format!(
                    "Model {} atom count {} differs from model 1 atom count {}; keeping model 1 only",
                    state_idx + 1,
                    rows.len(),
                    first_state.len()
                ));
                mismatch = true;
                break;
            }
            let mut bad_idx = None;
            for (i, (a, b)) in first_state.iter().zip(rows).enumerate() {
                if a.identity != b.identity {
                    bad_idx = Some(i + 1);
                    break;
                }
            }
            if let Some(row) = bad_idx {
                warnings.push(format!(
                    "Model {} topology/order differs at atom row {}; keeping model 1 only",
                    state_idx + 1,
                    row
                ));
                mismatch = true;
                break;
            }
            mol.coord_sets.push(rows.iter().map(|r| r.coord).collect());
        }

        if mismatch {
            mol.coord_sets.truncate(1);
        }
    } else {
        mol.coord_sets = vec![Vec::new()];
    }

    // Post-processing based on model/state 1 topology.
    if mol.bonds.is_empty() {
        infer_bonds(&mut mol);
    }
    mol.build_residues();
    mol.apply_default_representation();

    log::info!(
        "Loaded mmCIF: {} atoms, {} bonds, {} residues, {} state(s) from {} model(s)",
        mol.atoms.len(),
        mol.bonds.len(),
        mol.residues.len(),
        mol.state_count(),
        source_model_count
    );

    Ok(LoadResult {
        molecule: mol,
        warnings,
        source_model_count,
    })
}

fn is_atom_site_loop(tokens: &mut Tokenizer) -> bool {
    // Check if the next key starts with _atom_site.
    if let Some(key) = tokens.peek() {
        return key.starts_with("_atom_site.");
    }
    false
}

fn parse_atom_site_loop(
    tokens: &mut Tokenizer,
    states: &mut Vec<Vec<ParsedAtom>>,
) -> Result<(), String> {
    // 1. Read loop keys/columns.
    let mut columns: HashMap<String, usize> = HashMap::new();
    let mut col_idx = 0usize;

    while let Some(tok) = tokens.peek() {
        if tok.starts_with('_') {
            columns.insert(tok.to_string(), col_idx);
            tokens.next();
            col_idx += 1;
        } else {
            break;
        }
    }

    if col_idx == 0 {
        return Err("No columns found in _atom_site loop".into());
    }

    // 2. Resolve needed columns.
    let c_group = columns.get("_atom_site.group_PDB"); // ATOM/HETATM
    let c_symbol = columns.get("_atom_site.type_symbol");
    let c_atom = columns
        .get("_atom_site.auth_atom_id")
        .or(columns.get("_atom_site.label_atom_id"));
    let c_resn = columns
        .get("_atom_site.auth_comp_id")
        .or(columns.get("_atom_site.label_comp_id"));
    let c_chain = columns
        .get("_atom_site.auth_asym_id")
        .or(columns.get("_atom_site.label_asym_id"));
    let c_resi = columns
        .get("_atom_site.auth_seq_id")
        .or(columns.get("_atom_site.label_seq_id"));
    let c_ins = columns.get("_atom_site.pdbx_PDB_ins_code");
    let c_alt = columns.get("_atom_site.label_alt_id");
    let c_x = columns.get("_atom_site.Cartn_x");
    let c_y = columns.get("_atom_site.Cartn_y");
    let c_z = columns.get("_atom_site.Cartn_z");
    let c_occ = columns.get("_atom_site.occupancy");
    let c_b = columns.get("_atom_site.B_iso_or_equiv");
    let c_id = columns.get("_atom_site.id");
    let c_model = columns.get("_atom_site.pdbx_PDB_model_num");

    if c_x.is_none() || c_y.is_none() || c_z.is_none() {
        return Err("Missing coordinate columns in _atom_site loop".into());
    }

    // Raw model -> remapped contiguous state id (1..N in first-seen order).
    let mut model_to_state: HashMap<i32, usize> = HashMap::new();

    // 3. Read rows.
    loop {
        let peek = match tokens.peek() {
            Some(p) => p,
            None => break,
        };
        if peek.starts_with('_') || peek == "loop_" || peek.starts_with("data_") {
            break;
        }

        let mut row: Vec<&str> = Vec::with_capacity(col_idx);
        for _ in 0..col_idx {
            if let Some(val) = tokens.next() {
                row.push(val);
            } else {
                break;
            }
        }

        if row.len() != col_idx {
            break;
        }

        let group = c_group.map(|&i| row[i]).unwrap_or("ATOM");
        if group != "ATOM" && group != "HETATM" {
            continue;
        }
        let is_hetatm = group == "HETATM";

        let raw_model: i32 = c_model
            .and_then(|&i| row[i].parse().ok())
            .unwrap_or(1);
        let state_id = if let Some(&existing) = model_to_state.get(&raw_model) {
            existing
        } else {
            let next = model_to_state.len() + 1;
            model_to_state.insert(raw_model, next);
            next
        };

        while states.len() < state_id {
            states.push(Vec::new());
        }

        let x: f32 = c_x.and_then(|&i| row[i].parse().ok()).unwrap_or(0.0);
        let y: f32 = c_y.and_then(|&i| row[i].parse().ok()).unwrap_or(0.0);
        let z: f32 = c_z.and_then(|&i| row[i].parse().ok()).unwrap_or(0.0);

        let atom_name = c_atom.map(|&i| row[i]).unwrap_or("CA").to_string();
        let elem_sym = c_symbol.map(|&i| row[i]).unwrap_or("C").to_string();
        let resn = c_resn.map(|&i| row[i]).unwrap_or("UNK").to_string();
        let chain = c_chain
            .map(|&i| row[i].chars().next().unwrap_or(' '))
            .unwrap_or(' ');
        let resi = parse_cif_int(c_resi.map(|&i| row[i]).unwrap_or("0"));
        let ins_code = parse_cif_code(c_ins.map(|&i| row[i]).unwrap_or("?"), '\0');
        let alt = parse_cif_code(c_alt.map(|&i| row[i]).unwrap_or("?"), ' ');
        let occupancy: f32 = c_occ.and_then(|&i| row[i].parse().ok()).unwrap_or(1.0);
        let b_factor: f32 = c_b.and_then(|&i| row[i].parse().ok()).unwrap_or(0.0);
        let serial: u32 = c_id.and_then(|&i| row[i].parse().ok()).unwrap_or(0);

        let elem_data = element_by_symbol(&elem_sym);
        let atomic_number = elem_data
            .map(|e| ELEMENTS.iter().position(|x| std::ptr::eq(x, e)).unwrap_or(0) as u8)
            .unwrap_or(0);
        let vdw = elem_data.map(|e| e.vdw).unwrap_or(1.7);
        let color = elem_data.map(|e| e.color).unwrap_or([0.5, 0.5, 0.5]);

        let atom = AtomInfo {
            name: atom_name.clone(),
            element: atomic_number,
            elem_symbol: elem_sym,
            resn: resn.clone(),
            resi,
            ins_code,
            chain,
            alt,
            b_factor,
            occupancy,
            vdw,
            color,
            is_hetatm,
            serial,
            ..Default::default()
        };

        let identity = AtomIdentity {
            is_hetatm,
            chain,
            resn,
            resi,
            ins_code,
            name: atom_name,
            alt,
        };

        states[state_id - 1].push(ParsedAtom {
            identity,
            atom,
            coord: [x, y, z],
        });
    }

    Ok(())
}

fn parse_cif_int(s: &str) -> i32 {
    if s == "." || s == "?" {
        0
    } else {
        s.parse().unwrap_or(0)
    }
}

fn parse_cif_code(s: &str, default: char) -> char {
    if s == "." || s == "?" || s.is_empty() {
        default
    } else {
        s.chars().next().unwrap_or(default)
    }
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

#[cfg(test)]
mod tests {
    use std::path::Path;

    use super::parse_cif_string;

    fn cif_with_rows(rows: &[&str]) -> String {
        let mut s = String::from(
            "data_test\n\
             loop_\n\
             _atom_site.group_PDB\n\
             _atom_site.id\n\
             _atom_site.type_symbol\n\
             _atom_site.label_atom_id\n\
             _atom_site.label_comp_id\n\
             _atom_site.label_asym_id\n\
             _atom_site.label_entity_id\n\
             _atom_site.label_seq_id\n\
             _atom_site.Cartn_x\n\
             _atom_site.Cartn_y\n\
             _atom_site.Cartn_z\n\
             _atom_site.occupancy\n\
             _atom_site.B_iso_or_equiv\n\
             _atom_site.pdbx_PDB_model_num\n",
        );
        for row in rows {
            s.push_str(row);
            s.push('\n');
        }
        s
    }

    #[test]
    fn cif_multimodel_happy_path() {
        let cif = cif_with_rows(&[
            "ATOM 1 C CA ALA A 1 1 0.0 0.0 0.0 1.0 10.0 1",
            "ATOM 2 O O  ALA A 1 1 1.0 0.0 0.0 1.0 10.0 1",
            "ATOM 3 C CA ALA A 1 1 0.0 1.0 0.0 1.0 10.0 2",
            "ATOM 4 O O  ALA A 1 1 1.0 1.0 0.0 1.0 10.0 2",
            "ATOM 5 C CA ALA A 1 1 0.0 2.0 0.0 1.0 10.0 3",
            "ATOM 6 O O  ALA A 1 1 1.0 2.0 0.0 1.0 10.0 3",
        ]);

        let result = parse_cif_string(&cif, Path::new("multi.cif")).unwrap();
        assert_eq!(result.source_model_count, 3);
        assert_eq!(result.molecule.state_count(), 3);
        assert_eq!(result.molecule.atoms.len(), 2);
        assert!(result.warnings.is_empty());
    }

    #[test]
    fn cif_non_contiguous_model_ids_remap() {
        let cif = cif_with_rows(&[
            "ATOM 1 C CA ALA A 1 1 0.0 0.0 0.0 1.0 10.0 1",
            "ATOM 2 O O  ALA A 1 1 1.0 0.0 0.0 1.0 10.0 1",
            "ATOM 3 C CA ALA A 1 1 0.0 1.0 0.0 1.0 10.0 7",
            "ATOM 4 O O  ALA A 1 1 1.0 1.0 0.0 1.0 10.0 7",
        ]);

        let result = parse_cif_string(&cif, Path::new("noncontig.cif")).unwrap();
        assert_eq!(result.source_model_count, 2);
        assert_eq!(result.molecule.state_count(), 2);
        assert_eq!(result.molecule.coords_for_state(2)[0], [0.0, 1.0, 0.0]);
    }

    #[test]
    fn cif_mismatch_falls_back_to_model_1() {
        let cif = cif_with_rows(&[
            "ATOM 1 C CA ALA A 1 1 0.0 0.0 0.0 1.0 10.0 1",
            "ATOM 2 O O  ALA A 1 1 1.0 0.0 0.0 1.0 10.0 1",
            "ATOM 3 C XX ALA A 1 1 0.0 1.0 0.0 1.0 10.0 2",
            "ATOM 4 O O  ALA A 1 1 1.0 1.0 0.0 1.0 10.0 2",
        ]);

        let result = parse_cif_string(&cif, Path::new("mismatch.cif")).unwrap();
        assert_eq!(result.source_model_count, 2);
        assert_eq!(result.molecule.state_count(), 1);
        assert!(!result.warnings.is_empty());
    }
}
