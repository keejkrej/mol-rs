use std::collections::HashMap;
use std::path::Path;

use crate::core::atom::AtomInfo;
use crate::core::bond::BondInfo;
use crate::core::element::{element_by_symbol, ELEMENTS};
use crate::core::molecule::Molecule;
use crate::core::secondary_structure::SSType;
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
    serial: u32,
}

/// Parse a PDB file and return a Molecule + loader metadata.
pub fn load_pdb(path: &Path) -> Result<LoadResult, String> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| format!("Failed to read {}: {}", path.display(), e))?;
    parse_pdb_string(&content, path)
}

/// Parse PDB content from a string.
pub fn parse_pdb_string(content: &str, path: &Path) -> Result<LoadResult, String> {
    let name = path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown")
        .to_string();

    let mut mol = Molecule::new(name);
    let mut warnings = Vec::new();

    // Raw MODEL number -> remapped contiguous state index (1..N).
    let mut model_to_state: HashMap<i32, usize> = HashMap::new();
    let mut state_rows: Vec<Vec<ParsedAtom>> = Vec::new();
    let mut current_model_raw = 1;

    // Keep these from the first topology.
    let mut conect_lines: Vec<String> = Vec::new();
    let mut ss_assignments: Vec<(char, i32, i32, SSType)> = Vec::new();

    for line in content.lines() {
        let rec = if line.len() >= 6 { &line[..6] } else { line };

        match rec.trim() {
            "MODEL" => {
                if let Some(num) = parse_model_num(line) {
                    current_model_raw = num.max(1);
                }
            }
            "ENDMDL" => {}
            "ATOM" | "HETATM" => {
                if line.len() < 54 {
                    continue;
                }
                if let Some(parsed) = parse_atom_line(line, rec.starts_with("HETATM")) {
                    let state_id = if let Some(&existing) = model_to_state.get(&current_model_raw) {
                        existing
                    } else {
                        let next = model_to_state.len() + 1;
                        model_to_state.insert(current_model_raw, next);
                        next
                    };
                    while state_rows.len() < state_id {
                        state_rows.push(Vec::new());
                    }
                    state_rows[state_id - 1].push(parsed);
                }
            }
            "CONECT" => {
                conect_lines.push(line.to_string());
            }
            "HELIX" => {
                if line.len() >= 38 {
                    if let Some(ss) = parse_helix_line(line) {
                        ss_assignments.push(ss);
                    }
                }
            }
            "SHEET" => {
                if line.len() >= 38 {
                    if let Some(ss) = parse_sheet_line(line) {
                        ss_assignments.push(ss);
                    }
                }
            }
            _ => {}
        }
    }

    let source_model_count = state_rows.len();

    // Build molecule topology from state/model 1 and validate others.
    if let Some(first_state) = state_rows.first() {
        let mut serial_to_idx: HashMap<u32, usize> = HashMap::new();
        mol.atoms = Vec::with_capacity(first_state.len());
        for row in first_state {
            let idx = mol.atoms.len();
            serial_to_idx.insert(row.serial, idx);
            mol.atoms.push(row.atom.clone());
        }
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

        for line in &conect_lines {
            parse_conect_line(line, &serial_to_idx, &mut mol.bonds);
        }
    } else {
        mol.coord_sets = vec![Vec::new()];
    }

    // If no CONECT records, infer bonds from distances (state 1).
    if mol.bonds.is_empty() {
        infer_bonds(&mut mol);
    }

    // Apply secondary structure
    if !ss_assignments.is_empty() {
        mol.apply_ss(&ss_assignments);
    }

    // Build residue ranges and defaults
    mol.build_residues();
    mol.apply_default_representation();

    log::info!(
        "Loaded PDB: {} atoms, {} bonds, {} residues, {} state(s) from {} model(s)",
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

fn parse_model_num(line: &str) -> Option<i32> {
    if let Some(num_str) = line.get(10..14) {
        if let Ok(num) = num_str.trim().parse::<i32>() {
            return Some(num);
        }
    }
    line.split_whitespace().nth(1)?.parse().ok()
}

/// Parse a single ATOM/HETATM line.
fn parse_atom_line(line: &str, is_hetatm: bool) -> Option<ParsedAtom> {
    // PDB format columns (1-indexed):
    //  7-11  serial
    // 13-16  atom name
    // 17     alt loc
    // 18-20  resName
    // 22     chainID
    // 23-26  resSeq
    // 27     iCode
    // 31-38  x
    // 39-46  y
    // 47-54  z
    // 55-60  occupancy
    // 61-66  bfactor
    // 77-78  element symbol
    let bytes = line.as_bytes();
    let get = |start: usize, end: usize| -> &str {
        if end <= bytes.len() {
            std::str::from_utf8(&bytes[start..end]).unwrap_or("").trim()
        } else {
            ""
        }
    };

    let serial: u32 = get(6, 11).parse().ok()?;
    let atom_name = get(12, 16).to_string();
    let alt_raw = bytes.get(16).map(|&b| b as char).unwrap_or(' ');
    let alt = if alt_raw == ' ' { ' ' } else { alt_raw };
    let resn = get(17, 20).to_string();
    let chain = bytes.get(21).map(|&b| b as char).unwrap_or(' ');
    let resi: i32 = get(22, 26).parse().ok()?;
    let ins_raw = bytes.get(26).map(|&b| b as char).unwrap_or(' ');
    let ins_code = if ins_raw == ' ' { '\0' } else { ins_raw };

    let x: f32 = get(30, 38).parse().ok()?;
    let y: f32 = get(38, 46).parse().ok()?;
    let z: f32 = get(46, 54).parse().ok()?;

    let occupancy: f32 = get(54, 60).parse().unwrap_or(1.0);
    let b_factor: f32 = get(60, 66).parse().unwrap_or(0.0);

    // Element symbol: columns 77-78, or fall back to first non-digit char of atom name.
    let elem_sym = {
        let raw = get(76, 78);
        if raw.is_empty() {
            atom_name
                .chars()
                .find(|c| c.is_ascii_alphabetic())
                .map(|c| c.to_string())
                .unwrap_or_default()
        } else {
            raw.to_string()
        }
    };

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

    Some(ParsedAtom {
        identity,
        atom,
        coord: [x, y, z],
        serial,
    })
}

/// Parse CONECT record and add bonds (dedup later).
fn parse_conect_line(
    line: &str,
    serial_map: &HashMap<u32, usize>,
    bonds: &mut Vec<BondInfo>,
) {
    // CONECT columns: 7-11 = origin, then 12-16, 17-21, 22-26, 27-31 = bonded atoms
    let bytes = line.as_bytes();
    let get = |s: usize, e: usize| -> Option<u32> {
        if e <= bytes.len() {
            std::str::from_utf8(&bytes[s..e])
                .ok()?
                .trim()
                .parse()
                .ok()
        } else {
            None
        }
    };

    let origin_serial = match get(6, 11) {
        Some(s) => s,
        None => return,
    };
    let origin_idx = match serial_map.get(&origin_serial) {
        Some(&i) => i,
        None => return,
    };

    for &(s, e) in &[(11, 16), (16, 21), (21, 26), (26, 31)] {
        if let Some(partner_serial) = get(s, e) {
            if let Some(&partner_idx) = serial_map.get(&partner_serial) {
                // Only add bond in one direction (lower index first)
                if origin_idx < partner_idx {
                    let exists = bonds
                        .iter()
                        .any(|b| b.atom_a == origin_idx && b.atom_b == partner_idx);
                    if !exists {
                        bonds.push(BondInfo {
                            atom_a: origin_idx,
                            atom_b: partner_idx,
                            order: 1,
                        });
                    }
                }
            }
        }
    }
}

/// Parse HELIX record.
fn parse_helix_line(line: &str) -> Option<(char, i32, i32, SSType)> {
    let bytes = line.as_bytes();
    let chain = *bytes.get(19)? as char;
    let start: i32 = std::str::from_utf8(bytes.get(21..25)?)
        .ok()?
        .trim()
        .parse()
        .ok()?;
    let end: i32 = std::str::from_utf8(bytes.get(33..37)?)
        .ok()?
        .trim()
        .parse()
        .ok()?;
    Some((chain, start, end, SSType::Helix))
}

/// Parse SHEET record.
fn parse_sheet_line(line: &str) -> Option<(char, i32, i32, SSType)> {
    let bytes = line.as_bytes();
    let chain = *bytes.get(21)? as char;
    let start: i32 = std::str::from_utf8(bytes.get(22..26)?)
        .ok()?
        .trim()
        .parse()
        .ok()?;
    let end: i32 = std::str::from_utf8(bytes.get(33..37)?)
        .ok()?
        .trim()
        .parse()
        .ok()?;
    Some((chain, start, end, SSType::Sheet))
}

/// Infer covalent bonds based on inter-atomic distance.
/// Uses a simple distance cutoff: bonded if dist < (vdw_a + vdw_b) * 0.6
/// For efficiency, uses a spatial grid.
pub fn infer_bonds(mol: &mut Molecule) {
    let n = mol.atoms.len();
    if n == 0 {
        return;
    }

    let coords = mol.coords_for_state(1).to_vec();
    if coords.len() != n {
        return;
    }

    // Simple O(n*n) for small molecules, grid for larger
    let cutoff_scale = 0.6f32;

    if n <= 5000 {
        // Brute force for small structures
        for i in 0..n {
            for j in (i + 1)..n {
                let ai = &mol.atoms[i];
                let aj = &mol.atoms[j];
                let pi = coords[i];
                let pj = coords[j];

                let dx = pi[0] - pj[0];
                let dy = pi[1] - pj[1];
                let dz = pi[2] - pj[2];
                let dist_sq = dx * dx + dy * dy + dz * dz;

                let max_dist = (ai.vdw + aj.vdw) * cutoff_scale;
                let max_dist_sq = max_dist * max_dist;

                // Also enforce a minimum distance to avoid clashes
                if dist_sq < max_dist_sq && dist_sq > 0.16 {
                    mol.bonds.push(BondInfo {
                        atom_a: i,
                        atom_b: j,
                        order: 1,
                    });
                }
            }
        }
    } else {
        // Grid-based approach for large structures
        let cell_size = 2.5f32; // Å, covers most covalent bonds

        // Find bounding box
        let mut min = [f32::MAX; 3];
        let mut max = [f32::MIN; 3];
        for p in &coords {
            for k in 0..3 {
                min[k] = min[k].min(p[k]);
                max[k] = max[k].max(p[k]);
            }
        }

        let dims = [
            ((max[0] - min[0]) / cell_size) as usize + 1,
            ((max[1] - min[1]) / cell_size) as usize + 1,
            ((max[2] - min[2]) / cell_size) as usize + 1,
        ];

        let mut grid: HashMap<(usize, usize, usize), Vec<usize>> = HashMap::new();

        for (i, p) in coords.iter().enumerate() {
            let cx = ((p[0] - min[0]) / cell_size) as usize;
            let cy = ((p[1] - min[1]) / cell_size) as usize;
            let cz = ((p[2] - min[2]) / cell_size) as usize;
            grid.entry((cx, cy, cz)).or_default().push(i);
        }

        for (&(cx, cy, cz), atoms_in_cell) in &grid {
            // Check this cell and all 26 neighbors
            let x_range = cx.saturating_sub(1)..=(cx + 1).min(dims[0] - 1);
            for nx in x_range {
                let y_range = cy.saturating_sub(1)..=(cy + 1).min(dims[1] - 1);
                for ny in y_range.clone() {
                    let z_range = cz.saturating_sub(1)..=(cz + 1).min(dims[2] - 1);
                    for nz in z_range.clone() {
                        if let Some(neighbor_atoms) = grid.get(&(nx, ny, nz)) {
                            for &i in atoms_in_cell {
                                for &j in neighbor_atoms {
                                    if i >= j {
                                        continue;
                                    }
                                    let ai = &mol.atoms[i];
                                    let aj = &mol.atoms[j];
                                    let pi = coords[i];
                                    let pj = coords[j];

                                    let dx = pi[0] - pj[0];
                                    let dy = pi[1] - pj[1];
                                    let dz = pi[2] - pj[2];
                                    let dist_sq = dx * dx + dy * dy + dz * dz;

                                    let max_dist = (ai.vdw + aj.vdw) * cutoff_scale;
                                    if dist_sq < max_dist * max_dist && dist_sq > 0.16 {
                                        mol.bonds.push(BondInfo {
                                            atom_a: i,
                                            atom_b: j,
                                            order: 1,
                                        });
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use super::parse_pdb_string;

    fn atom_line(
        rec: &str,
        serial: u32,
        name: &str,
        resn: &str,
        chain: char,
        resi: i32,
        x: f32,
        y: f32,
        z: f32,
        elem: &str,
    ) -> String {
        format!(
            "{:<6}{:>5} {:<4}{:1}{:>3} {:1}{:>4}{:1}   {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}",
            rec, serial, name, "", resn, chain, resi, "", x, y, z, 1.0, 10.0, elem
        )
    }

    #[test]
    fn pdb_multimodel_happy_path() {
        let pdb = format!(
            "MODEL        1\n{}\n{}\nENDMDL\nMODEL        2\n{}\n{}\nENDMDL\n",
            atom_line("ATOM", 1, "CA", "ALA", 'A', 1, 0.0, 0.0, 0.0, "C"),
            atom_line("ATOM", 2, "O", "ALA", 'A', 1, 1.0, 0.0, 0.0, "O"),
            atom_line("ATOM", 1, "CA", "ALA", 'A', 1, 0.0, 1.0, 0.0, "C"),
            atom_line("ATOM", 2, "O", "ALA", 'A', 1, 1.0, 1.0, 0.0, "O")
        );

        let result = parse_pdb_string(&pdb, Path::new("multi.pdb")).unwrap();
        assert_eq!(result.source_model_count, 2);
        assert_eq!(result.molecule.state_count(), 2);
        assert_eq!(result.molecule.atoms.len(), 2);
        assert!(result.warnings.is_empty());
    }

    #[test]
    fn pdb_non_contiguous_model_ids_remap() {
        let pdb = format!(
            "MODEL        1\n{}\n{}\nENDMDL\nMODEL        7\n{}\n{}\nENDMDL\n",
            atom_line("ATOM", 1, "CA", "ALA", 'A', 1, 0.0, 0.0, 0.0, "C"),
            atom_line("ATOM", 2, "O", "ALA", 'A', 1, 1.0, 0.0, 0.0, "O"),
            atom_line("ATOM", 1, "CA", "ALA", 'A', 1, 0.0, 2.0, 0.0, "C"),
            atom_line("ATOM", 2, "O", "ALA", 'A', 1, 1.0, 2.0, 0.0, "O")
        );

        let result = parse_pdb_string(&pdb, Path::new("noncontig.pdb")).unwrap();
        assert_eq!(result.source_model_count, 2);
        assert_eq!(result.molecule.state_count(), 2);
        assert_eq!(result.molecule.coords_for_state(2)[0], [0.0, 2.0, 0.0]);
    }

    #[test]
    fn pdb_mismatch_falls_back_to_model_1() {
        let pdb = format!(
            "MODEL        1\n{}\n{}\nENDMDL\nMODEL        2\n{}\n{}\nENDMDL\n",
            atom_line("ATOM", 1, "CA", "ALA", 'A', 1, 0.0, 0.0, 0.0, "C"),
            atom_line("ATOM", 2, "O", "ALA", 'A', 1, 1.0, 0.0, 0.0, "O"),
            atom_line("ATOM", 1, "CB", "ALA", 'A', 1, 0.0, 1.0, 0.0, "C"),
            atom_line("ATOM", 2, "O", "ALA", 'A', 1, 1.0, 1.0, 0.0, "O")
        );

        let result = parse_pdb_string(&pdb, Path::new("mismatch.pdb")).unwrap();
        assert_eq!(result.source_model_count, 2);
        assert_eq!(result.molecule.state_count(), 1);
        assert!(!result.warnings.is_empty());
    }
}
