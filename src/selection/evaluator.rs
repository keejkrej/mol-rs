use crate::core::atom::{REP_ALL, REP_CARTOON, REP_LINES, REP_SPHERES, REP_STICKS};
use crate::core::element::{element_by_number, is_metal_atomic_number};
use crate::core::molecule::Molecule;
use crate::core::residue::{is_nucleic, is_protein, is_solvent};
use crate::core::secondary_structure::SSType;
use crate::selection::parser::{AtomProperty, CompareOp, Selector};

/// Evaluate a selection expression against a molecule.
/// Returns a Vec<bool> with one entry per atom — true means selected.
pub fn evaluate(sel: &Selector, mol: &Molecule) -> Vec<bool> {
    evaluate_with_coords(
        sel,
        mol,
        mol.coord_sets
            .first()
            .map_or(&[][..], |coords| coords.as_slice()),
    )
}

/// Evaluate a selection expression against a molecule using the provided coordinates.
/// Returns a Vec<bool> with one entry per atom — true means selected.
pub fn evaluate_with_coords(sel: &Selector, mol: &Molecule, coords: &[[f32; 3]]) -> Vec<bool> {
    let n = mol.atoms.len();
    match sel {
        Selector::All => vec![true; n],
        Selector::None => vec![false; n],
        Selector::Enabled => vec![mol.visible; n],
        Selector::Visible => mol
            .atoms
            .iter()
            .map(|a| mol.visible && a.vis_rep != 0)
            .collect(),
        Selector::Present => (0..n).map(|idx| idx < coords.len()).collect(),
        Selector::Bonded => bonded_atom_mask(mol),
        Selector::Chain(ch) => mol.atoms.iter().map(|a| a.chain == *ch).collect(),
        Selector::ChainPattern(pattern) => mol
            .atoms
            .iter()
            .map(|a| matches_alpha_list(&a.chain.to_string(), pattern))
            .collect(),
        Selector::Resi(lo, hi) => mol
            .atoms
            .iter()
            .map(|a| a.resi >= *lo && a.resi <= *hi)
            .collect(),
        Selector::ResiList(ranges) => mol
            .atoms
            .iter()
            .map(|a| in_numeric_ranges(a.resi, ranges))
            .collect(),
        Selector::Name(name) => mol
            .atoms
            .iter()
            .map(|a| matches_alpha_list(a.name.trim(), name))
            .collect(),
        Selector::Resn(resn) => mol
            .atoms
            .iter()
            .map(|a| matches_alpha_list(a.resn.trim(), resn))
            .collect(),
        Selector::Organic => mol.atoms.iter().map(|a| !a.is_hetatm).collect(),
        Selector::Inorganic => mol.atoms.iter().map(|a| a.is_hetatm).collect(),
        Selector::Serial(lo, hi) => mol
            .atoms
            .iter()
            .map(|a| {
                let serial = i32::try_from(a.serial).unwrap_or(0);
                serial >= *lo && serial <= *hi
            })
            .collect(),
        Selector::SerialList(ranges) => mol
            .atoms
            .iter()
            .map(|a| {
                let serial = i32::try_from(a.serial).unwrap_or(0);
                in_numeric_ranges(serial, ranges)
            })
            .collect(),
        Selector::Model(pattern) => {
            let selected = object_name_matches(pattern, &mol.name);
            vec![selected; n]
        }
        Selector::Rep(rep) => {
            let mask = rep_mask(rep);
            mol.atoms
                .iter()
                .map(|atom| mask != 0 && atom.vis_rep & mask != 0)
                .collect()
        }
        Selector::Color(color) => {
            let rgb = color_rgb(color);
            mol.atoms
                .iter()
                .map(|atom| rgb.is_some_and(|rgb| colors_match(atom.color, rgb)))
                .collect()
        }
        Selector::Property(property, op, value) => mol
            .atoms
            .iter()
            .enumerate()
            .map(|(idx, atom)| {
                property_value(*property, atom, idx, coords)
                    .is_some_and(|atom_value| compare_float(atom_value, *op, *value))
            })
            .collect(),
        Selector::Byres(inner) => {
            let inner_mask = evaluate_with_coords(inner, mol, coords);
            let mut residues: Vec<(char, i32, char)> = Vec::new();

            for (idx, selected) in inner_mask.iter().enumerate() {
                if *selected {
                    if let Some(atom) = mol.atoms.get(idx) {
                        let key = (atom.chain, atom.resi, atom.ins_code);
                        if !residues.contains(&key) {
                            residues.push(key);
                        }
                    }
                }
            }

            mol.atoms
                .iter()
                .map(|atom| residues.contains(&(atom.chain, atom.resi, atom.ins_code)))
                .collect()
        }
        Selector::Bychain(inner) => {
            let inner_mask = evaluate_with_coords(inner, mol, coords);
            let mut chains: Vec<char> = Vec::new();

            for (idx, selected) in inner_mask.iter().enumerate() {
                if *selected {
                    if let Some(atom) = mol.atoms.get(idx) {
                        if !chains.contains(&atom.chain) {
                            chains.push(atom.chain);
                        }
                    }
                }
            }

            mol.atoms
                .iter()
                .map(|atom| chains.contains(&atom.chain))
                .collect()
        }
        Selector::Byobject(inner) => {
            let inner_mask = evaluate_with_coords(inner, mol, coords);
            vec![inner_mask.iter().any(|selected| *selected); n]
        }
        Selector::Bymolecule(inner) => {
            let inner_mask = evaluate_with_coords(inner, mol, coords);
            bonded_component_mask(&inner_mask, mol)
        }
        Selector::First(inner) => {
            let inner_mask = evaluate_with_coords(inner, mol, coords);
            let mut out = vec![false; n];
            if let Some(idx) = inner_mask.iter().position(|selected| *selected) {
                out[idx] = true;
            }
            out
        }
        Selector::Last(inner) => {
            let inner_mask = evaluate_with_coords(inner, mol, coords);
            let mut out = vec![false; n];
            if let Some(idx) = inner_mask.iter().rposition(|selected| *selected) {
                out[idx] = true;
            }
            out
        }
        Selector::Neighbor(inner) => {
            let inner_mask = evaluate_with_coords(inner, mol, coords);
            let mut out = vec![false; n];

            if inner_mask.iter().any(|&x| x) {
                for bond in &mol.bonds {
                    if bond.atom_a < n && bond.atom_b < n {
                        if inner_mask[bond.atom_a] && !inner_mask[bond.atom_b] {
                            out[bond.atom_b] = true;
                        }
                        if inner_mask[bond.atom_b] && !inner_mask[bond.atom_a] {
                            out[bond.atom_a] = true;
                        }
                    }
                }
            }

            out
        }
        Selector::Around(distance, inner) => {
            let inner_mask = evaluate_with_coords(inner, mol, coords);
            distance_mask(*distance, inner, mol, coords)
                .into_iter()
                .enumerate()
                .map(|(idx, selected)| selected && !inner_mask.get(idx).copied().unwrap_or(false))
                .collect()
        }
        Selector::Within(distance, inner) => distance_mask(*distance, inner, mol, coords),
        Selector::Expand(distance, inner) => distance_mask(*distance, inner, mol, coords),
        Selector::Extend(count, inner) => {
            let inner_mask = evaluate_with_coords(inner, mol, coords);
            extend_bond_mask(&inner_mask, mol, *count)
        }
        Selector::Beyond(distance, inner) => distance_mask(*distance, inner, mol, coords)
            .into_iter()
            .enumerate()
            .map(|(idx, selected)| idx < coords.len() && !selected)
            .collect(),
        Selector::NearTo(distance, inner) => {
            let inner_mask = evaluate_with_coords(inner, mol, coords);
            distance_mask(*distance, inner, mol, coords)
                .into_iter()
                .enumerate()
                .map(|(idx, selected)| selected && !inner_mask.get(idx).copied().unwrap_or(false))
                .collect()
        }
        Selector::Elem(sym) => mol
            .atoms
            .iter()
            .map(|a| {
                let ed = element_by_number(a.element);
                matches_alpha_list(ed.symbol, sym)
            })
            .collect(),
        Selector::Alt(alt) => mol.atoms.iter().map(|a| a.alt == *alt).collect(),
        Selector::AltPattern(pattern) => mol
            .atoms
            .iter()
            .map(|a| matches_alpha_list(&a.alt.to_string(), pattern))
            .collect(),
        Selector::SS(ss) => mol
            .atoms
            .iter()
            .map(|a| matches_secondary_structure(a.ss_type, ss))
            .collect(),
        Selector::Hetatm => mol.atoms.iter().map(|a| a.is_hetatm).collect(),
        Selector::Hydrogen => mol
            .atoms
            .iter()
            .map(|a| {
                a.element == 1
                    || a.elem_symbol.eq_ignore_ascii_case("H")
                    || a.name.trim_start_matches(char::is_numeric).starts_with('H')
            })
            .collect(),
        Selector::Solvent => mol.atoms.iter().map(|a| is_solvent(&a.resn)).collect(),
        Selector::Polymer => mol
            .atoms
            .iter()
            .map(|a| is_protein(&a.resn) || is_nucleic(&a.resn))
            .collect(),
        Selector::Protein => mol.atoms.iter().map(|a| is_protein(&a.resn)).collect(),
        Selector::Nucleic => mol.atoms.iter().map(|a| is_nucleic(&a.resn)).collect(),
        Selector::Metals => mol
            .atoms
            .iter()
            .map(|a| is_metal_atomic_number(a.element))
            .collect(),
        Selector::Guide => mol.atoms.iter().map(is_guide_atom).collect(),
        Selector::And(left, right) => {
            let l = evaluate_with_coords(left, mol, coords);
            let r = evaluate_with_coords(right, mol, coords);
            l.iter().zip(r.iter()).map(|(a, b)| *a && *b).collect()
        }
        Selector::Or(left, right) => {
            let l = evaluate_with_coords(left, mol, coords);
            let r = evaluate_with_coords(right, mol, coords);
            l.iter().zip(r.iter()).map(|(a, b)| *a || *b).collect()
        }
        Selector::Not(inner) => {
            let v = evaluate_with_coords(inner, mol, coords);
            v.iter().map(|x| !x).collect()
        }
    }
}

/// Count how many atoms are selected.
pub fn count_selected(mask: &[bool]) -> usize {
    mask.iter().filter(|&&b| b).count()
}

fn matches_alpha_list(value: &str, pattern: &str) -> bool {
    pattern.split(['+', ',']).any(|item| {
        let item = item.trim();
        !item.is_empty() && wildcard_match_ci(item, value)
    })
}

fn in_numeric_ranges(value: i32, ranges: &[(i32, i32)]) -> bool {
    ranges.iter().any(|(lo, hi)| value >= *lo && value <= *hi)
}

fn object_name_matches(pattern: &str, name: &str) -> bool {
    let pattern = pattern.trim();
    if pattern.eq_ignore_ascii_case("all") || pattern == "*" {
        return true;
    }

    wildcard_match_ci(pattern, name)
}

fn wildcard_match_ci(pattern: &str, text: &str) -> bool {
    let pattern = pattern.to_ascii_lowercase();
    let text = text.to_ascii_lowercase();
    let pattern = pattern.as_bytes();
    let text = text.as_bytes();

    let (mut p, mut t) = (0usize, 0usize);
    let mut star = None;
    let mut star_text = 0usize;

    while t < text.len() {
        if p < pattern.len() && pattern[p] == text[t] {
            p += 1;
            t += 1;
        } else if p < pattern.len() && pattern[p] == b'*' {
            star = Some(p);
            p += 1;
            star_text = t;
        } else if let Some(star_idx) = star {
            p = star_idx + 1;
            star_text += 1;
            t = star_text;
        } else {
            return false;
        }
    }

    while p < pattern.len() && pattern[p] == b'*' {
        p += 1;
    }

    p == pattern.len()
}

fn property_value(
    property: AtomProperty,
    atom: &crate::core::atom::AtomInfo,
    idx: usize,
    coords: &[[f32; 3]],
) -> Option<f32> {
    match property {
        AtomProperty::BFactor => Some(atom.b_factor),
        AtomProperty::Occupancy => Some(atom.occupancy),
        AtomProperty::X => coords.get(idx).map(|coord| coord[0]),
        AtomProperty::Y => coords.get(idx).map(|coord| coord[1]),
        AtomProperty::Z => coords.get(idx).map(|coord| coord[2]),
    }
}

fn rep_mask(rep: &str) -> u32 {
    rep.split(['+', ','])
        .filter_map(|item| match item.trim().to_ascii_lowercase().as_str() {
            "lines" | "line" | "wire" | "wires" => Some(REP_LINES),
            "sticks" | "stick" => Some(REP_STICKS),
            "spheres" | "sphere" => Some(REP_SPHERES),
            "cartoon" | "ribbon" => Some(REP_CARTOON),
            "everything" | "all" => Some(REP_ALL),
            _ => None,
        })
        .fold(0, |mask, flag| mask | flag)
}

fn color_rgb(name: &str) -> Option<[f32; 3]> {
    match name.to_ascii_lowercase().as_str() {
        "red" => Some([1.0, 0.2, 0.2]),
        "green" => Some([0.2, 1.0, 0.2]),
        "blue" => Some([0.2, 0.2, 1.0]),
        "yellow" => Some([1.0, 1.0, 0.2]),
        "cyan" => Some([0.2, 1.0, 1.0]),
        "magenta" => Some([1.0, 0.2, 1.0]),
        "orange" => Some([1.0, 0.6, 0.2]),
        "white" => Some([1.0, 1.0, 1.0]),
        "gray" | "grey" => Some([0.5, 0.5, 0.5]),
        "pink" => Some([1.0, 0.65, 0.85]),
        "salmon" => Some([1.0, 0.6, 0.5]),
        "purple" => Some([0.6, 0.2, 0.8]),
        _ => None,
    }
}

fn colors_match(left: [f32; 3], right: [f32; 3]) -> bool {
    left.iter()
        .zip(right.iter())
        .all(|(left, right)| (left - right).abs() < 0.0001)
}

fn matches_secondary_structure(ss_type: SSType, pattern: &str) -> bool {
    let values: &[&str] = match ss_type {
        SSType::Helix => &["H", "HELIX"],
        SSType::Sheet => &["S", "SHEET"],
        SSType::Loop => &["L", "LOOP", "C", "COIL"],
    };

    values
        .iter()
        .any(|value| matches_alpha_list(value, pattern))
}

fn compare_float(left: f32, op: CompareOp, right: f32) -> bool {
    match op {
        CompareOp::Greater => left > right,
        CompareOp::Less => left < right,
        CompareOp::Equal => (left - right).abs() < 0.0001,
        CompareOp::GreaterEqual => left >= right,
        CompareOp::LessEqual => left <= right,
    }
}

fn distance_mask(
    distance: f32,
    inner: &Selector,
    mol: &Molecule,
    coords: &[[f32; 3]],
) -> Vec<bool> {
    let n = mol.atoms.len();
    let inner_mask = evaluate_with_coords(inner, mol, coords);
    let threshold2 = distance.max(0.0) * distance.max(0.0);
    let mut center_coords: Vec<[f32; 3]> = Vec::new();

    for (idx, selected) in inner_mask.iter().enumerate() {
        if *selected && idx < coords.len() {
            center_coords.push(coords[idx]);
        }
    }

    if center_coords.is_empty() {
        return vec![false; n];
    }

    mol.atoms
        .iter()
        .enumerate()
        .map(|(idx, _)| {
            if idx >= coords.len() {
                return false;
            }
            let p = coords[idx];
            center_coords.iter().any(|c| {
                let dx = p[0] - c[0];
                let dy = p[1] - c[1];
                let dz = p[2] - c[2];
                dx * dx + dy * dy + dz * dz <= threshold2
            })
        })
        .collect()
}

fn bonded_atom_mask(mol: &Molecule) -> Vec<bool> {
    let mut out = vec![false; mol.atoms.len()];

    for bond in mol.bonds.iter().filter(|bond| bond.order > 0) {
        if bond.atom_a < out.len() {
            out[bond.atom_a] = true;
        }
        if bond.atom_b < out.len() {
            out[bond.atom_b] = true;
        }
    }

    out
}

fn bonded_component_mask(inner_mask: &[bool], mol: &Molecule) -> Vec<bool> {
    let n = mol.atoms.len();
    let mut out = vec![false; n];
    let mut stack: Vec<usize> = inner_mask
        .iter()
        .enumerate()
        .filter_map(|(idx, selected)| (*selected && idx < n).then_some(idx))
        .collect();

    while let Some(idx) = stack.pop() {
        if out[idx] {
            continue;
        }
        out[idx] = true;

        for bond in mol.bonds.iter().filter(|bond| bond.order > 0) {
            let next = if bond.atom_a == idx {
                Some(bond.atom_b)
            } else if bond.atom_b == idx {
                Some(bond.atom_a)
            } else {
                None
            };

            if let Some(next) = next {
                if next < n && !out[next] {
                    stack.push(next);
                }
            }
        }
    }

    out
}

fn extend_bond_mask(inner_mask: &[bool], mol: &Molecule, count: usize) -> Vec<bool> {
    let n = mol.atoms.len();
    let mut out: Vec<bool> = inner_mask.iter().take(n).copied().collect();
    out.resize(n, false);

    for _ in 0..count {
        let previous = out.clone();
        for bond in mol.bonds.iter().filter(|bond| bond.order > 0) {
            if bond.atom_a < n && bond.atom_b < n {
                if previous[bond.atom_a] {
                    out[bond.atom_b] = true;
                }
                if previous[bond.atom_b] {
                    out[bond.atom_a] = true;
                }
            }
        }
    }

    out
}

fn is_guide_atom(atom: &crate::core::atom::AtomInfo) -> bool {
    let name = atom.name.trim();
    if is_protein(&atom.resn) {
        return name.eq_ignore_ascii_case("CA");
    }
    if is_nucleic(&atom.resn) {
        return matches!(
            name.to_ascii_uppercase().as_str(),
            "C4'" | "C4*" | "C3'" | "C3*"
        );
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::atom::{AtomInfo, REP_LINES};
    use crate::core::molecule::Molecule;

    fn test_molecule() -> Molecule {
        let mut mol = Molecule::new("test".into());
        mol.atoms.push(AtomInfo {
            name: "N".to_string(),
            element: 7,
            serial: 10,
            is_hetatm: false,
            resi: 1,
            resn: "ALA".to_string(),
            chain: 'A',
            ..AtomInfo::default()
        });
        mol.atoms.push(AtomInfo {
            name: "C".to_string(),
            element: 6,
            serial: 11,
            is_hetatm: true,
            alt: 'A',
            resi: 1,
            resn: "ATP".to_string(),
            chain: 'A',
            ss_type: SSType::Helix,
            ..AtomInfo::default()
        });
        mol.atoms.push(AtomInfo {
            name: "C".to_string(),
            element: 6,
            serial: 12,
            is_hetatm: false,
            resi: 2,
            resn: "GLY".to_string(),
            chain: 'B',
            ss_type: SSType::Sheet,
            ..AtomInfo::default()
        });
        mol.coord_sets = vec![
            vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [10.0, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [10.0, 0.0, 0.0]],
        ];
        mol
    }

    #[test]
    fn evaluate_organic_inorganic() {
        let mol = test_molecule();
        let organic = evaluate(&Selector::Organic, &mol);
        let inorganic = evaluate(&Selector::Inorganic, &mol);
        assert_eq!(organic, vec![true, false, true]);
        assert_eq!(inorganic, vec![false, true, false]);
    }

    #[test]
    fn evaluate_serial_range() {
        let mol = test_molecule();
        let serial = evaluate(&Selector::Serial(11, 12), &mol);
        assert_eq!(serial, vec![false, true, true]);
    }

    #[test]
    fn evaluate_chain_pattern_selector() {
        let mol = test_molecule();

        assert_eq!(
            evaluate(&Selector::ChainPattern("A+B".to_string()), &mol),
            vec![true, true, true]
        );
        assert_eq!(
            evaluate(&Selector::ChainPattern("A".to_string()), &mol),
            vec![true, true, false]
        );
        assert_eq!(
            evaluate(&Selector::ChainPattern("Z*".to_string()), &mol),
            vec![false, false, false]
        );
    }

    #[test]
    fn evaluate_plus_separated_numeric_lists() {
        let mol = test_molecule();

        assert_eq!(
            evaluate(&Selector::ResiList(vec![(2, 2), (4, 4)]), &mol),
            vec![false, false, true]
        );
        assert_eq!(
            evaluate(&Selector::SerialList(vec![(10, 10), (12, 12)]), &mol),
            vec![true, false, true]
        );
    }

    #[test]
    fn evaluate_model_selector() {
        let mut mol = test_molecule();
        mol.name = "Protein_A".to_string();

        assert_eq!(
            evaluate(&Selector::Model("protein_*".to_string()), &mol),
            vec![true, true, true]
        );
        assert_eq!(
            evaluate(&Selector::Model("ligand".to_string()), &mol),
            vec![false, false, false]
        );
    }

    #[test]
    fn evaluate_rep_selector() {
        let mut mol = test_molecule();
        mol.atoms[0].vis_rep = REP_LINES;
        mol.atoms[1].vis_rep = REP_STICKS;
        mol.atoms[2].vis_rep = REP_SPHERES | REP_CARTOON;

        assert_eq!(
            evaluate(&Selector::Rep("lines".to_string()), &mol),
            vec![true, false, false]
        );
        assert_eq!(
            evaluate(&Selector::Rep("wire".to_string()), &mol),
            vec![true, false, false]
        );
        assert_eq!(
            evaluate(&Selector::Rep("sticks+spheres".to_string()), &mol),
            vec![false, true, true]
        );
        assert_eq!(
            evaluate(&Selector::Rep("ribbon".to_string()), &mol),
            vec![false, false, true]
        );
        assert_eq!(
            evaluate(&Selector::Rep("everything".to_string()), &mol),
            vec![true, true, true]
        );
        assert_eq!(
            evaluate(&Selector::Rep("sticks,spheres".to_string()), &mol),
            vec![false, true, true]
        );
        assert_eq!(
            evaluate(&Selector::Rep("surface".to_string()), &mol),
            vec![false, false, false]
        );
    }

    #[test]
    fn evaluate_color_selector() {
        let mut mol = test_molecule();
        mol.atoms[0].color = [1.0, 0.2, 0.2];
        mol.atoms[1].color = [0.5, 0.5, 0.5];
        mol.atoms[2].color = [0.2, 0.2, 1.0];

        assert_eq!(
            evaluate(&Selector::Color("red".to_string()), &mol),
            vec![true, false, false]
        );
        assert_eq!(
            evaluate(&Selector::Color("grey".to_string()), &mol),
            vec![false, true, false]
        );
        assert_eq!(
            evaluate(&Selector::Color("unknown".to_string()), &mol),
            vec![false, false, false]
        );
    }

    #[test]
    fn evaluate_property_selectors() {
        let mut mol = test_molecule();
        mol.atoms[0].b_factor = 10.0;
        mol.atoms[1].b_factor = 20.0;
        mol.atoms[2].b_factor = 30.0;
        mol.atoms[0].occupancy = 0.25;
        mol.atoms[1].occupancy = 0.5;
        mol.atoms[2].occupancy = 1.0;

        assert_eq!(
            evaluate(
                &Selector::Property(AtomProperty::BFactor, CompareOp::Less, 25.0),
                &mol,
            ),
            vec![true, true, false]
        );
        assert_eq!(
            evaluate(
                &Selector::Property(AtomProperty::Occupancy, CompareOp::GreaterEqual, 0.5),
                &mol,
            ),
            vec![false, true, true]
        );
        assert_eq!(
            evaluate_with_coords(
                &Selector::Property(AtomProperty::X, CompareOp::LessEqual, 0.2),
                &mol,
                &mol.coord_sets[1],
            ),
            vec![true, true, false]
        );
        assert_eq!(
            evaluate_with_coords(
                &Selector::Property(AtomProperty::Z, CompareOp::Equal, 0.0),
                &mol,
                &mol.coord_sets[1],
            ),
            vec![true, true, true]
        );
    }

    #[test]
    fn evaluate_plus_separated_alpha_lists() {
        let mol = test_molecule();

        assert_eq!(
            evaluate(&Selector::Name("N+CA".to_string()), &mol),
            vec![true, false, false]
        );
        assert_eq!(
            evaluate(&Selector::Resn("ALA+GLY".to_string()), &mol),
            vec![true, false, true]
        );
        assert_eq!(
            evaluate(&Selector::Elem("C+N".to_string()), &mol),
            vec![true, true, true]
        );
        assert_eq!(
            evaluate(&Selector::Elem("O+S".to_string()), &mol),
            vec![false, false, false]
        );
    }

    #[test]
    fn evaluate_wildcard_alpha_lists() {
        let mol = test_molecule();

        assert_eq!(
            evaluate(&Selector::Name("C*".to_string()), &mol),
            vec![false, true, true]
        );
        assert_eq!(
            evaluate(&Selector::Resn("AL*+GL*".to_string()), &mol),
            vec![true, false, true]
        );
        assert_eq!(
            evaluate(&Selector::Elem("C*".to_string()), &mol),
            vec![false, true, true]
        );
    }

    #[test]
    fn evaluate_enabled_visible_selectors() {
        let mut mol = test_molecule();
        mol.atoms[0].vis_rep = REP_LINES;
        mol.atoms[1].vis_rep = 0;
        mol.atoms[2].vis_rep = REP_LINES;

        assert_eq!(evaluate(&Selector::Enabled, &mol), vec![true, true, true]);
        assert_eq!(evaluate(&Selector::Visible, &mol), vec![true, false, true]);

        mol.visible = false;
        assert_eq!(
            evaluate(&Selector::Enabled, &mol),
            vec![false, false, false]
        );
        assert_eq!(
            evaluate(&Selector::Visible, &mol),
            vec![false, false, false]
        );
    }

    #[test]
    fn evaluate_present_selector() {
        let mol = test_molecule();

        assert_eq!(
            evaluate_with_coords(&Selector::Present, &mol, &mol.coord_sets[0][..2]),
            vec![true, true, false]
        );
        assert_eq!(
            evaluate_with_coords(&Selector::Present, &mol, &[]),
            vec![false, false, false]
        );
    }

    #[test]
    fn evaluate_bonded_selector() {
        let mut mol = test_molecule();
        mol.bonds.push(crate::core::bond::BondInfo {
            atom_a: 0,
            atom_b: 1,
            order: 1,
        });

        assert_eq!(evaluate(&Selector::Bonded, &mol), vec![true, true, false]);
    }

    #[test]
    fn evaluate_hydrogen_and_solvent_selectors() {
        let mut mol = Molecule::new("test".into());
        mol.atoms.push(AtomInfo {
            name: "H1".to_string(),
            elem_symbol: "H".to_string(),
            element: 1,
            resn: "ALA".to_string(),
            ..AtomInfo::default()
        });
        mol.atoms.push(AtomInfo {
            name: "O".to_string(),
            elem_symbol: "O".to_string(),
            element: 8,
            resn: "HOH".to_string(),
            is_hetatm: true,
            ..AtomInfo::default()
        });
        mol.atoms.push(AtomInfo {
            name: "CA".to_string(),
            elem_symbol: "C".to_string(),
            element: 6,
            resn: "ALA".to_string(),
            ..AtomInfo::default()
        });

        assert_eq!(
            evaluate(&Selector::Hydrogen, &mol),
            vec![true, false, false]
        );
        assert_eq!(evaluate(&Selector::Solvent, &mol), vec![false, true, false]);
    }

    #[test]
    fn evaluate_polymer_selectors() {
        let mut mol = Molecule::new("test".into());
        mol.atoms.push(AtomInfo {
            resn: "ALA".to_string(),
            ..AtomInfo::default()
        });
        mol.atoms.push(AtomInfo {
            resn: "DA".to_string(),
            ..AtomInfo::default()
        });
        mol.atoms.push(AtomInfo {
            resn: "ATP".to_string(),
            is_hetatm: true,
            ..AtomInfo::default()
        });

        assert_eq!(evaluate(&Selector::Polymer, &mol), vec![true, true, false]);
        assert_eq!(evaluate(&Selector::Protein, &mol), vec![true, false, false]);
        assert_eq!(evaluate(&Selector::Nucleic, &mol), vec![false, true, false]);
    }

    #[test]
    fn evaluate_guide_selector() {
        let mut mol = Molecule::new("test".into());
        mol.atoms.push(AtomInfo {
            name: "CA".to_string(),
            resn: "ALA".to_string(),
            ..AtomInfo::default()
        });
        mol.atoms.push(AtomInfo {
            name: "CB".to_string(),
            resn: "ALA".to_string(),
            ..AtomInfo::default()
        });
        mol.atoms.push(AtomInfo {
            name: "C4'".to_string(),
            resn: "DA".to_string(),
            ..AtomInfo::default()
        });
        mol.atoms.push(AtomInfo {
            name: "C3*".to_string(),
            resn: "DG".to_string(),
            ..AtomInfo::default()
        });
        mol.atoms.push(AtomInfo {
            name: "CA".to_string(),
            resn: "ATP".to_string(),
            ..AtomInfo::default()
        });

        assert_eq!(
            evaluate(&Selector::Guide, &mol),
            vec![true, false, true, true, false]
        );
    }

    #[test]
    fn evaluate_metals_selector() {
        let mut mol = Molecule::new("test".into());
        mol.atoms.push(AtomInfo {
            element: 20,
            elem_symbol: "Ca".to_string(),
            ..AtomInfo::default()
        });
        mol.atoms.push(AtomInfo {
            element: 6,
            elem_symbol: "C".to_string(),
            ..AtomInfo::default()
        });
        mol.atoms.push(AtomInfo {
            element: 30,
            elem_symbol: "Zn".to_string(),
            ..AtomInfo::default()
        });

        assert_eq!(evaluate(&Selector::Metals, &mol), vec![true, false, true]);
    }

    #[test]
    fn evaluate_within_state_coords() {
        let mol = test_molecule();
        let sel = Selector::Within(0.2, Box::new(Selector::Serial(1, 10)));
        let first_state = evaluate_with_coords(&sel, &mol, mol.coord_sets.first().unwrap());
        let second_state = evaluate_with_coords(&sel, &mol, &mol.coord_sets[1]);
        assert_eq!(first_state, vec![true, false, false]);
        assert_eq!(second_state, vec![true, true, false]);
    }

    #[test]
    fn evaluate_beyond_and_near_to_selectors() {
        let mol = test_molecule();

        let beyond = Selector::Beyond(1.0, Box::new(Selector::Serial(10, 10)));
        assert_eq!(
            evaluate_with_coords(&beyond, &mol, &mol.coord_sets[0]),
            vec![false, true, true]
        );

        let near_to = Selector::NearTo(5.0, Box::new(Selector::Serial(10, 10)));
        assert_eq!(
            evaluate_with_coords(&near_to, &mol, &mol.coord_sets[0]),
            vec![false, true, false]
        );

        let beyond_empty = Selector::Beyond(1.0, Box::new(Selector::None));
        assert_eq!(
            evaluate_with_coords(&beyond_empty, &mol, &mol.coord_sets[0][..2]),
            vec![true, true, false]
        );
    }

    #[test]
    fn evaluate_around_expand_and_extend_selectors() {
        let mut mol = test_molecule();
        mol.bonds.push(crate::core::bond::BondInfo {
            atom_a: 0,
            atom_b: 1,
            order: 1,
        });
        mol.bonds.push(crate::core::bond::BondInfo {
            atom_a: 1,
            atom_b: 2,
            order: 1,
        });

        let around = Selector::Around(5.0, Box::new(Selector::Serial(10, 10)));
        assert_eq!(
            evaluate_with_coords(&around, &mol, &mol.coord_sets[0]),
            vec![false, true, false]
        );

        let expand = Selector::Expand(5.0, Box::new(Selector::Serial(10, 10)));
        assert_eq!(
            evaluate_with_coords(&expand, &mol, &mol.coord_sets[0]),
            vec![true, true, false]
        );

        let extend = Selector::Extend(2, Box::new(Selector::Serial(10, 10)));
        assert_eq!(evaluate(&extend, &mol), vec![true, true, true]);

        let extend_zero = Selector::Extend(0, Box::new(Selector::Serial(10, 10)));
        assert_eq!(evaluate(&extend_zero, &mol), vec![true, false, false]);
    }

    #[test]
    fn evaluate_logical_ops_preserve_coords() {
        let mol = test_molecule();
        let sel = Selector::And(
            Box::new(Selector::Serial(10, 12)),
            Box::new(Selector::Within(0.2, Box::new(Selector::Serial(10, 10)))),
        );
        let second_state = evaluate_with_coords(&sel, &mol, &mol.coord_sets[1]);
        assert_eq!(second_state, vec![true, true, false]);
    }

    #[test]
    fn evaluate_alt_selector() {
        let mol = test_molecule();
        assert_eq!(
            evaluate(&Selector::Alt('A'), &mol),
            vec![false, true, false]
        );
        assert_eq!(evaluate(&Selector::Alt(' '), &mol), vec![true, false, true]);
        assert_eq!(
            evaluate(&Selector::AltPattern("A+B".to_string()), &mol),
            vec![false, true, false]
        );
    }

    #[test]
    fn evaluate_ss_selector() {
        let mol = test_molecule();
        assert_eq!(
            evaluate(&Selector::SS("H".to_string()), &mol),
            vec![false, true, false]
        );
        assert_eq!(
            evaluate(&Selector::SS("SHEET".to_string()), &mol),
            vec![false, false, true]
        );
        assert_eq!(
            evaluate(&Selector::SS("LOOP".to_string()), &mol),
            vec![true, false, false]
        );
        assert_eq!(
            evaluate(&Selector::SS("H+S".to_string()), &mol),
            vec![false, true, true]
        );
        assert_eq!(
            evaluate(&Selector::SS("CO*".to_string()), &mol),
            vec![true, false, false]
        );
    }

    #[test]
    fn evaluate_neighbor_selector() {
        let mut mol = test_molecule();
        mol.bonds.push(crate::core::bond::BondInfo {
            atom_a: 0,
            atom_b: 1,
            order: 1,
        });
        mol.bonds.push(crate::core::bond::BondInfo {
            atom_a: 1,
            atom_b: 2,
            order: 1,
        });

        let sel = Selector::Neighbor(Box::new(Selector::Serial(10, 10)));
        assert_eq!(evaluate(&sel, &mol), vec![false, true, false]);

        let sel2 = Selector::Neighbor(Box::new(Selector::Serial(11, 12)));
        assert_eq!(evaluate(&sel2, &mol), vec![true, false, false]);
    }

    #[test]
    fn evaluate_byres_selector() {
        let mut mol = test_molecule();
        mol.atoms.push(AtomInfo {
            chain: 'A',
            resi: 1,
            ins_code: '\0',
            resn: "ALA".to_string(),
            name: "O".to_string(),
            ..AtomInfo::default()
        });

        let sel = Selector::Byres(Box::new(Selector::Serial(10, 10)));
        assert_eq!(evaluate(&sel, &mol), vec![true, true, false, true]);
    }

    #[test]
    fn evaluate_bycalpha_shape() {
        let mut mol = Molecule::new("test".into());
        mol.atoms.push(AtomInfo {
            serial: 1,
            chain: 'A',
            resi: 1,
            name: "N".to_string(),
            resn: "ALA".to_string(),
            ..AtomInfo::default()
        });
        mol.atoms.push(AtomInfo {
            serial: 2,
            chain: 'A',
            resi: 1,
            name: "CA".to_string(),
            resn: "ALA".to_string(),
            ..AtomInfo::default()
        });
        mol.atoms.push(AtomInfo {
            serial: 3,
            chain: 'A',
            resi: 2,
            name: "CA".to_string(),
            resn: "GLY".to_string(),
            ..AtomInfo::default()
        });

        let byca = Selector::And(
            Box::new(Selector::Byres(Box::new(Selector::Serial(1, 1)))),
            Box::new(Selector::Name("CA".to_string())),
        );
        assert_eq!(evaluate(&byca, &mol), vec![false, true, false]);
    }

    #[test]
    fn evaluate_bychain_byobject_and_bymolecule_selectors() {
        let mut mol = test_molecule();
        mol.bonds.push(crate::core::bond::BondInfo {
            atom_a: 0,
            atom_b: 1,
            order: 1,
        });

        let bychain = Selector::Bychain(Box::new(Selector::Serial(10, 10)));
        assert_eq!(evaluate(&bychain, &mol), vec![true, true, false]);

        let byobject = Selector::Byobject(Box::new(Selector::Serial(12, 12)));
        assert_eq!(evaluate(&byobject, &mol), vec![true, true, true]);

        let byobject_empty = Selector::Byobject(Box::new(Selector::None));
        assert_eq!(evaluate(&byobject_empty, &mol), vec![false, false, false]);

        let bymolecule = Selector::Bymolecule(Box::new(Selector::Serial(10, 10)));
        assert_eq!(evaluate(&bymolecule, &mol), vec![true, true, false]);

        let bymolecule_isolated = Selector::Bymolecule(Box::new(Selector::Serial(12, 12)));
        assert_eq!(
            evaluate(&bymolecule_isolated, &mol),
            vec![false, false, true]
        );
    }

    #[test]
    fn evaluate_first_last_selectors() {
        let mol = test_molecule();

        let first = Selector::First(Box::new(Selector::Chain('A')));
        assert_eq!(evaluate(&first, &mol), vec![true, false, false]);

        let last = Selector::Last(Box::new(Selector::Chain('A')));
        assert_eq!(evaluate(&last, &mol), vec![false, true, false]);

        let first_empty = Selector::First(Box::new(Selector::None));
        assert_eq!(evaluate(&first_empty, &mol), vec![false, false, false]);
    }
}
