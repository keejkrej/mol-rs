use crate::core::atom::{
    NAMED_SELECTION_PROPERTY_PREFIX, REP_ALL, REP_CARTOON, REP_LINES, REP_SPHERES, REP_STICKS,
};
use crate::core::element::{element_by_number, is_metal_atomic_number};
use crate::core::molecule::Molecule;
use crate::core::residue::{is_nucleic, is_protein, is_solvent, protein_one_letter};
use crate::core::secondary_structure::SSType;
use crate::selection::parser::{AtomProperty, CompareOp, CustomPropertyOp, Selector};

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
        Selector::Named(name) => {
            let key = named_selection_property_key(name);
            mol.atoms
                .iter()
                .map(|atom| atom.properties.contains_key(&key))
                .collect()
        }
        Selector::Identifier(name) => {
            let key = named_selection_property_key(name);
            let object_selected = object_name_matches(name, &mol.name);
            mol.atoms
                .iter()
                .map(|atom| object_selected || atom.properties.contains_key(&key))
                .collect()
        }
        Selector::Present => (0..n).map(|idx| idx < coords.len()).collect(),
        Selector::State(state) => state_mask(*state, n, mol, coords),
        Selector::Bonded => bonded_atom_mask(mol),
        Selector::Donors => (0..n).map(|idx| is_hbond_donor(idx, mol)).collect(),
        Selector::Acceptors => (0..n).map(|idx| is_hbond_acceptor(idx, mol)).collect(),
        Selector::Delocalized => (0..n).map(|idx| is_delocalized(idx, mol)).collect(),
        Selector::Flag(flag) => {
            let mask = 1u32 << *flag;
            mol.atoms
                .iter()
                .map(|atom| atom.flags & mask != 0)
                .collect()
        }
        Selector::Masked => mol.atoms.iter().map(|atom| atom.masked).collect(),
        Selector::Protected => mol.atoms.iter().map(|atom| atom.protected).collect(),
        Selector::Chain(ch) => mol.atoms.iter().map(|a| a.chain == *ch).collect(),
        Selector::ChainPattern(pattern) => mol
            .atoms
            .iter()
            .map(|a| matches_alpha_list(&a.chain.to_string(), pattern))
            .collect(),
        Selector::Segi(pattern) => mol
            .atoms
            .iter()
            .map(|a| matches_alpha_list(a.segi.trim(), pattern))
            .collect(),
        Selector::Resi(lo, hi, ins_lo, ins_hi) => mol
            .atoms
            .iter()
            .map(|a| matches_resi_range(a.resi, a.ins_code, *lo, *hi, *ins_lo, *ins_hi))
            .collect(),
        Selector::ResiList(ranges) => mol
            .atoms
            .iter()
            .map(|a| in_resi_ranges(a.resi, a.ins_code, ranges))
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
        Selector::Pepseq(pattern) => pepseq_mask(pattern, mol),
        Selector::TextType(text_type) => mol
            .atoms
            .iter()
            .map(|a| matches_alpha_list(a.text_type.trim(), text_type))
            .collect(),
        Selector::NumericType(ranges) => mol
            .atoms
            .iter()
            .map(|a| in_numeric_ranges(a.numeric_type, ranges))
            .collect(),
        Selector::Custom(custom) => mol
            .atoms
            .iter()
            .map(|a| matches_alpha_list(a.custom.trim(), custom))
            .collect(),
        Selector::Label(label) => mol
            .atoms
            .iter()
            .map(|a| matches_alpha_list(a.label.trim(), label))
            .collect(),
        Selector::Stereo(stereo) => mol
            .atoms
            .iter()
            .map(|a| matches_alpha_list(a.stereo.trim(), stereo))
            .collect(),
        Selector::Organic => mol.atoms.iter().map(|a| !a.is_hetatm).collect(),
        Selector::Inorganic => mol.atoms.iter().map(|a| a.is_hetatm).collect(),
        Selector::Serial(lo, hi) => mol
            .atoms
            .iter()
            .enumerate()
            .map(|(idx, a)| {
                let id = atom_id(a, idx);
                id >= *lo && id <= *hi
            })
            .collect(),
        Selector::SerialList(ranges) => mol
            .atoms
            .iter()
            .enumerate()
            .map(|(idx, a)| in_numeric_ranges(atom_id(a, idx), ranges))
            .collect(),
        Selector::Index(lo, hi) => mol
            .atoms
            .iter()
            .enumerate()
            .map(|(idx, _)| {
                let index = i32::try_from(idx + 1).unwrap_or(0);
                index >= *lo && index <= *hi
            })
            .collect(),
        Selector::IndexList(ranges) => mol
            .atoms
            .iter()
            .enumerate()
            .map(|(idx, _)| {
                let index = i32::try_from(idx + 1).unwrap_or(0);
                in_numeric_ranges(index, ranges)
            })
            .collect(),
        Selector::Rank(lo, hi) => mol
            .atoms
            .iter()
            .map(|a| {
                let rank = i32::try_from(a.rank).unwrap_or(0);
                rank >= *lo && rank <= *hi
            })
            .collect(),
        Selector::RankList(ranges) => mol
            .atoms
            .iter()
            .map(|a| {
                let rank = i32::try_from(a.rank).unwrap_or(0);
                in_numeric_ranges(rank, ranges)
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
        Selector::CartoonColor(color) => {
            let rgb = color_rgb(color);
            mol.atoms
                .iter()
                .map(|atom| {
                    rgb.is_some_and(|rgb| {
                        atom.cartoon_color
                            .is_some_and(|cartoon_color| colors_match(cartoon_color, rgb))
                    })
                })
                .collect()
        }
        Selector::RibbonColor(color) => {
            let rgb = color_rgb(color);
            mol.atoms
                .iter()
                .map(|atom| {
                    rgb.is_some_and(|rgb| {
                        atom.ribbon_color
                            .is_some_and(|ribbon_color| colors_match(ribbon_color, rgb))
                    })
                })
                .collect()
        }
        Selector::Property(property, op, value) => mol
            .atoms
            .iter()
            .enumerate()
            .map(|(idx, atom)| {
                property_value(*property, atom, idx, mol, coords)
                    .is_some_and(|atom_value| compare_float(atom_value, *op, *value))
            })
            .collect(),
        Selector::CustomProperty(property, op, value) => mol
            .atoms
            .iter()
            .map(|atom| {
                atom.properties
                    .get(property)
                    .is_some_and(|atom_value| custom_property_matches(atom_value, *op, value))
            })
            .collect(),
        Selector::Byres(inner) => {
            let inner_mask = evaluate_with_coords(inner, mol, coords);
            let mut residues: Vec<(&str, char, i32, char)> = Vec::new();

            for (idx, selected) in inner_mask.iter().enumerate() {
                if *selected {
                    if let Some(atom) = mol.atoms.get(idx) {
                        let key = (atom.segi.as_str(), atom.chain, atom.resi, atom.ins_code);
                        if !residues.contains(&key) {
                            residues.push(key);
                        }
                    }
                }
            }

            mol.atoms
                .iter()
                .map(|atom| {
                    residues.contains(&(atom.segi.as_str(), atom.chain, atom.resi, atom.ins_code))
                })
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
        Selector::Bysegment(inner) => {
            let inner_mask = evaluate_with_coords(inner, mol, coords);
            let mut segments: Vec<&str> = Vec::new();

            for (idx, selected) in inner_mask.iter().enumerate() {
                if *selected {
                    if let Some(atom) = mol.atoms.get(idx) {
                        let segi = atom.segi.as_str();
                        if !segments.contains(&segi) {
                            segments.push(segi);
                        }
                    }
                }
            }

            mol.atoms
                .iter()
                .map(|atom| segments.contains(&atom.segi.as_str()))
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
        Selector::Byring(inner) => {
            let inner_mask = evaluate_with_coords(inner, mol, coords);
            ring_mask(&inner_mask, mol)
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
        Selector::BoundTo(inner) => {
            let inner_mask = evaluate_with_coords(inner, mol, coords);
            let mut out = vec![false; n];

            if inner_mask.iter().any(|&x| x) {
                for bond in &mol.bonds {
                    if bond.atom_a < n && bond.atom_b < n {
                        if inner_mask[bond.atom_a] {
                            out[bond.atom_b] = true;
                        }
                        if inner_mask[bond.atom_b] {
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
        Selector::Gap(distance, inner) => gap_mask(*distance, inner, mol, coords),
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
        Selector::Hydrogen => mol.atoms.iter().map(is_hydrogen_atom).collect(),
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
        Selector::In(left, right) => {
            let left_mask = evaluate_with_coords(left, mol, coords);
            let right_mask = evaluate_with_coords(right, mol, coords);
            atom_identity_mask(&left_mask, &right_mask, mol, atom_in_matches)
        }
        Selector::Like(left, right) => {
            let left_mask = evaluate_with_coords(left, mol, coords);
            let right_mask = evaluate_with_coords(right, mol, coords);
            atom_identity_mask(&left_mask, &right_mask, mol, atom_like_matches)
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
        !item.is_empty()
            && (if item.contains(':') {
                matches_alpha_range_ci(item, value)
            } else {
                wildcard_match_ci(item, value)
            })
    })
}

fn matches_alpha_range_ci(range: &str, value: &str) -> bool {
    let mut parts = range.split(':');
    let low = parts.next();
    let high = parts.next();
    if low.is_none() || high.is_none() || parts.next().is_some() {
        return false;
    }

    let low = low.unwrap_or("");
    let high = high.unwrap_or("");
    if low.is_empty()
        || high.is_empty()
        || low.contains('*')
        || high.contains('*')
        || low.contains('+')
        || high.contains('+')
    {
        return false;
    }

    let value = value.to_ascii_lowercase();
    let low = low.to_ascii_lowercase();
    let high = high.to_ascii_lowercase();
    low <= value && value <= high
}

pub fn named_selection_property_key(name: &str) -> String {
    format!("{}{}", NAMED_SELECTION_PROPERTY_PREFIX, name)
}

fn atom_identity_mask(
    left_mask: &[bool],
    right_mask: &[bool],
    mol: &Molecule,
    matches: fn(&crate::core::atom::AtomInfo, &crate::core::atom::AtomInfo) -> bool,
) -> Vec<bool> {
    let right_atoms: Vec<&crate::core::atom::AtomInfo> = mol
        .atoms
        .iter()
        .enumerate()
        .filter_map(|(idx, atom)| {
            right_mask
                .get(idx)
                .copied()
                .unwrap_or(false)
                .then_some(atom)
        })
        .collect();

    mol.atoms
        .iter()
        .enumerate()
        .map(|(idx, atom)| {
            left_mask.get(idx).copied().unwrap_or(false)
                && right_atoms
                    .iter()
                    .any(|right_atom| matches(atom, right_atom))
        })
        .collect()
}

fn atom_in_matches(
    left: &crate::core::atom::AtomInfo,
    right: &crate::core::atom::AtomInfo,
) -> bool {
    left.resi == right.resi
        && left.chain.eq_ignore_ascii_case(&right.chain)
        && left.name.eq_ignore_ascii_case(&right.name)
        && left.ins_code.eq_ignore_ascii_case(&right.ins_code)
        && left.resn.eq_ignore_ascii_case(&right.resn)
        && left.segi.eq_ignore_ascii_case(&right.segi)
}

fn atom_like_matches(
    left: &crate::core::atom::AtomInfo,
    right: &crate::core::atom::AtomInfo,
) -> bool {
    left.resi == right.resi
        && left.name.eq_ignore_ascii_case(&right.name)
        && left.ins_code.eq_ignore_ascii_case(&right.ins_code)
}

#[derive(Debug, Clone)]
struct ResidueSpan {
    resn: String,
    start: usize,
    end: usize,
}

fn pepseq_mask(pattern: &str, mol: &Molecule) -> Vec<bool> {
    let n = mol.atoms.len();
    let pattern: Vec<char> = pattern
        .trim()
        .chars()
        .filter(|ch| !ch.is_whitespace())
        .map(|ch| ch.to_ascii_uppercase())
        .collect();

    if pattern.is_empty() {
        return vec![false; n];
    }

    let residues = residue_spans(mol);
    let mut out = vec![false; n];
    if residues.len() < pattern.len() {
        return out;
    }

    for start in 0..=residues.len() - pattern.len() {
        if !pepseq_matches_at(&residues, &pattern, start) {
            continue;
        }

        for (offset, pattern_char) in pattern.iter().enumerate() {
            if *pattern_char == '-' {
                continue;
            }
            let residue = &residues[start + offset];
            for idx in residue.start..residue.end.min(n) {
                out[idx] = true;
            }
        }
    }

    out
}

fn pepseq_matches_at(residues: &[ResidueSpan], pattern: &[char], start: usize) -> bool {
    pattern.iter().enumerate().all(|(offset, pattern_char)| {
        if *pattern_char == '-' || *pattern_char == '+' {
            return true;
        }
        protein_one_letter(&residues[start + offset].resn)
            .is_some_and(|residue_char| residue_char == *pattern_char)
    })
}

fn residue_spans(mol: &Molecule) -> Vec<ResidueSpan> {
    let mut residues = Vec::new();
    if mol.atoms.is_empty() {
        return residues;
    }

    let mut start = 0usize;
    for idx in 1..=mol.atoms.len() {
        let new_residue = if idx == mol.atoms.len() {
            true
        } else {
            let previous = &mol.atoms[idx - 1];
            let current = &mol.atoms[idx];
            previous.segi != current.segi
                || previous.chain != current.chain
                || previous.resi != current.resi
                || previous.ins_code != current.ins_code
                || previous.resn != current.resn
        };

        if new_residue {
            residues.push(ResidueSpan {
                resn: mol.atoms[start].resn.clone(),
                start,
                end: idx,
            });
            start = idx;
        }
    }

    residues
}

fn in_numeric_ranges(value: i32, ranges: &[(i32, i32)]) -> bool {
    ranges.iter().any(|(lo, hi)| value >= *lo && value <= *hi)
}

fn state_mask(state: isize, atom_count: usize, mol: &Molecule, coords: &[[f32; 3]]) -> Vec<bool> {
    let state_coords = if state == -1 {
        Some(coords)
    } else if state >= 1 {
        mol.coord_sets
            .get((state - 1) as usize)
            .map(|state_coords| state_coords.as_slice())
    } else {
        None
    };

    (0..atom_count)
        .map(|idx| state_coords.is_some_and(|state_coords| idx < state_coords.len()))
        .collect()
}

fn atom_id(atom: &crate::core::atom::AtomInfo, idx: usize) -> i32 {
    let id = if atom.serial == 0 {
        u32::try_from(idx + 1).unwrap_or(0)
    } else {
        atom.serial
    };
    i32::try_from(id).unwrap_or(0)
}

fn in_resi_ranges(
    value: i32,
    ins_code: char,
    ranges: &[(i32, i32, Option<char>, Option<char>)],
) -> bool {
    ranges.iter().any(|(lo, hi, ins_lo, ins_hi)| {
        matches_resi_range(value, ins_code, *lo, *hi, *ins_lo, *ins_hi)
    })
}

fn matches_resi_range(
    value: i32,
    ins_code: char,
    lo: i32,
    hi: i32,
    ins_lo: Option<char>,
    ins_hi: Option<char>,
) -> bool {
    if value < lo || value > hi {
        return false;
    }

    if lo == hi {
        return match (ins_lo, ins_hi) {
            (None, None) => true,
            (Some(start), Some(end)) => ins_code >= start && ins_code <= end,
            (Some(start), None) => ins_code >= start,
            (None, Some(end)) => ins_code <= end,
        };
    }

    if value == lo {
        if let Some(start) = ins_lo {
            if ins_code < start {
                return false;
            }
        }
    }

    if value == hi {
        if let Some(end) = ins_hi {
            if ins_code > end {
                return false;
            }
        }
    }

    true
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
    mol: &Molecule,
    coords: &[[f32; 3]],
) -> Option<f32> {
    match property {
        AtomProperty::BFactor => Some(atom.b_factor),
        AtomProperty::Occupancy => Some(atom.occupancy),
        AtomProperty::FormalCharge => Some(atom.formal_charge as f32),
        AtomProperty::PartialCharge => Some(atom.partial_charge),
        AtomProperty::Vdw => Some(atom.vdw),
        AtomProperty::ElecRadius => Some(atom.elec_radius),
        AtomProperty::Cartoon => Some(atom.cartoon as f32),
        AtomProperty::Geom => Some(atom.geom as f32),
        AtomProperty::Valence => Some(atom.valence as f32),
        AtomProperty::Reps => Some(atom.vis_rep as f32),
        AtomProperty::Protons => Some(atom.element as f32),
        AtomProperty::Flags => Some(atom.flags as f32),
        AtomProperty::ExplicitDegree => Some(explicit_degree(idx, mol)),
        AtomProperty::ExplicitValence => Some(explicit_valence(idx, mol)),
        AtomProperty::X => coords.get(idx).map(|coord| coord[0]),
        AtomProperty::Y => coords.get(idx).map(|coord| coord[1]),
        AtomProperty::Z => coords.get(idx).map(|coord| coord[2]),
    }
}

fn custom_property_matches(atom_value: &str, op: CustomPropertyOp, selector_value: &str) -> bool {
    match op {
        CustomPropertyOp::In => matches_alpha_list(atom_value.trim(), selector_value),
        CustomPropertyOp::Greater
        | CustomPropertyOp::Less
        | CustomPropertyOp::Equal
        | CustomPropertyOp::GreaterEqual
        | CustomPropertyOp::LessEqual => {
            let Some(left) = atom_value.trim().parse::<f32>().ok() else {
                return false;
            };
            let Some(right) = selector_value.parse::<f32>().ok() else {
                return false;
            };
            match op {
                CustomPropertyOp::Greater => left > right,
                CustomPropertyOp::Less => left < right,
                CustomPropertyOp::Equal => (left - right).abs() < 0.0001,
                CustomPropertyOp::GreaterEqual => left >= right,
                CustomPropertyOp::LessEqual => left <= right,
                CustomPropertyOp::In => false,
            }
        }
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

fn gap_mask(distance: f32, inner: &Selector, mol: &Molecule, coords: &[[f32; 3]]) -> Vec<bool> {
    let n = mol.atoms.len();
    let inner_mask = evaluate_with_coords(inner, mol, coords);
    let clearance = distance.max(0.0);
    let mut centers: Vec<(usize, [f32; 3], f32)> = Vec::new();

    for (idx, selected) in inner_mask.iter().enumerate() {
        if *selected && idx < coords.len() {
            let vdw = mol.atoms.get(idx).map_or(0.0, |atom| atom.vdw);
            centers.push((idx, coords[idx], vdw));
        }
    }

    (0..n)
        .map(|idx| {
            if idx >= coords.len() || inner_mask.get(idx).copied().unwrap_or(false) {
                return false;
            }

            let atom_vdw = mol.atoms.get(idx).map_or(0.0, |atom| atom.vdw);
            let p = coords[idx];
            !centers.iter().any(|(_, center, center_vdw)| {
                let threshold = clearance + atom_vdw + center_vdw;
                let dx = p[0] - center[0];
                let dy = p[1] - center[1];
                let dz = p[2] - center[2];
                dx * dx + dy * dy + dz * dz <= threshold * threshold
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

fn is_hydrogen_atom(atom: &crate::core::atom::AtomInfo) -> bool {
    atom.element == 1
        || atom.elem_symbol.eq_ignore_ascii_case("H")
        || atom
            .name
            .trim_start_matches(char::is_numeric)
            .starts_with('H')
}

fn is_hbond_donor(idx: usize, mol: &Molecule) -> bool {
    let Some(atom) = mol.atoms.get(idx) else {
        return false;
    };
    matches!(atom.element, 7 | 8 | 16) && has_bonded_hydrogen(idx, mol)
}

fn is_hbond_acceptor(idx: usize, mol: &Molecule) -> bool {
    let Some(atom) = mol.atoms.get(idx) else {
        return false;
    };

    if atom.formal_charge > 0 || is_hydrogen_atom(atom) {
        return false;
    }

    match atom.element {
        8 | 16 => true,
        7 => !has_bonded_hydrogen(idx, mol),
        _ => false,
    }
}

fn has_bonded_hydrogen(idx: usize, mol: &Molecule) -> bool {
    mol.bonds.iter().filter(|bond| bond.order > 0).any(|bond| {
        let other = if bond.atom_a == idx {
            Some(bond.atom_b)
        } else if bond.atom_b == idx {
            Some(bond.atom_a)
        } else {
            None
        };

        other
            .and_then(|other| mol.atoms.get(other))
            .is_some_and(is_hydrogen_atom)
    })
}

fn explicit_degree(idx: usize, mol: &Molecule) -> f32 {
    mol.bonds
        .iter()
        .filter(|bond| bond.order > 0 && (bond.atom_a == idx || bond.atom_b == idx))
        .count() as f32
}

fn explicit_valence(idx: usize, mol: &Molecule) -> f32 {
    mol.bonds
        .iter()
        .filter(|bond| bond.order > 0 && (bond.atom_a == idx || bond.atom_b == idx))
        .map(|bond| match bond.order {
            4 => 1.5,
            order => order as f32,
        })
        .sum()
}

fn is_delocalized(idx: usize, mol: &Molecule) -> bool {
    mol.bonds
        .iter()
        .any(|bond| bond.order == 4 && (bond.atom_a == idx || bond.atom_b == idx))
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

fn ring_mask(inner_mask: &[bool], mol: &Molecule) -> Vec<bool> {
    const MAX_RING_SIZE: usize = 7;

    let n = mol.atoms.len();
    let adjacency = bond_adjacency(mol);
    let mut out = vec![false; n];

    for start in inner_mask
        .iter()
        .enumerate()
        .filter_map(|(idx, selected)| (*selected && idx < n).then_some(idx))
    {
        let mut visited = vec![false; n];
        let mut path = vec![start];
        visited[start] = true;
        collect_rings_from(
            start,
            start,
            &adjacency,
            &mut visited,
            &mut path,
            &mut out,
            MAX_RING_SIZE,
        );
    }

    out
}

fn collect_rings_from(
    start: usize,
    current: usize,
    adjacency: &[Vec<usize>],
    visited: &mut [bool],
    path: &mut Vec<usize>,
    out: &mut [bool],
    max_ring_size: usize,
) {
    for &next in &adjacency[current] {
        if next == start {
            if path.len() >= 3 {
                for &idx in path.iter() {
                    out[idx] = true;
                }
            }
        } else if path.len() < max_ring_size && !visited[next] {
            visited[next] = true;
            path.push(next);
            collect_rings_from(start, next, adjacency, visited, path, out, max_ring_size);
            path.pop();
            visited[next] = false;
        }
    }
}

fn bond_adjacency(mol: &Molecule) -> Vec<Vec<usize>> {
    let n = mol.atoms.len();
    let mut adjacency = vec![Vec::new(); n];

    for bond in mol.bonds.iter().filter(|bond| bond.order > 0) {
        if bond.atom_a < n && bond.atom_b < n {
            adjacency[bond.atom_a].push(bond.atom_b);
            adjacency[bond.atom_b].push(bond.atom_a);
        }
    }

    adjacency
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
    use crate::core::atom::{AtomInfo, REP_LINES, REP_STICKS};
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
    fn evaluate_flags_masked_and_protected() {
        let mut mol = test_molecule();
        mol.atoms[0].flags = 1 << 25;
        mol.atoms[1].flags = (1 << 31) | (1 << 3);
        mol.atoms[2].flags = 1 << 2;
        mol.atoms[1].masked = true;
        mol.atoms[2].protected = true;

        assert_eq!(
            evaluate(&Selector::Flag(25), &mol),
            vec![true, false, false]
        );
        assert_eq!(
            evaluate(&Selector::Flag(31), &mol),
            vec![false, true, false]
        );
        assert_eq!(evaluate(&Selector::Flag(3), &mol), vec![false, true, false]);
        assert_eq!(evaluate(&Selector::Flag(2), &mol), vec![false, false, true]);
        assert_eq!(evaluate(&Selector::Masked, &mol), vec![false, true, false]);
        assert_eq!(
            evaluate(&Selector::Protected, &mol),
            vec![false, false, true]
        );
    }

    #[test]
    fn evaluate_serial_range() {
        let mut mol = test_molecule();
        mol.atoms[0].rank = 3;
        mol.atoms[1].rank = 1;
        mol.atoms[2].rank = 2;
        let serial = evaluate(&Selector::Serial(11, 12), &mol);
        assert_eq!(serial, vec![false, true, true]);
        let index = evaluate(&Selector::Index(1, 2), &mol);
        assert_eq!(index, vec![true, true, false]);
        let rank = evaluate(&Selector::Rank(1, 2), &mol);
        assert_eq!(rank, vec![false, true, true]);

        mol.atoms[1].serial = 0;
        let id_fallback = evaluate(&Selector::Serial(2, 2), &mol);
        assert_eq!(id_fallback, vec![false, true, false]);
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
        assert_eq!(
            evaluate(&Selector::ChainPattern("A:C".to_string()), &mol),
            vec![true, true, true]
        );
        assert_eq!(
            evaluate(&Selector::ChainPattern("A:C+Z*".to_string()), &mol),
            vec![true, true, true]
        );
    }

    #[test]
    fn evaluate_segi_selector() {
        let mut mol = test_molecule();
        mol.atoms[0].segi = "PROA".to_string();
        mol.atoms[1].segi = "PROA".to_string();
        mol.atoms[2].segi = "LIG".to_string();

        assert_eq!(
            evaluate(&Selector::Segi("PRO*".to_string()), &mol),
            vec![true, true, false]
        );
        assert_eq!(
            evaluate(&Selector::Segi("PROA+LIG".to_string()), &mol),
            vec![true, true, true]
        );
    }

    #[test]
    fn evaluate_plus_separated_numeric_lists() {
        let mut mol = test_molecule();
        mol.atoms[0].numeric_type = 10;
        mol.atoms[1].numeric_type = 20;
        mol.atoms[2].numeric_type = 25;

        assert_eq!(
            evaluate(
                &Selector::ResiList(vec![(2, 2, None, None), (4, 4, None, None)]),
                &mol
            ),
            vec![false, false, true]
        );
        assert_eq!(
            evaluate(&Selector::SerialList(vec![(10, 10), (12, 12)]), &mol),
            vec![true, false, true]
        );
        assert_eq!(
            evaluate(&Selector::IndexList(vec![(1, 1), (3, 3)]), &mol),
            vec![true, false, true]
        );
        assert_eq!(
            evaluate(&Selector::NumericType(vec![(10, 10), (20, 30)]), &mol),
            vec![true, true, true]
        );
        assert_eq!(
            evaluate(&Selector::NumericType(vec![(-9999, -9999)]), &mol),
            vec![false, false, false]
        );
    }

    #[test]
    fn evaluate_resi_with_insertion_code() {
        let mut mol = Molecule::new("test".into());
        mol.atoms.push(AtomInfo {
            resi: 9,
            ins_code: 'A',
            ..AtomInfo::default()
        });
        mol.atoms.push(AtomInfo {
            resi: 9,
            ins_code: 'B',
            ..AtomInfo::default()
        });
        mol.atoms.push(AtomInfo {
            resi: 9,
            ..AtomInfo::default()
        });
        mol.atoms.push(AtomInfo {
            resi: 10,
            ins_code: 'A',
            ..AtomInfo::default()
        });

        assert_eq!(
            evaluate(&Selector::Resi(9, 9, None, None), &mol),
            vec![true, true, true, false]
        );
        assert_eq!(
            evaluate(&Selector::Resi(9, 9, Some('A'), Some('A')), &mol),
            vec![true, false, false, false]
        );
        assert_eq!(
            evaluate(&Selector::Resi(9, 9, Some('B'), Some('B')), &mol),
            vec![false, true, false, false]
        );
        assert_eq!(
            evaluate(
                &Selector::ResiList(vec![(9, 9, Some('A'), Some('A')), (10, 10, None, None)]),
                &mol
            ),
            vec![true, false, false, true]
        );
        assert_eq!(
            evaluate(&Selector::Resi(9, 10, Some('A'), Some('A')), &mol),
            vec![true, true, false, true]
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
        mol.atoms[0].cartoon_color = Some([0.2, 0.2, 1.0]);
        mol.atoms[1].ribbon_color = Some([1.0, 0.2, 0.2]);
        mol.atoms[2].cartoon_color = Some([0.5, 0.5, 0.5]);

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
        assert_eq!(
            evaluate(&Selector::CartoonColor("blue".to_string()), &mol),
            vec![true, false, false]
        );
        assert_eq!(
            evaluate(&Selector::CartoonColor("grey".to_string()), &mol),
            vec![false, false, true]
        );
        assert_eq!(
            evaluate(&Selector::RibbonColor("red".to_string()), &mol),
            vec![false, true, false]
        );
        assert_eq!(
            evaluate(&Selector::RibbonColor("blue".to_string()), &mol),
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
        mol.atoms[0].formal_charge = -1;
        mol.atoms[1].formal_charge = 0;
        mol.atoms[2].formal_charge = 1;
        mol.atoms[0].partial_charge = -0.3;
        mol.atoms[1].partial_charge = 0.0;
        mol.atoms[2].partial_charge = 0.2;
        mol.atoms[0].vdw = 1.5;
        mol.atoms[1].vdw = 1.7;
        mol.atoms[2].vdw = 2.0;
        mol.atoms[0].elec_radius = 1.1;
        mol.atoms[1].elec_radius = 1.3;
        mol.atoms[2].elec_radius = 1.5;
        mol.atoms[0].cartoon = 0;
        mol.atoms[1].cartoon = 2;
        mol.atoms[2].cartoon = 3;
        mol.atoms[0].geom = 1;
        mol.atoms[1].geom = 3;
        mol.atoms[2].geom = 4;
        mol.atoms[0].valence = 1;
        mol.atoms[1].valence = 3;
        mol.atoms[2].valence = 4;
        mol.atoms[0].vis_rep = REP_LINES;
        mol.atoms[1].vis_rep = REP_STICKS;
        mol.atoms[2].vis_rep = REP_LINES | REP_STICKS;
        mol.atoms[0].flags = 0;
        mol.atoms[1].flags = 1 << 2;
        mol.atoms[2].flags = 1 << 3;
        mol.bonds.push(crate::core::bond::BondInfo {
            atom_a: 0,
            atom_b: 1,
            order: 1,
        });
        mol.bonds.push(crate::core::bond::BondInfo {
            atom_a: 1,
            atom_b: 2,
            order: 2,
        });

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
            evaluate(
                &Selector::Property(AtomProperty::FormalCharge, CompareOp::Equal, -1.0),
                &mol,
            ),
            vec![true, false, false]
        );
        assert_eq!(
            evaluate(
                &Selector::Property(AtomProperty::FormalCharge, CompareOp::GreaterEqual, 0.0),
                &mol,
            ),
            vec![false, true, true]
        );
        assert_eq!(
            evaluate(
                &Selector::Property(AtomProperty::PartialCharge, CompareOp::Less, -0.1),
                &mol,
            ),
            vec![true, false, false]
        );
        assert_eq!(
            evaluate(
                &Selector::Property(AtomProperty::PartialCharge, CompareOp::Greater, 0.1),
                &mol,
            ),
            vec![false, false, true]
        );
        assert_eq!(
            evaluate(
                &Selector::Property(AtomProperty::Vdw, CompareOp::LessEqual, 1.7),
                &mol,
            ),
            vec![true, true, false]
        );
        assert_eq!(
            evaluate(
                &Selector::Property(AtomProperty::ElecRadius, CompareOp::Greater, 1.2),
                &mol,
            ),
            vec![false, true, true]
        );
        assert_eq!(
            evaluate(
                &Selector::Property(AtomProperty::Cartoon, CompareOp::GreaterEqual, 2.0),
                &mol,
            ),
            vec![false, true, true]
        );
        assert_eq!(
            evaluate(
                &Selector::Property(AtomProperty::Geom, CompareOp::Equal, 3.0),
                &mol,
            ),
            vec![false, true, false]
        );
        assert_eq!(
            evaluate(
                &Selector::Property(AtomProperty::Valence, CompareOp::Less, 4.0),
                &mol,
            ),
            vec![true, true, false]
        );
        assert_eq!(
            evaluate(
                &Selector::Property(AtomProperty::Reps, CompareOp::Equal, REP_STICKS as f32),
                &mol,
            ),
            vec![false, true, false]
        );
        assert_eq!(
            evaluate(
                &Selector::Property(AtomProperty::Protons, CompareOp::Equal, 6.0),
                &mol,
            ),
            vec![false, true, true]
        );
        assert_eq!(
            evaluate(
                &Selector::Property(AtomProperty::Flags, CompareOp::Greater, 0.0),
                &mol,
            ),
            vec![false, true, true]
        );
        assert_eq!(
            evaluate(
                &Selector::Property(AtomProperty::ExplicitDegree, CompareOp::GreaterEqual, 2.0),
                &mol,
            ),
            vec![false, true, false]
        );
        assert_eq!(
            evaluate(
                &Selector::Property(AtomProperty::ExplicitValence, CompareOp::Equal, 3.0),
                &mol,
            ),
            vec![false, true, false]
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
    fn evaluate_hbond_and_delocalized_selectors() {
        let mut mol = Molecule::new("chem".to_string());
        mol.atoms = vec![
            AtomInfo {
                name: "N".to_string(),
                elem_symbol: "N".to_string(),
                element: 7,
                ..AtomInfo::default()
            },
            AtomInfo {
                name: "H".to_string(),
                elem_symbol: "H".to_string(),
                element: 1,
                ..AtomInfo::default()
            },
            AtomInfo {
                name: "O".to_string(),
                elem_symbol: "O".to_string(),
                element: 8,
                ..AtomInfo::default()
            },
            AtomInfo {
                name: "C1".to_string(),
                elem_symbol: "C".to_string(),
                element: 6,
                ..AtomInfo::default()
            },
            AtomInfo {
                name: "C2".to_string(),
                elem_symbol: "C".to_string(),
                element: 6,
                ..AtomInfo::default()
            },
        ];
        mol.coord_sets = vec![vec![[0.0; 3]; mol.atoms.len()]];
        mol.bonds = vec![
            crate::core::bond::BondInfo {
                atom_a: 0,
                atom_b: 1,
                order: 1,
            },
            crate::core::bond::BondInfo {
                atom_a: 3,
                atom_b: 4,
                order: 4,
            },
        ];

        assert_eq!(
            evaluate(&Selector::Donors, &mol),
            vec![true, false, false, false, false]
        );
        assert_eq!(
            evaluate(&Selector::Acceptors, &mol),
            vec![false, false, true, false, false]
        );
        assert_eq!(
            evaluate(&Selector::Delocalized, &mol),
            vec![false, false, false, true, true]
        );
        assert_eq!(
            evaluate(
                &Selector::Property(AtomProperty::ExplicitValence, CompareOp::Equal, 1.5),
                &mol,
            ),
            vec![false, false, false, true, true]
        );
    }

    #[test]
    fn evaluate_custom_property_selectors() {
        let mut mol = test_molecule();
        mol.atoms[0]
            .properties
            .insert("score".to_string(), "0.75".to_string());
        mol.atoms[1]
            .properties
            .insert("score".to_string(), "0.25".to_string());
        mol.atoms[2]
            .properties
            .insert("score".to_string(), "inactive".to_string());
        mol.atoms[0]
            .properties
            .insert("kind".to_string(), "ligand_core".to_string());
        mol.atoms[1]
            .properties
            .insert("kind".to_string(), "ligand_tail".to_string());
        mol.atoms[2]
            .properties
            .insert("kind".to_string(), "solvent".to_string());

        assert_eq!(
            evaluate(
                &Selector::CustomProperty(
                    "score".to_string(),
                    CustomPropertyOp::GreaterEqual,
                    "0.5".to_string(),
                ),
                &mol,
            ),
            vec![true, false, false]
        );
        assert_eq!(
            evaluate(
                &Selector::CustomProperty(
                    "score".to_string(),
                    CustomPropertyOp::Less,
                    "0.5".to_string(),
                ),
                &mol,
            ),
            vec![false, true, false]
        );
        assert_eq!(
            evaluate(
                &Selector::CustomProperty(
                    "kind".to_string(),
                    CustomPropertyOp::In,
                    "ligand_*".to_string(),
                ),
                &mol,
            ),
            vec![true, true, false]
        );
        assert_eq!(
            evaluate(
                &Selector::CustomProperty(
                    "missing".to_string(),
                    CustomPropertyOp::In,
                    "*".to_string(),
                ),
                &mol,
            ),
            vec![false, false, false]
        );
    }

    #[test]
    fn evaluate_plus_separated_alpha_lists() {
        let mut mol = test_molecule();
        mol.atoms[0].text_type = "CT".to_string();
        mol.atoms[1].text_type = "HC".to_string();
        mol.atoms[2].text_type = "OW".to_string();
        mol.atoms[0].custom = "ligand_core".to_string();
        mol.atoms[1].custom = "ligand_tail".to_string();
        mol.atoms[2].custom = "solvent".to_string();
        mol.atoms[0].label = "active_site".to_string();
        mol.atoms[1].label = "active_ligand".to_string();
        mol.atoms[2].label = "bulk".to_string();
        mol.atoms[0].stereo = "R".to_string();
        mol.atoms[1].stereo = "S".to_string();
        mol.atoms[2].stereo = "odd".to_string();

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
        assert_eq!(
            evaluate(&Selector::TextType("CT+HC".to_string()), &mol),
            vec![true, true, false]
        );
        assert_eq!(
            evaluate(&Selector::Custom("ligand_*".to_string()), &mol),
            vec![true, true, false]
        );
        assert_eq!(
            evaluate(&Selector::Label("active_*".to_string()), &mol),
            vec![true, true, false]
        );
        assert_eq!(
            evaluate(&Selector::Stereo("R+S".to_string()), &mol),
            vec![true, true, false]
        );
    }

    #[test]
    fn evaluate_pepseq_selector() {
        let mut mol = Molecule::new("pep".to_string());
        let residues = [("ALA", 1), ("GLY", 2), ("SER", 3), ("ALA", 4)];
        for (resn, resi) in residues {
            mol.atoms.push(AtomInfo {
                name: "N".to_string(),
                element: 7,
                resn: resn.to_string(),
                resi,
                chain: 'A',
                ..AtomInfo::default()
            });
            mol.atoms.push(AtomInfo {
                name: "CA".to_string(),
                element: 6,
                resn: resn.to_string(),
                resi,
                chain: 'A',
                ..AtomInfo::default()
            });
        }
        mol.coord_sets = vec![vec![[0.0; 3]; mol.atoms.len()]];

        assert_eq!(
            evaluate(&Selector::Pepseq("AG".to_string()), &mol),
            vec![true, true, true, true, false, false, false, false]
        );
        assert_eq!(
            evaluate(&Selector::Pepseq("A+S".to_string()), &mol),
            vec![true, true, true, true, true, true, false, false]
        );
        assert_eq!(
            evaluate(&Selector::Pepseq("A-S".to_string()), &mol),
            vec![true, true, false, false, true, true, false, false]
        );
        assert_eq!(
            evaluate(&Selector::Pepseq("WG".to_string()), &mol),
            vec![false; 8]
        );
    }

    #[test]
    fn evaluate_wildcard_alpha_lists() {
        let mut mol = test_molecule();
        mol.atoms[0].text_type = "CT".to_string();
        mol.atoms[1].text_type = "HC".to_string();
        mol.atoms[2].text_type = "OW".to_string();

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
        assert_eq!(
            evaluate(&Selector::TextType("?C".to_string()), &mol),
            vec![false, false, false]
        );
        assert_eq!(
            evaluate(&Selector::TextType("*C".to_string()), &mol),
            vec![false, true, false]
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
    fn evaluate_named_selection_selector() {
        let mut mol = test_molecule();
        mol.name = "obj".to_string();
        let key = named_selection_property_key("stored");
        mol.atoms[0].properties.insert(key.clone(), "1".to_string());
        mol.atoms[2].properties.insert(key, "1".to_string());

        assert_eq!(
            evaluate(&Selector::Named("stored".to_string()), &mol),
            vec![true, false, true]
        );
        assert_eq!(
            evaluate(&Selector::Named("missing".to_string()), &mol),
            vec![false, false, false]
        );
        assert_eq!(
            evaluate(&Selector::Identifier("stored".to_string()), &mol),
            vec![true, false, true]
        );
        assert_eq!(
            evaluate(&Selector::Identifier("obj".to_string()), &mol),
            vec![true, true, true]
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
    fn evaluate_state_selector() {
        let mol = test_molecule();

        assert_eq!(evaluate(&Selector::State(1), &mol), vec![true, true, true]);
        assert_eq!(evaluate(&Selector::State(2), &mol), vec![true, true, true]);
        assert_eq!(
            evaluate_with_coords(&Selector::State(-1), &mol, &mol.coord_sets[0][..2]),
            vec![true, true, false]
        );
        assert_eq!(
            evaluate(&Selector::State(3), &mol),
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
    fn evaluate_gap_selector_uses_vdw_clearance() {
        let mut mol = test_molecule();
        mol.atoms[0].vdw = 1.0;
        mol.atoms[1].vdw = 1.0;
        mol.atoms[2].vdw = 1.0;
        mol.coord_sets = vec![vec![[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [4.5, 0.0, 0.0]]];

        let gap = Selector::Gap(1.0, Box::new(Selector::Serial(10, 10)));
        assert_eq!(evaluate(&gap, &mol), vec![false, false, true]);

        let no_extra_gap = Selector::Gap(0.0, Box::new(Selector::Serial(10, 10)));
        assert_eq!(evaluate(&no_extra_gap, &mol), vec![false, true, true]);
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
    fn evaluate_in_like_and_subtract_operators() {
        let mut mol = Molecule::new("identity".to_string());
        mol.atoms = vec![
            AtomInfo {
                name: "CA".to_string(),
                resn: "ALA".to_string(),
                resi: 1,
                chain: 'A',
                segi: "S1".to_string(),
                serial: 10,
                ..AtomInfo::default()
            },
            AtomInfo {
                name: "CA".to_string(),
                resn: "ALA".to_string(),
                resi: 2,
                chain: 'A',
                segi: "S1".to_string(),
                serial: 11,
                ..AtomInfo::default()
            },
            AtomInfo {
                name: "CA".to_string(),
                resn: "ALA".to_string(),
                resi: 1,
                chain: 'B',
                segi: "S2".to_string(),
                serial: 12,
                ..AtomInfo::default()
            },
            AtomInfo {
                name: "CA".to_string(),
                resn: "ALA".to_string(),
                resi: 1,
                chain: 'A',
                segi: "S1".to_string(),
                serial: 13,
                ..AtomInfo::default()
            },
        ];
        mol.coord_sets = vec![vec![[0.0; 3]; mol.atoms.len()]];

        let in_sel = Selector::In(
            Box::new(Selector::SerialList(vec![(10, 12)])),
            Box::new(Selector::Serial(13, 13)),
        );
        assert_eq!(evaluate(&in_sel, &mol), vec![true, false, false, false]);

        let like_sel = Selector::Like(
            Box::new(Selector::SerialList(vec![(10, 12)])),
            Box::new(Selector::Serial(13, 13)),
        );
        assert_eq!(evaluate(&like_sel, &mol), vec![true, false, true, false]);

        let subtract_sel = Selector::And(
            Box::new(Selector::All),
            Box::new(Selector::Not(Box::new(Selector::Chain('A')))),
        );
        assert_eq!(
            evaluate(&subtract_sel, &mol),
            vec![false, false, true, false]
        );
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

        let bound_to = Selector::BoundTo(Box::new(Selector::Serial(11, 12)));
        assert_eq!(evaluate(&bound_to, &mol), vec![true, true, true]);
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
    fn evaluate_byres_keeps_segments_separate() {
        let mut mol = Molecule::new("segments".to_string());
        mol.atoms.push(AtomInfo {
            serial: 1,
            chain: 'A',
            segi: "A".to_string(),
            resi: 1,
            ..AtomInfo::default()
        });
        mol.atoms.push(AtomInfo {
            serial: 2,
            chain: 'A',
            segi: "B".to_string(),
            resi: 1,
            ..AtomInfo::default()
        });

        assert_eq!(
            evaluate(&Selector::Byres(Box::new(Selector::Serial(1, 1))), &mol),
            vec![true, false]
        );
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

        mol.atoms[0].segi = "PROA".to_string();
        mol.atoms[1].segi = "PROA".to_string();
        mol.atoms[2].segi = "LIG".to_string();
        let bysegment = Selector::Bysegment(Box::new(Selector::Serial(10, 10)));
        assert_eq!(evaluate(&bysegment, &mol), vec![true, true, false]);

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
    fn evaluate_byring_selector() {
        let mut mol = Molecule::new("rings".into());
        for serial in 1..=10 {
            mol.atoms.push(AtomInfo {
                serial,
                ..AtomInfo::default()
            });
        }

        let bonds = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (0, 6),
            (6, 7),
            (7, 8),
            (8, 6),
            (8, 9),
        ];
        mol.bonds = bonds
            .iter()
            .map(|&(atom_a, atom_b)| crate::core::bond::BondInfo {
                atom_a,
                atom_b,
                order: 1,
            })
            .collect();

        let benzene_like = Selector::Byring(Box::new(Selector::Serial(1, 1)));
        assert_eq!(
            evaluate(&benzene_like, &mol),
            vec![true, true, true, true, true, true, false, false, false, false]
        );

        let three_membered = Selector::Byring(Box::new(Selector::Serial(8, 8)));
        assert_eq!(
            evaluate(&three_membered, &mol),
            vec![false, false, false, false, false, false, true, true, true, false]
        );

        let acyclic = Selector::Byring(Box::new(Selector::Serial(10, 10)));
        assert_eq!(evaluate(&acyclic, &mol), vec![false; 10]);
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
