use crate::core::element::element_by_number;
use crate::core::molecule::Molecule;
use crate::selection::parser::Selector;

/// Evaluate a selection expression against a molecule.
/// Returns a Vec<bool> with one entry per atom — true means selected.
pub fn evaluate(sel: &Selector, mol: &Molecule) -> Vec<bool> {
    let n = mol.atoms.len();
    match sel {
        Selector::All => vec![true; n],
        Selector::None => vec![false; n],
        Selector::Chain(ch) => mol.atoms.iter().map(|a| a.chain == *ch).collect(),
        Selector::Resi(lo, hi) => mol.atoms.iter().map(|a| a.resi >= *lo && a.resi <= *hi).collect(),
        Selector::Name(name) => {
            mol.atoms.iter().map(|a| a.name.trim().eq_ignore_ascii_case(name)).collect()
        }
        Selector::Resn(resn) => {
            mol.atoms.iter().map(|a| a.resn.trim().eq_ignore_ascii_case(resn)).collect()
        }
        Selector::Elem(sym) => {
            mol.atoms.iter().map(|a| {
                let ed = element_by_number(a.element);
                ed.symbol.eq_ignore_ascii_case(sym)
            }).collect()
        }
        Selector::Hetatm => mol.atoms.iter().map(|a| a.is_hetatm).collect(),
        Selector::And(left, right) => {
            let l = evaluate(left, mol);
            let r = evaluate(right, mol);
            l.iter().zip(r.iter()).map(|(a, b)| *a && *b).collect()
        }
        Selector::Or(left, right) => {
            let l = evaluate(left, mol);
            let r = evaluate(right, mol);
            l.iter().zip(r.iter()).map(|(a, b)| *a || *b).collect()
        }
        Selector::Not(inner) => {
            let v = evaluate(inner, mol);
            v.iter().map(|x| !x).collect()
        }
    }
}

/// Count how many atoms are selected.
pub fn count_selected(mask: &[bool]) -> usize {
    mask.iter().filter(|&&b| b).count()
}
