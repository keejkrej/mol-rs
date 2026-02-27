use crate::core::element::element_by_number;
use crate::core::molecule::Molecule;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColorScheme {
    ByElement,
    ByChain,
}

/// A fixed palette for chain coloring.
const CHAIN_COLORS: &[[f32; 3]] = &[
    [0.2, 0.8, 0.2],  // green
    [0.3, 0.5, 1.0],  // blue
    [1.0, 0.6, 0.2],  // orange
    [0.9, 0.2, 0.6],  // magenta
    [0.2, 0.9, 0.9],  // cyan
    [0.8, 0.8, 0.2],  // yellow
    [0.6, 0.3, 0.8],  // purple
    [0.9, 0.4, 0.4],  // salmon
    [0.4, 0.9, 0.5],  // lime
    [0.5, 0.7, 0.9],  // light blue
];

/// Apply a color scheme to all atoms in a molecule.
pub fn apply_color_scheme(mol: &mut Molecule, scheme: ColorScheme) {
    match scheme {
        ColorScheme::ByElement => {
            for atom in &mut mol.atoms {
                let elem = element_by_number(atom.element);
                atom.color = elem.color;
            }
        }
        ColorScheme::ByChain => {
            // Collect unique chains in order
            let mut chains: Vec<char> = Vec::new();
            for atom in &mol.atoms {
                if !chains.contains(&atom.chain) {
                    chains.push(atom.chain);
                }
            }
            for atom in &mut mol.atoms {
                let idx = chains.iter().position(|&c| c == atom.chain).unwrap_or(0);
                atom.color = CHAIN_COLORS[idx % CHAIN_COLORS.len()];
            }
        }
    }
}
