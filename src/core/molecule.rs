use std::cmp::Ordering;

use super::atom::{AtomInfo, REP_CARTOON, REP_LINES, REP_STICKS};
use super::bond::BondInfo;
use super::residue::{is_nucleic, is_protein, is_solvent, ResidueRange};
use super::secondary_structure::SSType;

/// A loaded molecular object, analogous to PyMOL's ObjectMolecule.
#[derive(Debug, Clone)]
pub struct Molecule {
    pub name: String,
    pub atoms: Vec<AtomInfo>,
    /// State-specific coordinates. coord_sets[0] is state 1.
    /// Each state vector aligns with `atoms` by index.
    pub coord_sets: Vec<Vec<[f32; 3]>>,
    pub bonds: Vec<BondInfo>,
    /// Precomputed residue groupings (filled after loading).
    pub residues: Vec<ResidueRange>,
    /// Whether this object is visible.
    pub visible: bool,
}

impl Molecule {
    pub fn new(name: String) -> Self {
        Self {
            name,
            atoms: Vec::new(),
            coord_sets: vec![Vec::new()],
            bonds: Vec::new(),
            residues: Vec::new(),
            visible: true,
        }
    }

    pub fn state_count(&self) -> usize {
        self.coord_sets.len()
    }

    pub fn coords_for_state(&self, state_1_based: usize) -> &[[f32; 3]] {
        if self.coord_sets.is_empty() {
            return &[];
        }
        let idx = state_1_based
            .saturating_sub(1)
            .min(self.coord_sets.len() - 1);
        &self.coord_sets[idx]
    }

    /// Remove atoms selected by mask and keep coordinates, bonds, and residues aligned.
    pub fn remove_atoms(&mut self, mask: &[bool]) -> usize {
        let old_len = self.atoms.len();
        if old_len == 0 {
            return 0;
        }

        let mut old_to_new = vec![None; old_len];
        let mut kept_atoms = Vec::with_capacity(old_len);
        let mut removed = 0usize;

        for (idx, atom) in self.atoms.iter().cloned().enumerate() {
            if mask.get(idx).copied().unwrap_or(false) {
                removed += 1;
            } else {
                old_to_new[idx] = Some(kept_atoms.len());
                kept_atoms.push(atom);
            }
        }

        if removed == 0 {
            return 0;
        }

        self.atoms = kept_atoms;

        for coords in &mut self.coord_sets {
            let mut kept_coords = Vec::with_capacity(coords.len().min(self.atoms.len()));
            for (idx, coord) in coords.iter().copied().enumerate() {
                if idx < old_len && old_to_new[idx].is_some() {
                    kept_coords.push(coord);
                }
            }
            *coords = kept_coords;
        }

        self.bonds = self
            .bonds
            .iter()
            .filter_map(|bond| {
                let atom_a = old_to_new.get(bond.atom_a).copied().flatten()?;
                let atom_b = old_to_new.get(bond.atom_b).copied().flatten()?;
                Some(BondInfo {
                    atom_a,
                    atom_b,
                    order: bond.order,
                })
            })
            .collect();

        self.build_residues();
        removed
    }

    /// Sort atoms in a PyMOL-like identifier order and keep topology aligned.
    pub fn sort_atoms(&mut self) -> bool {
        let old_len = self.atoms.len();
        if old_len < 2 {
            return false;
        }

        let mut order: Vec<usize> = (0..old_len).collect();
        order.sort_by(|&a, &b| {
            atom_sort_cmp(&self.atoms[a], &self.atoms[b]).then_with(|| a.cmp(&b))
        });

        if order
            .iter()
            .enumerate()
            .all(|(new_idx, old_idx)| new_idx == *old_idx)
        {
            return false;
        }

        let mut old_to_new = vec![0usize; old_len];
        for (new_idx, old_idx) in order.iter().copied().enumerate() {
            old_to_new[old_idx] = new_idx;
        }

        self.atoms = order
            .iter()
            .map(|old_idx| self.atoms[*old_idx].clone())
            .collect();

        for coords in &mut self.coord_sets {
            let old_coords = coords.clone();
            *coords = order
                .iter()
                .filter_map(|old_idx| old_coords.get(*old_idx).copied())
                .collect();
        }

        for bond in &mut self.bonds {
            bond.atom_a = old_to_new[bond.atom_a];
            bond.atom_b = old_to_new[bond.atom_b];
        }
        self.bonds.sort_by_key(|bond| {
            (
                bond.atom_a.min(bond.atom_b),
                bond.atom_a.max(bond.atom_b),
                bond.order,
            )
        });

        self.build_residues();
        true
    }

    /// Compute the centroid of atom coordinates in a specific state.
    pub fn centroid_for_state(&self, state_1_based: usize) -> [f32; 3] {
        let coords = self.coords_for_state(state_1_based);
        if coords.is_empty() {
            return [0.0; 3];
        }
        let mut c = [0.0f32; 3];
        for p in coords {
            c[0] += p[0];
            c[1] += p[1];
            c[2] += p[2];
        }
        let n = coords.len() as f32;
        [c[0] / n, c[1] / n, c[2] / n]
    }

    /// Compute the maximum distance from centroid in a specific state.
    pub fn radius_for_state(&self, state_1_based: usize) -> f32 {
        let coords = self.coords_for_state(state_1_based);
        let c = self.centroid_for_state(state_1_based);
        coords
            .iter()
            .map(|p| {
                let dx = p[0] - c[0];
                let dy = p[1] - c[1];
                let dz = p[2] - c[2];
                dx * dx + dy * dy + dz * dz
            })
            .fold(0.0f32, f32::max)
            .sqrt()
    }

    /// Build residue ranges from the atom array.
    /// Assumes atoms are already sorted by chain/resi as they come from PDB.
    pub fn build_residues(&mut self) {
        self.residues.clear();
        if self.atoms.is_empty() {
            return;
        }

        let mut start = 0usize;
        let mut ca_idx: Option<usize> = None;

        for i in 1..=self.atoms.len() {
            let new_residue = if i == self.atoms.len() {
                true
            } else {
                let prev = &self.atoms[i - 1];
                let curr = &self.atoms[i];
                prev.chain != curr.chain || prev.resi != curr.resi || prev.ins_code != curr.ins_code
            };

            // Track CA
            if i > start && i <= self.atoms.len() {
                let a = &self.atoms[i - 1];
                if a.name.trim() == "CA" {
                    ca_idx = Some(i - 1);
                }
            }

            if new_residue {
                let first = &self.atoms[start];
                self.residues.push(ResidueRange {
                    chain: first.chain,
                    resn: first.resn.clone(),
                    resi: first.resi,
                    ins_code: first.ins_code,
                    ss_type: first.ss_type,
                    atom_start: start,
                    atom_end: i,
                    ca_index: ca_idx,
                });
                start = i;
                ca_idx = None;
            }
        }
    }

    /// Apply secondary structure assignments (from HELIX/SHEET records) to atoms.
    pub fn apply_ss(&mut self, assignments: &[(char, i32, i32, SSType)]) {
        for atom in &mut self.atoms {
            for &(chain, start, end, ss) in assignments {
                if atom.chain == chain && atom.resi >= start && atom.resi <= end {
                    atom.ss_type = ss;
                    break;
                }
            }
        }
    }

    /// Apply a smart default representation based on residue type.
    /// - Proteins: Cartoon
    /// - Nucleic acids: Sticks (until cartoon is supported)
    /// - Ligands/HETATM: Sticks
    /// - Water: Lines (or hidden?)
    /// - Others: Lines
    pub fn apply_default_representation(&mut self) {
        for atom in &mut self.atoms {
            atom.vis_rep = 0; // Clear default

            if atom.is_hetatm {
                if is_solvent(&atom.resn) {
                    // Water - keep it simple
                    atom.vis_rep = REP_LINES;
                } else {
                    // Ligands / ions
                    atom.vis_rep = REP_STICKS;
                }
            } else {
                let is_prot = is_protein(&atom.resn);
                let is_nuc = is_nucleic(&atom.resn);

                if is_prot {
                    atom.vis_rep = REP_CARTOON;
                } else if is_nuc {
                    atom.vis_rep = REP_STICKS;
                } else {
                    atom.vis_rep = REP_LINES;
                }
            }
        }
    }
}

fn atom_sort_cmp(a: &AtomInfo, b: &AtomInfo) -> Ordering {
    a.chain
        .cmp(&b.chain)
        .then_with(|| a.is_hetatm.cmp(&b.is_hetatm))
        .then_with(|| a.resi.cmp(&b.resi))
        .then_with(|| a.ins_code.cmp(&b.ins_code))
        .then_with(|| a.resn.cmp(&b.resn))
        .then_with(|| atom_name_sort_cmp(&a.name, &b.name))
        .then_with(|| a.alt.cmp(&b.alt))
}

fn atom_name_sort_cmp(a: &str, b: &str) -> Ordering {
    let a_trimmed = a.trim();
    let b_trimmed = b.trim();
    let a_without_digit = a_trimmed
        .strip_prefix(|ch: char| ch.is_ascii_digit())
        .unwrap_or(a_trimmed);
    let b_without_digit = b_trimmed
        .strip_prefix(|ch: char| ch.is_ascii_digit())
        .unwrap_or(b_trimmed);

    a_without_digit
        .cmp(b_without_digit)
        .then_with(|| a_trimmed.cmp(b_trimmed))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atom(chain: char, resi: i32, name: &str) -> AtomInfo {
        AtomInfo {
            chain,
            resi,
            name: name.to_string(),
            ..AtomInfo::default()
        }
    }

    #[test]
    fn remove_atoms_keeps_topology_aligned() {
        let mut mol = Molecule::new("remove".to_string());
        mol.atoms = vec![
            atom('A', 1, "N"),
            atom('A', 1, "CA"),
            atom('A', 2, "N"),
            atom('A', 2, "CA"),
        ];
        mol.coord_sets = vec![
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [3.0, 0.0, 0.0],
            ],
            vec![
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
                [2.0, 1.0, 0.0],
                [3.0, 1.0, 0.0],
            ],
        ];
        mol.bonds = vec![
            BondInfo {
                atom_a: 0,
                atom_b: 1,
                order: 1,
            },
            BondInfo {
                atom_a: 1,
                atom_b: 2,
                order: 1,
            },
            BondInfo {
                atom_a: 2,
                atom_b: 3,
                order: 1,
            },
        ];
        mol.build_residues();

        let removed = mol.remove_atoms(&[false, true, false, false]);

        assert_eq!(removed, 1);
        assert_eq!(
            mol.atoms
                .iter()
                .map(|a| a.name.as_str())
                .collect::<Vec<_>>(),
            vec!["N", "N", "CA"]
        );
        assert_eq!(
            mol.coord_sets,
            vec![
                vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
                vec![[0.0, 1.0, 0.0], [2.0, 1.0, 0.0], [3.0, 1.0, 0.0]],
            ]
        );
        assert_eq!(mol.bonds.len(), 1);
        assert_eq!(mol.bonds[0].atom_a, 1);
        assert_eq!(mol.bonds[0].atom_b, 2);
        assert_eq!(mol.residues.len(), 2);
        assert_eq!(mol.residues[1].atom_start, 1);
        assert_eq!(mol.residues[1].atom_end, 3);
    }

    #[test]
    fn remove_atoms_noop_for_empty_mask() {
        let mut mol = Molecule::new("remove".to_string());
        mol.atoms = vec![atom('A', 1, "N")];
        mol.coord_sets = vec![vec![[0.0, 0.0, 0.0]]];
        mol.build_residues();

        assert_eq!(mol.remove_atoms(&[]), 0);
        assert_eq!(mol.atoms.len(), 1);
        assert_eq!(mol.coord_sets[0].len(), 1);
        assert_eq!(mol.residues.len(), 1);
    }

    #[test]
    fn sort_atoms_keeps_coordinates_bonds_and_residues_aligned() {
        let mut mol = Molecule::new("sort".to_string());
        mol.atoms = vec![
            atom('A', 2, "O"),
            atom('A', 1, "CA"),
            atom('A', 1, "N"),
            atom('B', 1, "C"),
        ];
        mol.coord_sets = vec![
            vec![
                [20.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [11.0, 0.0, 0.0],
                [30.0, 0.0, 0.0],
            ],
            vec![
                [20.0, 1.0, 0.0],
                [10.0, 1.0, 0.0],
                [11.0, 1.0, 0.0],
                [30.0, 1.0, 0.0],
            ],
        ];
        mol.bonds = vec![
            BondInfo {
                atom_a: 0,
                atom_b: 1,
                order: 1,
            },
            BondInfo {
                atom_a: 2,
                atom_b: 3,
                order: 1,
            },
        ];
        mol.build_residues();

        assert!(mol.sort_atoms());
        assert_eq!(
            mol.atoms
                .iter()
                .map(|a| (a.chain, a.resi, a.name.as_str()))
                .collect::<Vec<_>>(),
            vec![('A', 1, "CA"), ('A', 1, "N"), ('A', 2, "O"), ('B', 1, "C")]
        );
        assert_eq!(
            mol.coord_sets,
            vec![
                vec![
                    [10.0, 0.0, 0.0],
                    [11.0, 0.0, 0.0],
                    [20.0, 0.0, 0.0],
                    [30.0, 0.0, 0.0],
                ],
                vec![
                    [10.0, 1.0, 0.0],
                    [11.0, 1.0, 0.0],
                    [20.0, 1.0, 0.0],
                    [30.0, 1.0, 0.0],
                ],
            ]
        );
        let bond_pairs = mol
            .bonds
            .iter()
            .map(|bond| (bond.atom_a.min(bond.atom_b), bond.atom_a.max(bond.atom_b)))
            .collect::<Vec<_>>();
        assert_eq!(bond_pairs, vec![(0, 2), (1, 3)]);
        assert_eq!(mol.residues.len(), 3);
        assert_eq!(mol.residues[0].atom_start, 0);
        assert_eq!(mol.residues[0].atom_end, 2);
        assert!(!mol.sort_atoms());
    }
}
