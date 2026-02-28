use super::atom::{AtomInfo, REP_CARTOON, REP_LINES, REP_STICKS};
use super::bond::BondInfo;
use super::residue::{is_nucleic, is_protein, ResidueRange};
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
                if atom.resn == "HOH" || atom.resn == "WAT" || atom.resn == "DOD" {
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
