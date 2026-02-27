use super::atom::AtomInfo;
use super::bond::BondInfo;
use super::residue::ResidueRange;
use super::secondary_structure::SSType;

/// A loaded molecular object, analogous to PyMOL's ObjectMolecule.
#[derive(Debug, Clone)]
pub struct Molecule {
    pub name: String,
    pub atoms: Vec<AtomInfo>,
    /// Coordinates as flat [x, y, z] per atom; coords[i] is the position of atoms[i].
    pub coords: Vec<[f32; 3]>,
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
            coords: Vec::new(),
            bonds: Vec::new(),
            residues: Vec::new(),
            visible: true,
        }
    }

    /// Compute the centroid of all atom coordinates.
    pub fn centroid(&self) -> [f32; 3] {
        if self.coords.is_empty() {
            return [0.0; 3];
        }
        let mut c = [0.0f32; 3];
        for p in &self.coords {
            c[0] += p[0];
            c[1] += p[1];
            c[2] += p[2];
        }
        let n = self.coords.len() as f32;
        [c[0] / n, c[1] / n, c[2] / n]
    }

    /// Compute the maximum distance from the centroid to any atom.
    pub fn radius(&self) -> f32 {
        let c = self.centroid();
        self.coords
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
}
