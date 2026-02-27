use super::secondary_structure::SSType;

/// A contiguous range of atoms belonging to one residue.
#[derive(Debug, Clone)]
pub struct ResidueRange {
    pub chain: char,
    pub resn: String,
    pub resi: i32,
    pub ins_code: char,
    pub ss_type: SSType,
    /// Index of the first atom in the molecule's atom array.
    pub atom_start: usize,
    /// One-past-the-end index.
    pub atom_end: usize,
    /// Index of the CA atom within the molecule's atom array, if present.
    pub ca_index: Option<usize>,
}
