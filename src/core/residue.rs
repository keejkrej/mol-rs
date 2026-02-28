use super::secondary_structure::SSType;

/// A contiguous range of atoms belonging to one residue.
#[derive(Debug, Clone)]
#[allow(dead_code)]
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

pub fn is_protein(resn: &str) -> bool {
    matches!(
        resn.trim().to_uppercase().as_str(),
        "ALA" | "ARG" | "ASN" | "ASP" | "CYS" | "GLN" | "GLU" | "GLY" | "HIS" | "ILE" |
        "LEU" | "LYS" | "MET" | "PHE" | "PRO" | "SER" | "THR" | "TRP" | "TYR" | "VAL" |
        "ASX" | "GLX" | "UNK" | "MSE"
    )
}

pub fn is_nucleic(resn: &str) -> bool {
    matches!(
        resn.trim().to_uppercase().as_str(),
        "A" | "C" | "G" | "T" | "U" | "DA" | "DC" | "DG" | "DT" | "DU"
    )
}
