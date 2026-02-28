pub mod pdb;
pub mod cif;

use crate::core::molecule::Molecule;

pub struct LoadResult {
    pub molecule: Molecule,
    pub warnings: Vec<String>,
    pub source_model_count: usize,
}
