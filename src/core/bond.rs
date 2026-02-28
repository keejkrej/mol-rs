#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct BondInfo {
    /// Index of the first atom.
    pub atom_a: usize,
    /// Index of the second atom.
    pub atom_b: usize,
    /// Bond order (1 = single, 2 = double, 3 = triple, 4 = aromatic).
    pub order: u8,
}
