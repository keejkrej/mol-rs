use super::secondary_structure::SSType;

/// Representation visibility bitmask flags.
pub const REP_LINES: u32 = 1 << 0;
pub const REP_STICKS: u32 = 1 << 1;
pub const REP_SPHERES: u32 = 1 << 2;
pub const REP_CARTOON: u32 = 1 << 3;

/// Per-atom information, modeled after PyMOL's AtomInfoType.
#[derive(Debug, Clone)]
pub struct AtomInfo {
    /// Atom name, e.g. "CA", "N", "O" (PDB columns 13-16, trimmed).
    pub name: String,
    /// Atomic number (1 = H, 6 = C, 7 = N, 8 = O, etc.).
    pub element: u8,
    /// Element symbol string, e.g. "C", "N", "O".
    pub elem_symbol: String,
    /// Residue name, e.g. "ALA", "GLY".
    pub resn: String,
    /// Residue sequence number.
    pub resi: i32,
    /// Insertion code (PDB column 27); '\0' if none.
    pub ins_code: char,
    /// Chain identifier.
    pub chain: char,
    /// Alternate location indicator.
    pub alt: char,
    /// Secondary structure type assignment.
    pub ss_type: SSType,
    /// Isotropic temperature factor.
    pub b_factor: f32,
    /// Occupancy.
    pub occupancy: f32,
    /// Van der Waals radius (from element table).
    pub vdw: f32,
    /// Display color [r, g, b] in 0.0..1.0.
    pub color: [f32; 3],
    /// Bitmask of currently visible representations.
    pub vis_rep: u32,
    /// True if this atom came from a HETATM record.
    pub is_hetatm: bool,
    /// PDB serial number.
    pub serial: u32,
}

impl Default for AtomInfo {
    fn default() -> Self {
        Self {
            name: String::new(),
            element: 0,
            elem_symbol: String::new(),
            resn: String::new(),
            resi: 0,
            ins_code: '\0',
            chain: ' ',
            alt: ' ',
            ss_type: SSType::Loop,
            b_factor: 0.0,
            occupancy: 1.0,
            vdw: 1.7,
            color: [0.2, 1.0, 0.2], // default carbon green
            vis_rep: REP_LINES,
            is_hetatm: false,
            serial: 0,
        }
    }
}
