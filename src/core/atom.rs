use std::collections::HashMap;

use super::secondary_structure::SSType;

/// Representation visibility bitmask flags.
pub const REP_LINES: u32 = 1 << 0;
pub const REP_STICKS: u32 = 1 << 1;
pub const REP_SPHERES: u32 = 1 << 2;
pub const REP_CARTOON: u32 = 1 << 3;
pub const REP_ALL: u32 = REP_LINES | REP_STICKS | REP_SPHERES | REP_CARTOON;
pub const NAMED_SELECTION_PROPERTY_PREFIX: &str = "__selection:";

/// Per-atom information, modeled after PyMOL's AtomInfoType.
#[derive(Debug, Clone)]
#[allow(dead_code)]
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
    /// Segment identifier (PDB columns 73-76), used by PyMOL's segi selector.
    pub segi: String,
    /// Alternate location indicator.
    pub alt: char,
    /// Secondary structure type assignment.
    pub ss_type: SSType,
    /// Isotropic temperature factor.
    pub b_factor: f32,
    /// Occupancy.
    pub occupancy: f32,
    /// Formal charge.
    pub formal_charge: i8,
    /// Partial charge.
    pub partial_charge: f32,
    /// Force-field text type.
    pub text_type: String,
    /// Force-field numeric type.
    pub numeric_type: i32,
    /// Custom atom property string.
    pub custom: String,
    /// Arbitrary named atom properties, used by PyMOL-style `p.<name>` selectors.
    pub properties: HashMap<String, String>,
    /// Atom label text.
    pub label: String,
    /// Stereochemistry label (`R`, `S`, `odd`, `even`, or `?`).
    pub stereo: String,
    /// Van der Waals radius (from element table).
    pub vdw: f32,
    /// Electrostatic radius.
    pub elec_radius: f32,
    /// Per-atom cartoon type override; 0 means automatic/default.
    pub cartoon: i8,
    /// Molecular geometry code.
    pub geom: i8,
    /// Atom valence code.
    pub valence: i8,
    /// Display color [r, g, b] in 0.0..1.0.
    pub color: [f32; 3],
    /// Explicit cartoon color override; `None` uses normal atom color/settings.
    pub cartoon_color: Option<[f32; 3]>,
    /// Explicit ribbon color override; `None` uses normal atom color/settings.
    pub ribbon_color: Option<[f32; 3]>,
    /// Bitmask of currently visible representations.
    pub vis_rep: u32,
    /// True if this atom came from a HETATM record.
    pub is_hetatm: bool,
    /// PDB serial number.
    pub serial: u32,
    /// Original 1-based atom load order, used by PyMOL's rank selector.
    pub rank: u32,
    /// PyMOL atom flags bitmask.
    pub flags: u32,
    /// True if this atom is masked from picking/selection.
    pub masked: bool,
    /// True if this atom is protected from movement.
    pub protected: bool,
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
            segi: String::new(),
            alt: ' ',
            ss_type: SSType::Loop,
            b_factor: 0.0,
            occupancy: 1.0,
            formal_charge: 0,
            partial_charge: 0.0,
            text_type: "??".to_string(),
            numeric_type: -9999,
            custom: String::new(),
            properties: HashMap::new(),
            label: String::new(),
            stereo: String::new(),
            vdw: 1.7,
            elec_radius: 0.0,
            cartoon: 0,
            geom: 0,
            valence: 0,
            color: [0.2, 1.0, 0.2], // default carbon green
            cartoon_color: None,
            ribbon_color: None,
            vis_rep: REP_LINES,
            is_hetatm: false,
            serial: 0,
            rank: 0,
            flags: 0,
            masked: false,
            protected: false,
        }
    }
}
