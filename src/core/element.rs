/// Static element data: symbol, VDW radius, CPK color, atomic mass.
#[derive(Debug, Clone, Copy)]
pub struct ElementData {
    pub symbol: &'static str,
    pub name: &'static str,
    pub vdw: f32,
    pub mass: f32,
    pub color: [f32; 3],
}

/// Lookup element data by atomic number (protons). Returns a default for unknown elements.
pub fn element_by_number(atomic_number: u8) -> &'static ElementData {
    ELEMENTS
        .get(atomic_number as usize)
        .unwrap_or(&ELEMENTS[0])
}

/// Lookup element data by symbol string (case-insensitive first letter uppercase).
pub fn element_by_symbol(sym: &str) -> Option<&'static ElementData> {
    let sym = sym.trim();
    ELEMENTS.iter().skip(1).find(|e| {
        e.symbol.eq_ignore_ascii_case(sym)
    })
}

// CPK-ish color scheme, VDW radii from PyMOL's ElementTable (layer0/Element.cpp).
pub static ELEMENTS: &[ElementData] = &[
    // 0 - Lone pair / dummy
    ElementData { symbol: "LP", name: "Dummy", vdw: 0.5, mass: 0.0, color: [0.5, 0.5, 0.5] },
    // 1 - H
    ElementData { symbol: "H",  name: "Hydrogen",  vdw: 1.20, mass: 1.008,  color: [0.9, 0.9, 0.9] },
    // 2 - He
    ElementData { symbol: "He", name: "Helium",    vdw: 1.40, mass: 4.003,  color: [0.85, 1.0, 1.0] },
    // 3 - Li
    ElementData { symbol: "Li", name: "Lithium",   vdw: 1.82, mass: 6.941,  color: [0.8, 0.5, 1.0] },
    // 4 - Be
    ElementData { symbol: "Be", name: "Beryllium", vdw: 1.53, mass: 9.012,  color: [0.76, 1.0, 0.0] },
    // 5 - B
    ElementData { symbol: "B",  name: "Boron",     vdw: 1.92, mass: 10.81,  color: [1.0, 0.71, 0.71] },
    // 6 - C
    ElementData { symbol: "C",  name: "Carbon",    vdw: 1.70, mass: 12.011, color: [0.2, 1.0, 0.2] },
    // 7 - N
    ElementData { symbol: "N",  name: "Nitrogen",  vdw: 1.55, mass: 14.007, color: [0.2, 0.2, 1.0] },
    // 8 - O
    ElementData { symbol: "O",  name: "Oxygen",    vdw: 1.52, mass: 15.999, color: [1.0, 0.2, 0.2] },
    // 9 - F
    ElementData { symbol: "F",  name: "Fluorine",  vdw: 1.47, mass: 18.998, color: [0.56, 0.88, 0.31] },
    // 10 - Ne
    ElementData { symbol: "Ne", name: "Neon",      vdw: 1.54, mass: 20.180, color: [0.7, 0.89, 0.96] },
    // 11 - Na
    ElementData { symbol: "Na", name: "Sodium",    vdw: 2.27, mass: 22.990, color: [0.67, 0.36, 0.95] },
    // 12 - Mg
    ElementData { symbol: "Mg", name: "Magnesium", vdw: 1.73, mass: 24.305, color: [0.54, 1.0, 0.0] },
    // 13 - Al
    ElementData { symbol: "Al", name: "Aluminum",  vdw: 1.84, mass: 26.982, color: [0.75, 0.65, 0.65] },
    // 14 - Si
    ElementData { symbol: "Si", name: "Silicon",   vdw: 2.10, mass: 28.086, color: [0.94, 0.78, 0.63] },
    // 15 - P
    ElementData { symbol: "P",  name: "Phosphorus",vdw: 1.80, mass: 30.974, color: [1.0, 0.5, 0.0] },
    // 16 - S
    ElementData { symbol: "S",  name: "Sulfur",    vdw: 1.80, mass: 32.065, color: [0.9, 0.78, 0.2] },
    // 17 - Cl
    ElementData { symbol: "Cl", name: "Chlorine",  vdw: 1.75, mass: 35.453, color: [0.12, 0.94, 0.12] },
    // 18 - Ar
    ElementData { symbol: "Ar", name: "Argon",     vdw: 1.88, mass: 39.948, color: [0.5, 0.82, 0.89] },
    // 19 - K
    ElementData { symbol: "K",  name: "Potassium", vdw: 2.75, mass: 39.098, color: [0.56, 0.25, 0.83] },
    // 20 - Ca
    ElementData { symbol: "Ca", name: "Calcium",   vdw: 2.31, mass: 40.078, color: [0.24, 1.0, 0.0] },
    // 21 - Sc
    ElementData { symbol: "Sc", name: "Scandium",  vdw: 2.11, mass: 44.956, color: [0.9, 0.9, 0.9] },
    // 22 - Ti
    ElementData { symbol: "Ti", name: "Titanium",  vdw: 1.87, mass: 47.867, color: [0.75, 0.76, 0.78] },
    // 23 - V
    ElementData { symbol: "V",  name: "Vanadium",  vdw: 1.79, mass: 50.942, color: [0.65, 0.65, 0.67] },
    // 24 - Cr
    ElementData { symbol: "Cr", name: "Chromium",  vdw: 1.89, mass: 51.996, color: [0.54, 0.6, 0.78] },
    // 25 - Mn
    ElementData { symbol: "Mn", name: "Manganese", vdw: 1.97, mass: 54.938, color: [0.61, 0.48, 0.78] },
    // 26 - Fe
    ElementData { symbol: "Fe", name: "Iron",      vdw: 1.94, mass: 55.845, color: [0.88, 0.4, 0.2] },
    // 27 - Co
    ElementData { symbol: "Co", name: "Cobalt",    vdw: 1.92, mass: 58.933, color: [0.94, 0.56, 0.63] },
    // 28 - Ni
    ElementData { symbol: "Ni", name: "Nickel",    vdw: 1.63, mass: 58.693, color: [0.31, 0.82, 0.31] },
    // 29 - Cu
    ElementData { symbol: "Cu", name: "Copper",    vdw: 1.40, mass: 63.546, color: [0.78, 0.5, 0.2] },
    // 30 - Zn
    ElementData { symbol: "Zn", name: "Zinc",      vdw: 1.39, mass: 65.38,  color: [0.49, 0.5, 0.69] },
    // 31 - Ga
    ElementData { symbol: "Ga", name: "Gallium",   vdw: 1.87, mass: 69.723, color: [0.76, 0.56, 0.56] },
    // 32 - Ge
    ElementData { symbol: "Ge", name: "Germanium", vdw: 2.11, mass: 72.64,  color: [0.4, 0.56, 0.56] },
    // 33 - As
    ElementData { symbol: "As", name: "Arsenic",   vdw: 1.85, mass: 74.922, color: [0.74, 0.5, 0.89] },
    // 34 - Se
    ElementData { symbol: "Se", name: "Selenium",  vdw: 1.90, mass: 78.96,  color: [1.0, 0.63, 0.0] },
    // 35 - Br
    ElementData { symbol: "Br", name: "Bromine",   vdw: 1.85, mass: 79.904, color: [0.65, 0.16, 0.16] },
    // 36 - Kr
    ElementData { symbol: "Kr", name: "Krypton",   vdw: 2.02, mass: 83.798, color: [0.36, 0.72, 0.82] },
];
