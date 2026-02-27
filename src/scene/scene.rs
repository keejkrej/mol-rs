use crate::core::molecule::Molecule;
use crate::render::camera::Camera;
use crate::scene::color::{ColorScheme, apply_color_scheme};

pub struct Scene {
    pub molecules: Vec<Molecule>,
    pub camera: Camera,
    pub color_scheme: ColorScheme,
    /// True when geometry buffers need rebuilding.
    pub geometry_dirty: bool,
    /// Background color [r, g, b].
    pub bg_color: [f32; 3],
}

impl Default for Scene {
    fn default() -> Self {
        Self {
            molecules: Vec::new(),
            camera: Camera::default(),
            color_scheme: ColorScheme::ByElement,
            geometry_dirty: false,
            bg_color: [0.0, 0.0, 0.0],
        }
    }
}

impl Scene {
    /// Add a molecule to the scene and adjust camera to fit.
    pub fn add_molecule(&mut self, mut mol: Molecule) {
        apply_color_scheme(&mut mol, self.color_scheme);
        let center = mol.centroid();
        let radius = mol.radius();
        self.molecules.push(mol);

        // Fit camera to the newly loaded molecule
        self.camera.reset_to_fit(center, radius);
        self.geometry_dirty = true;
    }

    /// Recolor all molecules with the given scheme.
    pub fn set_color_scheme(&mut self, scheme: ColorScheme) {
        self.color_scheme = scheme;
        for mol in &mut self.molecules {
            apply_color_scheme(mol, scheme);
        }
        self.geometry_dirty = true;
    }
}
