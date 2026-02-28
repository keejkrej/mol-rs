use crate::core::molecule::Molecule;
use crate::render::camera::Camera;
use crate::scene::color::{ColorScheme, apply_color_scheme};

#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct Measurement {
    pub p1: [f32; 3],
    pub p2: [f32; 3],
    pub distance: f32,
    pub label: String,
}

pub struct Scene {
    pub molecules: Vec<Molecule>,
    pub measurements: Vec<Measurement>,
    pub camera: Camera,
    pub color_scheme: ColorScheme,
    pub current_state: usize,
    pub all_states: bool,
    /// True when geometry buffers need rebuilding.
    pub geometry_dirty: bool,
    /// Background color [r, g, b].
    pub bg_color: [f32; 3],
}

impl Default for Scene {
    fn default() -> Self {
        Self {
            molecules: Vec::new(),
            measurements: Vec::new(),
            camera: Camera::default(),
            color_scheme: ColorScheme::ByElement,
            current_state: 1,
            all_states: false,
            geometry_dirty: false,
            bg_color: [0.0, 0.0, 0.0],
        }
    }
}

impl Scene {
    pub fn max_state_count(&self) -> usize {
        self.molecules
            .iter()
            .map(|m| m.state_count())
            .max()
            .unwrap_or(1)
            .max(1)
    }

    pub fn set_state_clamped(&mut self, state: usize) {
        let max_state = self.max_state_count();
        let clamped = state.clamp(1, max_state);
        if self.current_state != clamped {
            self.current_state = clamped;
            self.geometry_dirty = true;
        }
    }

    pub fn next_state(&mut self) {
        self.set_state_clamped(self.current_state + 1);
    }

    pub fn prev_state(&mut self) {
        self.set_state_clamped(self.current_state.saturating_sub(1));
    }

    /// Add a molecule to the scene and adjust camera to fit.
    pub fn add_molecule(&mut self, mut mol: Molecule) {
        let requested_state = self.current_state;
        apply_color_scheme(&mut mol, self.color_scheme);
        let center = mol.centroid_for_state(requested_state);
        let radius = mol.radius_for_state(requested_state);
        self.molecules.push(mol);
        self.set_state_clamped(requested_state);

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
