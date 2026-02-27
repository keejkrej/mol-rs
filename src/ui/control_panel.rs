use egui::Ui;

use crate::scene::color::ColorScheme;
use crate::scene::scene::Scene;
use crate::selection::{parse_selection, evaluate, evaluator::count_selected};

/// Persistent state for the control panel.
pub struct ControlPanelState {
    pub selection_input: String,
    pub selection_count: Option<usize>,
    pub selection_error: Option<String>,
}

impl Default for ControlPanelState {
    fn default() -> Self {
        Self {
            selection_input: String::new(),
            selection_count: None,
            selection_error: None,
        }
    }
}

/// Draw the control panel.
pub fn control_panel(ui: &mut Ui, scene: &mut Scene, state: &mut ControlPanelState) {
    ui.heading("Display");
    ui.separator();

    // Color scheme
    ui.label("Color scheme:");
    let mut scheme = scene.color_scheme;
    ui.horizontal(|ui| {
        ui.radio_value(&mut scheme, ColorScheme::ByElement, "Element");
        ui.radio_value(&mut scheme, ColorScheme::ByChain, "Chain");
    });
    if scheme != scene.color_scheme {
        scene.set_color_scheme(scheme);
    }

    ui.add_space(8.0);

    // Background color
    ui.label("Background:");
    let mut bg = scene.bg_color;
    if ui.color_edit_button_rgb(&mut bg).changed() {
        scene.bg_color = bg;
    }

    ui.add_space(8.0);

    // Selection input
    ui.label("Selection:");
    let sel_response = ui.text_edit_singleline(&mut state.selection_input);
    if sel_response.changed() {
        // Live-update selection count
        let input = state.selection_input.trim();
        if input.is_empty() {
            state.selection_count = None;
            state.selection_error = None;
        } else {
            match parse_selection(input) {
                Ok(sel) => {
                    let mut total = 0;
                    for mol in &scene.molecules {
                        let mask = evaluate(&sel, mol);
                        total += count_selected(&mask);
                    }
                    state.selection_count = Some(total);
                    state.selection_error = None;
                }
                Err(e) => {
                    state.selection_count = None;
                    state.selection_error = Some(e);
                }
            }
        }
    }
    if let Some(count) = state.selection_count {
        ui.label(format!("{} atoms selected", count));
    }
    if let Some(err) = &state.selection_error {
        ui.colored_label(egui::Color32::from_rgb(255, 100, 100), err);
    }

    ui.add_space(8.0);

    // Reset view
    if ui.button("Reset View").clicked() {
        if let Some(mol) = scene.molecules.first() {
            let c = mol.centroid();
            let r = mol.radius();
            scene.camera.reset_to_fit(c, r);
        }
    }

    ui.add_space(8.0);

    // Info
    if let Some(mol) = scene.molecules.first() {
        ui.label(format!("Atoms: {}", mol.atoms.len()));
        ui.label(format!("Bonds: {}", mol.bonds.len()));
        ui.label(format!("Residues: {}", mol.residues.len()));
    }
}
