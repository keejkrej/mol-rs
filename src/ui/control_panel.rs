use egui::Ui;

use crate::scene::color::ColorScheme;
use crate::scene::scene::Scene;

/// Draw the control panel.
pub fn control_panel(ui: &mut Ui, scene: &mut Scene) {
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
