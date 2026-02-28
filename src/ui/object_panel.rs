use egui::Ui;

use crate::core::atom::{REP_CARTOON, REP_LINES, REP_SPHERES, REP_STICKS};
use crate::scene::scene::Scene;

/// Draw the object list panel.
pub fn object_panel(ui: &mut Ui, scene: &mut Scene) {
    ui.heading("Objects");
    ui.separator();

    if scene.molecules.is_empty() {
        ui.label("No molecules loaded.\nUse File > Open to load a PDB/CIF.");
        return;
    }

    let mut dirty = false;

    for (_i, mol) in scene.molecules.iter_mut().enumerate() {
        ui.horizontal(|ui| {
            if ui.checkbox(&mut mol.visible, "").changed() {
                dirty = true;
            }
            ui.label(format!("{} ({} states)", mol.name, mol.state_count()));
        });

        ui.horizontal(|ui| {
            ui.label("  ");
            // Rep toggle buttons
            let reps = [
                ("L", REP_LINES),
                ("S", REP_STICKS),
                ("Sp", REP_SPHERES),
                ("C", REP_CARTOON),
            ];
            for (label, flag) in &reps {
                let any_on = mol.atoms.iter().any(|a| a.vis_rep & flag != 0);
                let mut on = any_on;
                if ui.toggle_value(&mut on, *label).changed() {
                    for atom in &mut mol.atoms {
                        if on {
                            atom.vis_rep |= flag;
                        } else {
                            atom.vis_rep &= !flag;
                        }
                    }
                    dirty = true;
                }
            }
        });

        ui.add_space(4.0);
    }

    if dirty {
        scene.geometry_dirty = true;
    }
}
