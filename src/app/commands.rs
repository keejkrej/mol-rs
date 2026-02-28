use std::path::PathBuf;

use crate::core::atom::{REP_CARTOON, REP_LINES, REP_SPHERES, REP_STICKS};
use crate::scene::scene::Measurement;
use crate::selection::{evaluate, evaluator::count_selected, parse_selection};

use super::MolApp;

impl MolApp {
    /// Parse "rep_name, selection_expr" from a comma-separated argument string.
    /// If no comma, the entire string is the rep name and selection defaults to "all".
    fn parse_rep_selection(args: &str) -> (String, String) {
        if let Some(comma) = args.find(',') {
            let rep = args[..comma].trim().to_lowercase();
            let sel = args[comma + 1..].trim().to_string();
            (rep, sel)
        } else {
            (args.trim().to_lowercase(), String::new())
        }
    }

    fn rep_flag(name: &str) -> Option<u32> {
        match name {
            "lines" | "line" => Some(REP_LINES),
            "sticks" | "stick" => Some(REP_STICKS),
            "spheres" | "sphere" => Some(REP_SPHERES),
            "cartoon" => Some(REP_CARTOON),
            _ => None,
        }
    }

    fn parse_color(name: &str) -> Option<[f32; 3]> {
        match name {
            "red" => Some([1.0, 0.2, 0.2]),
            "green" => Some([0.2, 1.0, 0.2]),
            "blue" => Some([0.2, 0.2, 1.0]),
            "yellow" => Some([1.0, 1.0, 0.2]),
            "cyan" => Some([0.2, 1.0, 1.0]),
            "magenta" => Some([1.0, 0.2, 1.0]),
            "orange" => Some([1.0, 0.6, 0.2]),
            "white" => Some([1.0, 1.0, 1.0]),
            "gray" | "grey" => Some([0.5, 0.5, 0.5]),
            "pink" => Some([1.0, 0.65, 0.85]),
            "salmon" => Some([1.0, 0.6, 0.5]),
            "purple" => Some([0.6, 0.2, 0.8]),
            _ => None,
        }
    }

    pub(super) fn handle_command(&mut self, cmd: &str) {
        let parts: Vec<&str> = cmd.splitn(2, char::is_whitespace).collect();
        let verb = parts[0].to_lowercase();
        let args = if parts.len() > 1 { parts[1].trim() } else { "" };

        match verb.as_str() {
            "load" => {
                if args.is_empty() {
                    self.command_line.log("Usage: load <filename>");
                } else {
                    let path = PathBuf::from(args);
                    self.load_file(path);
                }
            }
            "show" => {
                let (rep_name, sel_str) = Self::parse_rep_selection(args);
                let flag = match Self::rep_flag(&rep_name) {
                    Some(f) => f,
                    None => {
                        self.command_line.log(format!(
                            "Unknown rep: '{}'. Use lines/sticks/spheres/cartoon",
                            rep_name
                        ));
                        return;
                    }
                };
                let sel = match parse_selection(&sel_str) {
                    Ok(s) => s,
                    Err(e) => {
                        self.command_line.log(format!("Selection error: {}", e));
                        return;
                    }
                };
                let mut total = 0usize;
                for mol in &mut self.scene.molecules {
                    let mask = evaluate(&sel, mol);
                    for (i, atom) in mol.atoms.iter_mut().enumerate() {
                        if mask[i] {
                            atom.vis_rep |= flag;
                            total += 1;
                        }
                    }
                }
                self.scene.geometry_dirty = true;
                self.command_line
                    .log(format!("show {}: {} atoms", rep_name, total));
            }
            "hide" => {
                let (rep_name, sel_str) = Self::parse_rep_selection(args);
                let flag = match Self::rep_flag(&rep_name) {
                    Some(f) => f,
                    None => {
                        self.command_line.log(format!(
                            "Unknown rep: '{}'. Use lines/sticks/spheres/cartoon",
                            rep_name
                        ));
                        return;
                    }
                };
                let sel = match parse_selection(&sel_str) {
                    Ok(s) => s,
                    Err(e) => {
                        self.command_line.log(format!("Selection error: {}", e));
                        return;
                    }
                };
                let mut total = 0usize;
                for mol in &mut self.scene.molecules {
                    let mask = evaluate(&sel, mol);
                    for (i, atom) in mol.atoms.iter_mut().enumerate() {
                        if mask[i] {
                            atom.vis_rep &= !flag;
                            total += 1;
                        }
                    }
                }
                self.scene.geometry_dirty = true;
                self.command_line
                    .log(format!("hide {}: {} atoms", rep_name, total));
            }
            "color" => {
                // color <color_name>, <selection>
                let (color_name, sel_str) = Self::parse_rep_selection(args);
                let rgb = match Self::parse_color(&color_name) {
                    Some(c) => c,
                    None => {
                        self.command_line.log(format!(
                            "Unknown color: '{}'. Try: red green blue yellow cyan magenta orange white gray pink salmon purple",
                            color_name
                        ));
                        return;
                    }
                };
                let sel = match parse_selection(&sel_str) {
                    Ok(s) => s,
                    Err(e) => {
                        self.command_line.log(format!("Selection error: {}", e));
                        return;
                    }
                };
                let mut total = 0usize;
                for mol in &mut self.scene.molecules {
                    let mask = evaluate(&sel, mol);
                    for (i, atom) in mol.atoms.iter_mut().enumerate() {
                        if mask[i] {
                            atom.color = rgb;
                            total += 1;
                        }
                    }
                }
                self.scene.geometry_dirty = true;
                self.command_line
                    .log(format!("color {}: {} atoms", color_name, total));
            }
            "select" => {
                // select <selection> — just counts matching atoms
                let sel = match parse_selection(args) {
                    Ok(s) => s,
                    Err(e) => {
                        self.command_line.log(format!("Selection error: {}", e));
                        return;
                    }
                };
                let mut total = 0usize;
                for mol in &self.scene.molecules {
                    let mask = evaluate(&sel, mol);
                    total += count_selected(&mask);
                }
                self.command_line.log(format!("Selected {} atoms", total));
            }
            "distance" | "dist" => {
                let (sel1_str, sel2_str) = if let Some(comma) = args.find(',') {
                    (args[..comma].trim(), args[comma + 1..].trim())
                } else {
                    self.command_line.log("Usage: distance <sel1>, <sel2>");
                    return;
                };

                let sel1 = match parse_selection(sel1_str) {
                    Ok(s) => s,
                    Err(e) => {
                        self.command_line.log(format!("Selection 1 error: {}", e));
                        return;
                    }
                };
                let sel2 = match parse_selection(sel2_str) {
                    Ok(s) => s,
                    Err(e) => {
                        self.command_line.log(format!("Selection 2 error: {}", e));
                        return;
                    }
                };

                // Compute centroids
                let mut c1 = [0.0f32; 3];
                let mut n1 = 0.0f32;
                let mut c2 = [0.0f32; 3];
                let mut n2 = 0.0f32;

                for mol in &self.scene.molecules {
                    let mask1 = evaluate(&sel1, mol);
                    let mask2 = evaluate(&sel2, mol);
                    for (i, p) in mol.coords.iter().enumerate() {
                        if mask1[i] {
                            c1[0] += p[0];
                            c1[1] += p[1];
                            c1[2] += p[2];
                            n1 += 1.0;
                        }
                        if mask2[i] {
                            c2[0] += p[0];
                            c2[1] += p[1];
                            c2[2] += p[2];
                            n2 += 1.0;
                        }
                    }
                }

                if n1 == 0.0 || n2 == 0.0 {
                    self.command_line.log("One or both selections are empty.");
                    return;
                }

                let p1 = [c1[0] / n1, c1[1] / n1, c1[2] / n1];
                let p2 = [c2[0] / n2, c2[1] / n2, c2[2] / n2];
                let dist =
                    ((p1[0] - p2[0]).powi(2) + (p1[1] - p2[1]).powi(2) + (p1[2] - p2[2]).powi(2))
                        .sqrt();

                self.scene.measurements.push(Measurement {
                    p1,
                    p2,
                    distance: dist,
                    label: format!("{:.2} Å", dist),
                });
                self.command_line.log(format!("Distance: {:.2} Å", dist));
            }
            "png" => {
                if args.is_empty() {
                    self.command_line.log("Usage: png <filename>");
                } else {
                    self.screenshot_requested = Some(PathBuf::from(args));
                    self.command_line.log("Screenshot requested...");
                }
            }
            "reset" => {
                if let Some(mol) = self.scene.molecules.first() {
                    let c = mol.centroid();
                    let r = mol.radius();
                    self.scene.camera.reset_to_fit(c, r);
                }
                self.scene.measurements.clear();
                self.command_line.log("View reset.");
            }
            "bg_color" | "bg" => {
                let color_name = args.trim().to_lowercase();
                if let Some(rgb) = Self::parse_color(&color_name) {
                    self.scene.bg_color = rgb;
                    self.command_line
                        .log(format!("Background set to {}", color_name));
                } else {
                    self.command_line
                        .log(format!("Unknown color: '{}'", color_name));
                }
            }
            "help" => {
                self.command_line.log("Commands:");
                self.command_line
                    .log("  load <file>             — Load a PDB file");
                self.command_line.log(
                    "  show <rep>[, <sel>]     — Show representation (lines/sticks/spheres/cartoon)",
                );
                self.command_line
                    .log("  hide <rep>[, <sel>]     — Hide representation");
                self.command_line
                    .log("  color <color>[, <sel>]  — Color atoms");
                self.command_line
                    .log("  select <sel>            — Count matching atoms");
                self.command_line
                    .log("  distance <s1>, <s2>     — Measure distance");
                self.command_line
                    .log("  png <file>              — Save screenshot");
                self.command_line
                    .log("  bg_color <color>        — Set background color");
                self.command_line
                    .log("  reset                   — Reset camera view");
                self.command_line.log(
                    "Selections: chain A, resi 1-50, name CA, resn ALA, elem C, hetatm, all, not/and/or, ()",
                );
            }
            _ => {
                self.command_line
                    .log(format!("Unknown command: '{}'. Type 'help' for usage.", verb));
            }
        }
    }
}
