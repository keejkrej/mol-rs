use egui::{CentralPanel, SidePanel, TopBottomPanel};

use crate::ui::control_panel;
use crate::ui::object_panel;

use super::screenshot::save_screenshot;
use super::MolApp;

impl eframe::App for MolApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Handle pending file load
        if let Some(path) = self.pending_file.take() {
            self.load_file(path);
        }

        // Menu bar
        TopBottomPanel::top("menu_bar").show(ctx, |ui| {
            egui::menu::bar(ui, |ui| {
                ui.menu_button("File", |ui| {
                    if ui.button("Open PDB...").clicked() {
                        self.open_file_requested = true;
                        ui.close_menu();
                    }
                    if ui.button("Quit").clicked() {
                        ctx.send_viewport_cmd(egui::ViewportCommand::Close);
                    }
                });
            });
        });

        // File dialog (non-blocking)
        if self.open_file_requested {
            self.open_file_requested = false;
            if let Some(path) = rfd::FileDialog::new()
                .add_filter("PDB/CIF files", &["pdb", "ent", "cif"])
                .add_filter("All files", &["*"])
                .pick_file()
            {
                self.pending_file = Some(path);
            }
        }

        // Bottom panel: command line
        TopBottomPanel::bottom("command_panel")
            .resizable(true)
            .min_height(60.0)
            .show(ctx, |ui| {
                self.command_line.draw_output(ui);
                if let Some(cmd) = self.command_line.draw(ui) {
                    self.handle_command(&cmd);
                }
            });

        // Left sidebar: object list + controls
        SidePanel::left("left_panel")
            .default_width(200.0)
            .resizable(true)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    object_panel::object_panel(ui, &mut self.scene);
                    ui.add_space(16.0);
                    control_panel::control_panel(ui, &mut self.scene, &mut self.control_panel_state);
                });
            });

        // Central panel: 3D viewport
        CentralPanel::default().show(ctx, |ui| {
            let available = ui.available_size();
            let vp_width = (available.x as u32).max(1);
            let vp_height = (available.y as u32).max(1);

            // Handle mouse input for camera control
            let response = ui.allocate_rect(
                egui::Rect::from_min_size(ui.cursor().min, available),
                egui::Sense::click_and_drag(),
            );

            if response.dragged_by(egui::PointerButton::Primary) {
                let delta = response.drag_delta();
                self.scene.camera.rotate(delta.x, delta.y);
            }
            if response.dragged_by(egui::PointerButton::Middle) {
                let delta = response.drag_delta();
                self.scene.camera.pan(delta.x, delta.y);
            }

            let scroll = ui.input(|i| i.raw_scroll_delta.y);
            if scroll != 0.0 {
                self.scene.camera.zoom(scroll * 0.1);
            }

            // --- Render 3D scene to offscreen texture ---
            self.ensure_renderer();

            // Handle screenshot if requested (blocks for 1 frame)
            if let Some(path) = self.screenshot_requested.take() {
                if let Some(rs) = &self.render_state {
                    if let Some(offscreen) = &self.offscreen {
                        save_screenshot(&rs.device, &rs.queue, &offscreen.texture, &path);
                        self.command_line
                            .log(format!("Saved screenshot to {}", path.display()));
                    }
                }
            }

            // Clone Arc to avoid borrow conflicts with &mut self
            if let Some(rs) = self.render_state.clone() {
                let device = &rs.device;
                let queue = &rs.queue;

                {
                    let mut egui_renderer = rs.renderer.write();
                    self.ensure_offscreen(device, &mut egui_renderer, vp_width, vp_height);
                }

                if let (Some(renderer), Some(offscreen)) = (&mut self.renderer, &self.offscreen) {
                    // Rebuild geometry if dirty
                    if self.scene.geometry_dirty {
                        renderer.update_geometry(device, &self.scene.molecules);
                        self.scene.geometry_dirty = false;
                    }

                    // Ensure depth texture matches
                    renderer.ensure_depth_texture(device, offscreen.width, offscreen.height);

                    // Update camera uniforms
                    let aspect = offscreen.width as f32 / offscreen.height as f32;
                    renderer.update_uniforms(
                        queue,
                        &self.scene.camera,
                        aspect,
                        offscreen.width,
                        offscreen.height,
                    );

                    // Render to offscreen texture
                    let mut encoder =
                        device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
                            label: Some("mol_encoder"),
                        });

                    renderer.paint(&mut encoder, &offscreen.view);

                    queue.submit(std::iter::once(encoder.finish()));

                    // Display the offscreen texture as an egui image
                    let rect = response.rect;
                    let uv =
                        egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0));
                    ui.painter()
                        .image(offscreen.egui_tex_id, rect, uv, egui::Color32::WHITE);

                    // Draw measurements overlay
                    if !self.scene.measurements.is_empty() {
                        let view = self.scene.camera.view_matrix();
                        let proj = self.scene.camera.projection_matrix(aspect);
                        let mvp = proj * view;
                        let size = rect.size();
                        let min = rect.min;

                        for m in &self.scene.measurements {
                            let p1 = glam::Vec3::from(m.p1);
                            let p2 = glam::Vec3::from(m.p2);
                            let v1 = mvp * glam::Vec4::new(p1.x, p1.y, p1.z, 1.0);
                            let v2 = mvp * glam::Vec4::new(p2.x, p2.y, p2.z, 1.0);

                            // Simple clipping check (w > 0)
                            if v1.w > 0.0 && v2.w > 0.0 {
                                let ndc1 = v1.truncate() / v1.w;
                                let ndc2 = v2.truncate() / v2.w;

                                let s1 = egui::pos2(
                                    min.x + (ndc1.x + 1.0) * 0.5 * size.x,
                                    min.y + (1.0 - ndc1.y) * 0.5 * size.y,
                                );
                                let s2 = egui::pos2(
                                    min.x + (ndc2.x + 1.0) * 0.5 * size.x,
                                    min.y + (1.0 - ndc2.y) * 0.5 * size.y,
                                );

                                ui.painter().line_segment(
                                    [s1, s2],
                                    egui::Stroke::new(1.0, egui::Color32::YELLOW),
                                );
                                let mid = s1 + (s2 - s1) * 0.5;
                                ui.painter().text(
                                    mid,
                                    egui::Align2::CENTER_CENTER,
                                    &m.label,
                                    egui::FontId::proportional(14.0),
                                    egui::Color32::YELLOW,
                                );
                            }
                        }
                    }
                }
            }

            // Request continuous repaint for smooth interaction
            ctx.request_repaint();
        });
    }
}
