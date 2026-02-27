use std::path::PathBuf;

use eframe::egui_wgpu;
use egui::{CentralPanel, SidePanel, TopBottomPanel};

use crate::core::atom::{REP_LINES, REP_STICKS, REP_SPHERES, REP_CARTOON};
use crate::io::{cif, pdb};
use crate::render::renderer::MolRenderer;
use crate::scene::scene::{Measurement, Scene};
use crate::selection::{parse_selection, evaluate, evaluator::count_selected};
use crate::ui::command_line::CommandLine;
use crate::ui::control_panel::{self, ControlPanelState};
use crate::ui::object_panel;

pub struct MolApp {
    pub scene: Scene,
    pub command_line: CommandLine,
    pub control_panel_state: ControlPanelState,
    /// wgpu render state obtained from creation context.
    render_state: Option<egui_wgpu::RenderState>,
    /// The renderer is created lazily on first frame (needs wgpu device).
    renderer: Option<MolRenderer>,
    /// Offscreen color texture + egui texture id for displaying the 3D viewport.
    offscreen: Option<OffscreenTarget>,
    /// Track if we need to open a file dialog.
    open_file_requested: bool,
    /// Pending file to load (from file dialog).
    pending_file: Option<PathBuf>,
    /// Pending screenshot path.
    screenshot_requested: Option<PathBuf>,
}

struct OffscreenTarget {
    #[allow(dead_code)]
    texture: wgpu::Texture,
    view: wgpu::TextureView,
    egui_tex_id: egui::TextureId,
    width: u32,
    height: u32,
}

impl MolApp {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        let render_state = cc.wgpu_render_state.clone();
        Self {
            scene: Scene::default(),
            command_line: CommandLine::default(),
            control_panel_state: ControlPanelState::default(),
            render_state,
            renderer: None,
            offscreen: None,
            open_file_requested: false,
            pending_file: None,
            screenshot_requested: None,
        }
    }

    fn ensure_renderer(&mut self) {
        if self.renderer.is_none() {
            if let Some(rs) = &self.render_state {
                let format = wgpu::TextureFormat::Bgra8UnormSrgb;
                self.renderer = Some(MolRenderer::new(&rs.device, format));
            }
        }
    }

    fn ensure_offscreen(
        &mut self,
        device: &wgpu::Device,
        egui_renderer: &mut egui_wgpu::Renderer,
        width: u32,
        height: u32,
    ) {
        let w = width.max(1);
        let h = height.max(1);

        let needs_recreate = match &self.offscreen {
            Some(o) => o.width != w || o.height != h,
            None => true,
        };

        if !needs_recreate {
            return;
        }

        // Clean up old texture registration
        if let Some(old) = self.offscreen.take() {
            egui_renderer.free_texture(&old.egui_tex_id);
        }

        let texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("offscreen_color"),
            size: wgpu::Extent3d {
                width: w,
                height: h,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Bgra8UnormSrgb,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::TEXTURE_BINDING,
            view_formats: &[],
        });

        let view = texture.create_view(&wgpu::TextureViewDescriptor::default());

        let egui_tex_id = egui_renderer.register_native_texture(
            device,
            &view,
            wgpu::FilterMode::Linear,
        );

        self.offscreen = Some(OffscreenTarget {
            texture,
            view,
            egui_tex_id,
            width: w,
            height: h,
        });
    }

    fn load_file(&mut self, path: PathBuf) {
        let ext = path.extension().and_then(|s| s.to_str()).unwrap_or("").to_lowercase();
        let res = if ext == "cif" || ext == "mmcif" {
            cif::load_cif(&path)
        } else {
            pdb::load_pdb(&path)
        };

        match res {
            Ok(mol) => {
                let name = mol.name.clone();
                let atoms = mol.atoms.len();
                let bonds = mol.bonds.len();
                self.scene.add_molecule(mol);
                self.command_line
                    .log(format!("Loaded '{}': {} atoms, {} bonds", name, atoms, bonds));
            }
            Err(e) => {
                self.command_line.log(format!("Error: {}", e));
            }
        }
    }

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

    fn handle_command(&mut self, cmd: &str) {
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
                        self.command_line.log(format!("Unknown rep: '{}'. Use lines/sticks/spheres/cartoon", rep_name));
                        return;
                    }
                };
                let sel = match parse_selection(&sel_str) {
                    Ok(s) => s,
                    Err(e) => { self.command_line.log(format!("Selection error: {}", e)); return; }
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
                self.command_line.log(format!("show {}: {} atoms", rep_name, total));
            }
            "hide" => {
                let (rep_name, sel_str) = Self::parse_rep_selection(args);
                let flag = match Self::rep_flag(&rep_name) {
                    Some(f) => f,
                    None => {
                        self.command_line.log(format!("Unknown rep: '{}'. Use lines/sticks/spheres/cartoon", rep_name));
                        return;
                    }
                };
                let sel = match parse_selection(&sel_str) {
                    Ok(s) => s,
                    Err(e) => { self.command_line.log(format!("Selection error: {}", e)); return; }
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
                self.command_line.log(format!("hide {}: {} atoms", rep_name, total));
            }
            "color" => {
                // color <color_name>, <selection>
                let (color_name, sel_str) = Self::parse_rep_selection(args);
                let rgb = match Self::parse_color(&color_name) {
                    Some(c) => c,
                    None => {
                        self.command_line.log(format!("Unknown color: '{}'. Try: red green blue yellow cyan magenta orange white gray pink salmon purple", color_name));
                        return;
                    }
                };
                let sel = match parse_selection(&sel_str) {
                    Ok(s) => s,
                    Err(e) => { self.command_line.log(format!("Selection error: {}", e)); return; }
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
                self.command_line.log(format!("color {}: {} atoms", color_name, total));
            }
            "select" => {
                // select <selection> — just counts matching atoms
                let sel = match parse_selection(args) {
                    Ok(s) => s,
                    Err(e) => { self.command_line.log(format!("Selection error: {}", e)); return; }
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
                    Err(e) => { self.command_line.log(format!("Selection 1 error: {}", e)); return; }
                };
                let sel2 = match parse_selection(sel2_str) {
                    Ok(s) => s,
                    Err(e) => { self.command_line.log(format!("Selection 2 error: {}", e)); return; }
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
                            c1[0] += p[0]; c1[1] += p[1]; c1[2] += p[2];
                            n1 += 1.0;
                        }
                        if mask2[i] {
                            c2[0] += p[0]; c2[1] += p[1]; c2[2] += p[2];
                            n2 += 1.0;
                        }
                    }
                }

                if n1 == 0.0 || n2 == 0.0 {
                    self.command_line.log("One or both selections are empty.");
                    return;
                }

                let p1 = [c1[0]/n1, c1[1]/n1, c1[2]/n1];
                let p2 = [c2[0]/n2, c2[1]/n2, c2[2]/n2];
                let dist = ((p1[0]-p2[0]).powi(2) + (p1[1]-p2[1]).powi(2) + (p1[2]-p2[2]).powi(2)).sqrt();

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
                    self.command_line.log(format!("Background set to {}", color_name));
                } else {
                    self.command_line.log(format!("Unknown color: '{}'", color_name));
                }
            }
            "help" => {
                self.command_line.log("Commands:");
                self.command_line.log("  load <file>             — Load a PDB file");
                self.command_line.log("  show <rep>[, <sel>]     — Show representation (lines/sticks/spheres/cartoon)");
                self.command_line.log("  hide <rep>[, <sel>]     — Hide representation");
                self.command_line.log("  color <color>[, <sel>]  — Color atoms");
                self.command_line.log("  select <sel>            — Count matching atoms");
                self.command_line.log("  distance <s1>, <s2>     — Measure distance");
                self.command_line.log("  png <file>              — Save screenshot");
                self.command_line.log("  bg_color <color>        — Set background color");
                self.command_line.log("  reset                   — Reset camera view");
                self.command_line.log("Selections: chain A, resi 1-50, name CA, resn ALA, elem C, hetatm, all, not/and/or, ()");
            }
            _ => {
                self.command_line
                    .log(format!("Unknown command: '{}'. Type 'help' for usage.", verb));
            }
        }
    }
}

fn save_screenshot(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    texture: &wgpu::Texture,
    path: &std::path::Path,
) {
    let size = texture.size();
    let width = size.width;
    let height = size.height;
    // Align to 256 bytes
    let bytes_per_pixel = 4;
    let unpadded_bytes_per_row = width * bytes_per_pixel;
    let align = 256;
    let padded_bytes_per_row = (unpadded_bytes_per_row + align - 1) / align * align;

    let buffer_size = (padded_bytes_per_row * height) as u64;

    let buffer = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("screenshot_buffer"),
        size: buffer_size,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });

    let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor::default());
    encoder.copy_texture_to_buffer(
        wgpu::ImageCopyTexture {
            texture,
            mip_level: 0,
            origin: wgpu::Origin3d::ZERO,
            aspect: wgpu::TextureAspect::All,
        },
        wgpu::ImageCopyBuffer {
            buffer: &buffer,
            layout: wgpu::ImageDataLayout {
                offset: 0,
                bytes_per_row: Some(padded_bytes_per_row),
                rows_per_image: Some(height),
            },
        },
        size,
    );
    queue.submit(Some(encoder.finish()));

    let slice = buffer.slice(..);
    let (tx, rx) = std::sync::mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |res| {
        tx.send(res).unwrap();
    });
    device.poll(wgpu::Maintain::Wait);
    rx.recv().unwrap().unwrap();

    let data = slice.get_mapped_range();
    // Convert BGRA (wgpu) to RGBA (image) and remove padding
    let mut pixels: Vec<u8> = Vec::with_capacity((width * height * 4) as usize);
    for row in 0..height {
        let start = (row * padded_bytes_per_row) as usize;
        let end = start + unpadded_bytes_per_row as usize;
        let row_data = &data[start..end];
        for chunk in row_data.chunks(4) {
            // BGRA -> RGBA
            pixels.push(chunk[2]);
            pixels.push(chunk[1]);
            pixels.push(chunk[0]);
            pixels.push(chunk[3]);
        }
    }
    
    if let Err(e) = image::save_buffer(
        path,
        &pixels,
        width,
        height,
        image::ColorType::Rgba8,
    ) {
        eprintln!("Failed to save screenshot: {}", e);
    }
    
    drop(data);
    buffer.unmap();
}

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
                        self.command_line.log(format!("Saved screenshot to {}", path.display()));
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

                if let (Some(renderer), Some(offscreen)) =
                    (&mut self.renderer, &self.offscreen)
                {
                    // Rebuild geometry if dirty
                    if self.scene.geometry_dirty {
                        renderer.update_geometry(device, &self.scene.molecules);
                        self.scene.geometry_dirty = false;
                    }

                    // Ensure depth texture matches
                    renderer.ensure_depth_texture(device, offscreen.width, offscreen.height);

                    // Update camera uniforms
                    let aspect = offscreen.width as f32 / offscreen.height as f32;
                    renderer.update_uniforms(queue, &self.scene.camera, aspect, offscreen.width, offscreen.height);

                    // Render to offscreen texture
                    let mut encoder =
                        device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
                            label: Some("mol_encoder"),
                        });

                    renderer.paint(
                        &mut encoder,
                        &offscreen.view,
                    );

                    queue.submit(std::iter::once(encoder.finish()));

                    // Display the offscreen texture as an egui image
                    let rect = response.rect;
                    let uv = egui::Rect::from_min_max(
                        egui::pos2(0.0, 0.0),
                        egui::pos2(1.0, 1.0),
                    );
                    ui.painter().image(
                        offscreen.egui_tex_id,
                        rect,
                        uv,
                        egui::Color32::WHITE,
                    );

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

                                ui.painter().line_segment([s1, s2], egui::Stroke::new(1.0, egui::Color32::YELLOW));
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
