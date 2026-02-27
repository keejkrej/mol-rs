use std::path::PathBuf;

use eframe::egui_wgpu;
use egui::{CentralPanel, SidePanel, TopBottomPanel};

use crate::core::atom::{REP_LINES, REP_STICKS, REP_SPHERES, REP_CARTOON};
use crate::io::pdb;
use crate::render::renderer::MolRenderer;
use crate::scene::scene::Scene;
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
        match pdb::load_pdb(&path) {
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
            "reset" => {
                if let Some(mol) = self.scene.molecules.first() {
                    let c = mol.centroid();
                    let r = mol.radius();
                    self.scene.camera.reset_to_fit(c, r);
                }
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
                .add_filter("PDB files", &["pdb", "ent"])
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
                }
            }

            // Request continuous repaint for smooth interaction
            ctx.request_repaint();
        });
    }
}
