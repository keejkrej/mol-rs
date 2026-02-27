use std::path::PathBuf;

use eframe::egui_wgpu;
use egui::{CentralPanel, SidePanel, TopBottomPanel};

use crate::io::pdb;
use crate::render::renderer::MolRenderer;
use crate::scene::scene::Scene;
use crate::ui::command_line::CommandLine;
use crate::ui::control_panel;
use crate::ui::object_panel;

pub struct MolApp {
    pub scene: Scene,
    pub command_line: CommandLine,
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

    fn handle_command(&mut self, cmd: &str) {
        let parts: Vec<&str> = cmd.split_whitespace().collect();
        if parts.is_empty() {
            return;
        }

        match parts[0] {
            "load" => {
                if parts.len() < 2 {
                    self.command_line.log("Usage: load <filename>");
                } else {
                    let path = PathBuf::from(parts[1..].join(" "));
                    self.load_file(path);
                }
            }
            "reset" => {
                if let Some(mol) = self.scene.molecules.first() {
                    let c = mol.centroid();
                    let r = mol.radius();
                    self.scene.camera.reset_to_fit(c, r);
                }
                self.command_line.log("View reset.");
            }
            "help" => {
                self.command_line.log("Commands: load <file>, reset, help");
            }
            _ => {
                self.command_line
                    .log(format!("Unknown command: '{}'", parts[0]));
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
                    control_panel::control_panel(ui, &mut self.scene);
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
                    renderer.update_uniforms(queue, &self.scene.camera, aspect);

                    // Render to offscreen texture
                    let mut encoder =
                        device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
                            label: Some("mol_encoder"),
                        });

                    renderer.paint(
                        device,
                        queue,
                        &mut encoder,
                        &offscreen.view,
                        offscreen.width,
                        offscreen.height,
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
