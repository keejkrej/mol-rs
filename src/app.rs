mod commands;
mod screenshot;
mod update;

use std::path::PathBuf;

use eframe::egui_wgpu;

use crate::io::{cif, pdb};
use crate::io::LoadResult;
use crate::render::renderer::MolRenderer;
use crate::scene::scene::Scene;
use crate::ui::command_line::CommandLine;
use crate::ui::control_panel::ControlPanelState;

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
            Ok(LoadResult {
                molecule,
                warnings,
                source_model_count,
            }) => {
                let name = molecule.name.clone();
                let atoms = molecule.atoms.len();
                let bonds = molecule.bonds.len();
                let states = molecule.state_count();
                self.scene.add_molecule(molecule);
                self.command_line
                    .log(format!(
                        "Loaded '{}': {} atoms, {} bonds, {} state(s) from {} model(s)",
                        name, atoms, bonds, states, source_model_count
                    ));
                for warning in warnings {
                    self.command_line.log(format!("Warning: {}", warning));
                }
            }
            Err(e) => {
                self.command_line.log(format!("Error: {}", e));
            }
        }
    }
}
