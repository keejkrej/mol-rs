use bytemuck::{Pod, Zeroable};
use glam::Mat4;
use wgpu::util::DeviceExt;

use crate::core::atom::REP_LINES;
use crate::core::molecule::Molecule;
use crate::render::camera::Camera;

// ── GPU types ───────────────────────────────────────────────────────────

#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable)]
pub struct Uniforms {
    pub view_proj: [[f32; 4]; 4],
    pub eye_pos: [f32; 4],
    pub light_dir: [f32; 4],
    pub ambient: [f32; 4],
}

#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable)]
pub struct LineVertex {
    pub position: [f32; 3],
    pub color: [f32; 3],
}

// ── Renderer resource that lives in egui_wgpu callback resources ────────

pub struct MolRenderer {
    pipeline: wgpu::RenderPipeline,
    uniform_buffer: wgpu::Buffer,
    uniform_bind_group: wgpu::BindGroup,
    vertex_buffer: Option<wgpu::Buffer>,
    num_vertices: u32,
    depth_texture: Option<wgpu::TextureView>,
    depth_size: (u32, u32),
}

impl MolRenderer {
    pub fn new(device: &wgpu::Device, target_format: wgpu::TextureFormat) -> Self {
        // Load shader
        let shader_src = include_str!("shaders/line.wgsl");
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("line_shader"),
            source: wgpu::ShaderSource::Wgsl(shader_src.into()),
        });

        // Uniform buffer + bind group layout
        let uniform_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("uniform_buffer"),
            contents: bytemuck::cast_slice(&[Uniforms {
                view_proj: Mat4::IDENTITY.to_cols_array_2d(),
                eye_pos: [0.0; 4],
                light_dir: [0.0, -1.0, -1.0, 0.0],
                ambient: [0.3, 0.3, 0.3, 1.0],
            }]),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });

        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("uniform_bgl"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            }],
        });

        let uniform_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("uniform_bg"),
            layout: &bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: uniform_buffer.as_entire_binding(),
            }],
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("line_pipeline_layout"),
            bind_group_layouts: &[&bind_group_layout],
            push_constant_ranges: &[],
        });

        let vertex_buffers = [wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<LineVertex>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &[
                wgpu::VertexAttribute {
                    offset: 0,
                    shader_location: 0,
                    format: wgpu::VertexFormat::Float32x3,
                },
                wgpu::VertexAttribute {
                    offset: 12,
                    shader_location: 1,
                    format: wgpu::VertexFormat::Float32x3,
                },
            ],
        }];

        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("line_pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: Some("vs_main"),
                buffers: &vertex_buffers,
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: Some("fs_main"),
                targets: &[Some(wgpu::ColorTargetState {
                    format: target_format,
                    blend: Some(wgpu::BlendState::REPLACE),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
                compilation_options: Default::default(),
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::LineList,
                strip_index_format: None,
                front_face: wgpu::FrontFace::Ccw,
                cull_mode: None,
                polygon_mode: wgpu::PolygonMode::Fill,
                unclipped_depth: false,
                conservative: false,
            },
            depth_stencil: Some(wgpu::DepthStencilState {
                format: wgpu::TextureFormat::Depth32Float,
                depth_write_enabled: true,
                depth_compare: wgpu::CompareFunction::Less,
                stencil: wgpu::StencilState::default(),
                bias: wgpu::DepthBiasState::default(),
            }),
            multisample: wgpu::MultisampleState::default(),
            multiview: None,
            cache: None,
        });

        Self {
            pipeline,
            uniform_buffer,
            uniform_bind_group,
            vertex_buffer: None,
            num_vertices: 0,
            depth_texture: None,
            depth_size: (0, 0),
        }
    }

    /// Rebuild the line vertex buffer from the current scene molecules.
    pub fn update_geometry(&mut self, device: &wgpu::Device, molecules: &[Molecule]) {
        let mut vertices: Vec<LineVertex> = Vec::new();

        for mol in molecules {
            if !mol.visible {
                continue;
            }
            for bond in &mol.bonds {
                let a = &mol.atoms[bond.atom_a];
                let b = &mol.atoms[bond.atom_b];

                // Only draw if both atoms have lines rep enabled
                if (a.vis_rep & REP_LINES) == 0 || (b.vis_rep & REP_LINES) == 0 {
                    continue;
                }

                let pa = mol.coords[bond.atom_a];
                let pb = mol.coords[bond.atom_b];

                // Midpoint for two-color bond
                let mid = [
                    (pa[0] + pb[0]) * 0.5,
                    (pa[1] + pb[1]) * 0.5,
                    (pa[2] + pb[2]) * 0.5,
                ];

                // First half: atom A color
                vertices.push(LineVertex { position: pa, color: a.color });
                vertices.push(LineVertex { position: mid, color: a.color });

                // Second half: atom B color
                vertices.push(LineVertex { position: mid, color: b.color });
                vertices.push(LineVertex { position: pb, color: b.color });
            }
        }

        self.num_vertices = vertices.len() as u32;

        if vertices.is_empty() {
            self.vertex_buffer = None;
            return;
        }

        self.vertex_buffer = Some(device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("line_vertices"),
            contents: bytemuck::cast_slice(&vertices),
            usage: wgpu::BufferUsages::VERTEX,
        }));
    }

    /// Update the uniform buffer with current camera matrices.
    pub fn update_uniforms(&self, queue: &wgpu::Queue, camera: &Camera, aspect: f32) {
        let view = camera.view_matrix();
        let proj = camera.projection_matrix(aspect);
        let view_proj = proj * view;
        let eye = camera.eye_position();

        let uniforms = Uniforms {
            view_proj: view_proj.to_cols_array_2d(),
            eye_pos: [eye.x, eye.y, eye.z, 1.0],
            light_dir: [0.3, -0.8, -0.5, 0.0],
            ambient: [0.3, 0.3, 0.3, 1.0],
        };

        queue.write_buffer(&self.uniform_buffer, 0, bytemuck::cast_slice(&[uniforms]));
    }

    /// Ensure the depth texture matches the given viewport size.
    pub fn ensure_depth_texture(&mut self, device: &wgpu::Device, width: u32, height: u32) {
        if self.depth_size == (width, height) && self.depth_texture.is_some() {
            return;
        }
        let w = width.max(1);
        let h = height.max(1);
        let tex = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("depth_texture"),
            size: wgpu::Extent3d { width: w, height: h, depth_or_array_layers: 1 },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Depth32Float,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            view_formats: &[],
        });
        self.depth_texture = Some(tex.create_view(&wgpu::TextureViewDescriptor::default()));
        self.depth_size = (w, h);
    }

    /// Render the molecule geometry. Called from the egui paint callback.
    pub fn paint(
        &self,
        _device: &wgpu::Device,
        _queue: &wgpu::Queue,
        encoder: &mut wgpu::CommandEncoder,
        color_view: &wgpu::TextureView,
        _width: u32,
        _height: u32,
    ) {
        let depth_view = match &self.depth_texture {
            Some(v) => v,
            None => return,
        };

        if self.num_vertices == 0 {
            return;
        }

        let mut rpass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
            label: Some("mol_render_pass"),
            color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                view: color_view,
                resolve_target: None,
                ops: wgpu::Operations {
                    load: wgpu::LoadOp::Clear(wgpu::Color {
                        r: 0.0,
                        g: 0.0,
                        b: 0.0,
                        a: 1.0,
                    }),
                    store: wgpu::StoreOp::Store,
                },
            })],
            depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                view: depth_view,
                depth_ops: Some(wgpu::Operations {
                    load: wgpu::LoadOp::Clear(1.0),
                    store: wgpu::StoreOp::Store,
                }),
                stencil_ops: None,
            }),
            timestamp_writes: None,
            occlusion_query_set: None,
        });

        rpass.set_pipeline(&self.pipeline);
        rpass.set_bind_group(0, &self.uniform_bind_group, &[]);
        if let Some(vb) = &self.vertex_buffer {
            rpass.set_vertex_buffer(0, vb.slice(..));
            rpass.draw(0..self.num_vertices, 0..1);
        }
    }
}
