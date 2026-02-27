use bytemuck::{Pod, Zeroable};
use glam::Mat4;
use wgpu::util::DeviceExt;

use crate::core::atom::REP_LINES;
use crate::core::molecule::Molecule;
use crate::render::camera::Camera;
use crate::render::rep_spheres::SphereRep;
use crate::render::rep_sticks::StickRep;

// ── Unified GPU uniform struct (matches all three shaders) ──────────────

#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable)]
pub struct Uniforms {
    pub view_proj: [[f32; 4]; 4],
    pub view: [[f32; 4]; 4],
    pub proj: [[f32; 4]; 4],
    pub eye_pos: [f32; 4],
    pub light_dir: [f32; 4],
    pub viewport_size: [f32; 4],
}

#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable)]
pub struct LineVertex {
    pub position: [f32; 3],
    pub color: [f32; 3],
}

// ── Main renderer ───────────────────────────────────────────────────────

pub struct MolRenderer {
    // Shared
    uniform_buffer: wgpu::Buffer,
    uniform_bind_group: wgpu::BindGroup,
    depth_texture: Option<wgpu::TextureView>,
    depth_size: (u32, u32),

    // Lines
    line_pipeline: wgpu::RenderPipeline,
    line_vertex_buffer: Option<wgpu::Buffer>,
    line_vertex_count: u32,

    // Spheres
    sphere_pipeline: wgpu::RenderPipeline,
    sphere_rep: SphereRep,

    // Sticks
    stick_pipeline: wgpu::RenderPipeline,
    stick_rep: StickRep,
}

impl MolRenderer {
    pub fn new(device: &wgpu::Device, target_format: wgpu::TextureFormat) -> Self {
        // ── Shared uniform buffer + bind group ──────────────────────────
        let uniform_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("uniform_buffer"),
            contents: bytemuck::cast_slice(&[Uniforms {
                view_proj: Mat4::IDENTITY.to_cols_array_2d(),
                view: Mat4::IDENTITY.to_cols_array_2d(),
                proj: Mat4::IDENTITY.to_cols_array_2d(),
                eye_pos: [0.0; 4],
                light_dir: [0.3, -0.8, -0.5, 0.0],
                viewport_size: [1280.0, 800.0, 0.0, 0.0],
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
            label: Some("shared_pipeline_layout"),
            bind_group_layouts: &[&bind_group_layout],
            push_constant_ranges: &[],
        });

        let depth_stencil = wgpu::DepthStencilState {
            format: wgpu::TextureFormat::Depth32Float,
            depth_write_enabled: true,
            depth_compare: wgpu::CompareFunction::Less,
            stencil: wgpu::StencilState::default(),
            bias: wgpu::DepthBiasState::default(),
        };

        // ── Line pipeline ───────────────────────────────────────────────
        let line_shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("line_shader"),
            source: wgpu::ShaderSource::Wgsl(include_str!("shaders/line.wgsl").into()),
        });

        let line_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("line_pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &line_shader,
                entry_point: Some("vs_main"),
                buffers: &[wgpu::VertexBufferLayout {
                    array_stride: std::mem::size_of::<LineVertex>() as u64,
                    step_mode: wgpu::VertexStepMode::Vertex,
                    attributes: &[
                        wgpu::VertexAttribute { offset: 0, shader_location: 0, format: wgpu::VertexFormat::Float32x3 },
                        wgpu::VertexAttribute { offset: 12, shader_location: 1, format: wgpu::VertexFormat::Float32x3 },
                    ],
                }],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &line_shader,
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
                cull_mode: None,
                ..Default::default()
            },
            depth_stencil: Some(depth_stencil.clone()),
            multisample: wgpu::MultisampleState::default(),
            multiview: None,
            cache: None,
        });

        // ── Sphere pipeline ─────────────────────────────────────────────
        let sphere_shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("sphere_shader"),
            source: wgpu::ShaderSource::Wgsl(include_str!("shaders/sphere.wgsl").into()),
        });

        use crate::render::rep_spheres::SphereInstance;
        let sphere_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("sphere_pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &sphere_shader,
                entry_point: Some("vs_main"),
                buffers: &[wgpu::VertexBufferLayout {
                    array_stride: std::mem::size_of::<SphereInstance>() as u64,
                    step_mode: wgpu::VertexStepMode::Instance,
                    attributes: &[
                        wgpu::VertexAttribute { offset: 0, shader_location: 0, format: wgpu::VertexFormat::Float32x3 }, // center
                        wgpu::VertexAttribute { offset: 12, shader_location: 1, format: wgpu::VertexFormat::Float32 },   // radius
                        wgpu::VertexAttribute { offset: 16, shader_location: 2, format: wgpu::VertexFormat::Float32x3 }, // color
                    ],
                }],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &sphere_shader,
                entry_point: Some("fs_main"),
                targets: &[Some(wgpu::ColorTargetState {
                    format: target_format,
                    blend: Some(wgpu::BlendState::REPLACE),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
                compilation_options: Default::default(),
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                cull_mode: None,
                ..Default::default()
            },
            depth_stencil: Some(depth_stencil.clone()),
            multisample: wgpu::MultisampleState::default(),
            multiview: None,
            cache: None,
        });

        // ── Stick pipeline ──────────────────────────────────────────────
        let stick_shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("cylinder_shader"),
            source: wgpu::ShaderSource::Wgsl(include_str!("shaders/cylinder.wgsl").into()),
        });

        use crate::render::rep_sticks::CylinderVertex;
        let stick_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("stick_pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &stick_shader,
                entry_point: Some("vs_main"),
                buffers: &[wgpu::VertexBufferLayout {
                    array_stride: std::mem::size_of::<CylinderVertex>() as u64,
                    step_mode: wgpu::VertexStepMode::Vertex,
                    attributes: &[
                        wgpu::VertexAttribute { offset: 0, shader_location: 0, format: wgpu::VertexFormat::Float32x3 },  // position
                        wgpu::VertexAttribute { offset: 12, shader_location: 1, format: wgpu::VertexFormat::Float32x3 }, // normal
                        wgpu::VertexAttribute { offset: 24, shader_location: 2, format: wgpu::VertexFormat::Float32x3 }, // color
                    ],
                }],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &stick_shader,
                entry_point: Some("fs_main"),
                targets: &[Some(wgpu::ColorTargetState {
                    format: target_format,
                    blend: Some(wgpu::BlendState::REPLACE),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
                compilation_options: Default::default(),
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                cull_mode: Some(wgpu::Face::Back),
                ..Default::default()
            },
            depth_stencil: Some(depth_stencil),
            multisample: wgpu::MultisampleState::default(),
            multiview: None,
            cache: None,
        });

        Self {
            uniform_buffer,
            uniform_bind_group,
            depth_texture: None,
            depth_size: (0, 0),
            line_pipeline,
            line_vertex_buffer: None,
            line_vertex_count: 0,
            sphere_pipeline,
            sphere_rep: SphereRep::new(),
            stick_pipeline,
            stick_rep: StickRep::new(),
        }
    }

    /// Rebuild all representation geometry buffers.
    pub fn update_geometry(&mut self, device: &wgpu::Device, molecules: &[Molecule]) {
        // Lines
        let mut line_verts: Vec<LineVertex> = Vec::new();
        for mol in molecules {
            if !mol.visible { continue; }
            for bond in &mol.bonds {
                let a = &mol.atoms[bond.atom_a];
                let b = &mol.atoms[bond.atom_b];
                if (a.vis_rep & REP_LINES) == 0 || (b.vis_rep & REP_LINES) == 0 { continue; }
                let pa = mol.coords[bond.atom_a];
                let pb = mol.coords[bond.atom_b];
                let mid = [(pa[0]+pb[0])*0.5, (pa[1]+pb[1])*0.5, (pa[2]+pb[2])*0.5];
                line_verts.push(LineVertex { position: pa, color: a.color });
                line_verts.push(LineVertex { position: mid, color: a.color });
                line_verts.push(LineVertex { position: mid, color: b.color });
                line_verts.push(LineVertex { position: pb, color: b.color });
            }
        }
        self.line_vertex_count = line_verts.len() as u32;
        self.line_vertex_buffer = if line_verts.is_empty() {
            None
        } else {
            Some(device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("line_vertices"),
                contents: bytemuck::cast_slice(&line_verts),
                usage: wgpu::BufferUsages::VERTEX,
            }))
        };

        // Spheres
        self.sphere_rep.update(device, molecules);

        // Sticks
        self.stick_rep.update(device, molecules);
    }

    /// Update the uniform buffer with current camera matrices.
    pub fn update_uniforms(&self, queue: &wgpu::Queue, camera: &Camera, aspect: f32, width: u32, height: u32) {
        let view = camera.view_matrix();
        let proj = camera.projection_matrix(aspect);
        let view_proj = proj * view;
        let eye = camera.eye_position();

        let uniforms = Uniforms {
            view_proj: view_proj.to_cols_array_2d(),
            view: view.to_cols_array_2d(),
            proj: proj.to_cols_array_2d(),
            eye_pos: [eye.x, eye.y, eye.z, 1.0],
            light_dir: [0.3, -0.8, -0.5, 0.0],
            viewport_size: [width as f32, height as f32, 0.0, 0.0],
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

    /// Render all representations.
    pub fn paint(
        &self,
        encoder: &mut wgpu::CommandEncoder,
        color_view: &wgpu::TextureView,
    ) {
        let depth_view = match &self.depth_texture {
            Some(v) => v,
            None => return,
        };

        let has_anything = self.line_vertex_count > 0
            || self.sphere_rep.instance_count > 0
            || self.stick_rep.index_count > 0;

        if !has_anything {
            // Still clear the framebuffer
            let _rpass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("clear_pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: color_view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color { r: 0.0, g: 0.0, b: 0.0, a: 1.0 }),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                timestamp_writes: None,
                occlusion_query_set: None,
            });
            return;
        }

        let mut rpass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
            label: Some("mol_render_pass"),
            color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                view: color_view,
                resolve_target: None,
                ops: wgpu::Operations {
                    load: wgpu::LoadOp::Clear(wgpu::Color { r: 0.0, g: 0.0, b: 0.0, a: 1.0 }),
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

        // Draw lines
        if self.line_vertex_count > 0 {
            if let Some(vb) = &self.line_vertex_buffer {
                rpass.set_pipeline(&self.line_pipeline);
                rpass.set_bind_group(0, &self.uniform_bind_group, &[]);
                rpass.set_vertex_buffer(0, vb.slice(..));
                rpass.draw(0..self.line_vertex_count, 0..1);
            }
        }

        // Draw sticks
        if self.stick_rep.index_count > 0 {
            if let (Some(vb), Some(ib)) = (&self.stick_rep.vertex_buffer, &self.stick_rep.index_buffer) {
                rpass.set_pipeline(&self.stick_pipeline);
                rpass.set_bind_group(0, &self.uniform_bind_group, &[]);
                rpass.set_vertex_buffer(0, vb.slice(..));
                rpass.set_index_buffer(ib.slice(..), wgpu::IndexFormat::Uint32);
                rpass.draw_indexed(0..self.stick_rep.index_count, 0, 0..1);
            }
        }

        // Draw spheres (instanced, 6 verts per instance = 1 billboard quad)
        if self.sphere_rep.instance_count > 0 {
            if let Some(ib) = &self.sphere_rep.instance_buffer {
                rpass.set_pipeline(&self.sphere_pipeline);
                rpass.set_bind_group(0, &self.uniform_bind_group, &[]);
                rpass.set_vertex_buffer(0, ib.slice(..));
                rpass.draw(0..6, 0..self.sphere_rep.instance_count);
            }
        }
    }
}
