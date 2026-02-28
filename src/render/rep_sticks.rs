use bytemuck::{Pod, Zeroable};
use glam::Vec3;
use wgpu::util::DeviceExt;

use crate::core::atom::REP_STICKS;
use crate::core::molecule::Molecule;

const CYLINDER_SEGMENTS: usize = 12;
const STICK_RADIUS: f32 = 0.15;

/// Vertex for static cylinder mesh (model space).
#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable)]
pub struct CylinderVertex {
    pub position: [f32; 3],
    pub normal: [f32; 3],
}

/// Per-instance data.
#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable)]
pub struct StickInstance {
    pub start: [f32; 3],
    pub _pad1: f32,
    pub end: [f32; 3],
    pub _pad2: f32,
    pub color: [f32; 3],
    pub radius: f32,
}

/// Holds GPU resources for stick rendering.
pub struct StickRep {
    pub vertex_buffer: wgpu::Buffer,
    pub index_buffer: wgpu::Buffer,
    pub index_count: u32,

    pub instance_buffer: Option<wgpu::Buffer>,
    pub instance_count: u32,
}

impl StickRep {
    pub fn new(device: &wgpu::Device) -> Self {
        let (vertices, indices) = build_unit_cylinder(CYLINDER_SEGMENTS);

        let vertex_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("stick_mesh_vertices"),
            contents: bytemuck::cast_slice(&vertices),
            usage: wgpu::BufferUsages::VERTEX,
        });

        let index_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("stick_mesh_indices"),
            contents: bytemuck::cast_slice(&indices),
            usage: wgpu::BufferUsages::INDEX,
        });

        Self {
            vertex_buffer,
            index_buffer,
            index_count: indices.len() as u32,
            instance_buffer: None,
            instance_count: 0,
        }
    }

    /// Rebuild instances from molecules.
    pub fn update(
        &mut self,
        device: &wgpu::Device,
        molecules: &[Molecule],
        current_state: usize,
        all_states: bool,
    ) {
        let mut instances: Vec<StickInstance> = Vec::new();

        for mol in molecules {
            if !mol.visible {
                continue;
            }
            let (start, end) = if all_states {
                (1, mol.state_count())
            } else {
                (current_state, current_state)
            };
            for state in start..=end {
                let coords = mol.coords_for_state(state);
                for bond in &mol.bonds {
                    let a = &mol.atoms[bond.atom_a];
                    let b = &mol.atoms[bond.atom_b];

                    if (a.vis_rep & REP_STICKS) == 0 || (b.vis_rep & REP_STICKS) == 0 {
                        continue;
                    }

                    let pa = Vec3::from(coords[bond.atom_a]);
                    let pb = Vec3::from(coords[bond.atom_b]);
                    let mid = (pa + pb) * 0.5;

                    // First half: atom A color
                    instances.push(StickInstance {
                        start: pa.to_array(),
                        _pad1: 0.0,
                        end: mid.to_array(),
                        _pad2: 0.0,
                        color: a.color,
                        radius: STICK_RADIUS,
                    });

                    // Second half: atom B color
                    instances.push(StickInstance {
                        start: mid.to_array(),
                        _pad1: 0.0,
                        end: pb.to_array(),
                        _pad2: 0.0,
                        color: b.color,
                        radius: STICK_RADIUS,
                    });
                }
            }
        }

        self.instance_count = instances.len() as u32;

        if instances.is_empty() {
            self.instance_buffer = None;
            return;
        }

        self.instance_buffer = Some(device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("stick_instances"),
            contents: bytemuck::cast_slice(&instances),
            usage: wgpu::BufferUsages::VERTEX,
        }));
    }
}

/// Generate a unit cylinder mesh (radius 1, height 1, along Y).
fn build_unit_cylinder(segments: usize) -> (Vec<CylinderVertex>, Vec<u32>) {
    let mut vertices = Vec::new();
    let mut indices = Vec::new();

    let radius = 1.0;
    // We will align the cylinder along Y axis from (0,0,0) to (0,1,0).
    // Or simpler: just create it along Y, and transform later.
    // Let's create it from y=0 to y=1.

    // 1. Side rings
    // Ring at y=0
    for i in 0..segments {
        let angle = (i as f32) * std::f32::consts::TAU / (segments as f32);
        let x = angle.cos() * radius;
        let z = angle.sin() * radius;
        
        vertices.push(CylinderVertex {
            position: [x, 0.0, z],
            normal: [x, 0.0, z], // Normal points out
        });
    }
    // Ring at y=1
    for i in 0..segments {
        let angle = (i as f32) * std::f32::consts::TAU / (segments as f32);
        let x = angle.cos() * radius;
        let z = angle.sin() * radius;
        
        vertices.push(CylinderVertex {
            position: [x, 1.0, z],
            normal: [x, 0.0, z], // Normal points out
        });
    }

    // Side indices
    for i in 0..segments {
        let bottom = i as u32;
        let top = (segments + i) as u32;
        let bottom_next = ((i + 1) % segments) as u32;
        let top_next = (segments + (i + 1) % segments) as u32;

        indices.push(bottom);
        indices.push(top);
        indices.push(bottom_next);

        indices.push(bottom_next);
        indices.push(top);
        indices.push(top_next);
    }

    // 2. End caps
    // Bottom cap (y=0) - needs different vertices for flat shading normals
    let bottom_center_idx = vertices.len() as u32;
    vertices.push(CylinderVertex {
        position: [0.0, 0.0, 0.0],
        normal: [0.0, -1.0, 0.0],
    });
    let bottom_ring_start = vertices.len() as u32;
    for i in 0..segments {
        let angle = (i as f32) * std::f32::consts::TAU / (segments as f32);
        let x = angle.cos() * radius;
        let z = angle.sin() * radius;
        vertices.push(CylinderVertex {
            position: [x, 0.0, z],
            normal: [0.0, -1.0, 0.0],
        });
    }
    for i in 0..segments {
        let i0 = bottom_center_idx;
        let i1 = bottom_ring_start + i as u32;
        let i2 = bottom_ring_start + ((i + 1) % segments) as u32;
        // Clockwise for bottom looking from bottom, so CounterClockwise looking from outside?
        // Standard is CCW.
        // Bottom cap normals point -Y.
        indices.push(i0);
        indices.push(i2);
        indices.push(i1);
    }

    // Top cap (y=1)
    let top_center_idx = vertices.len() as u32;
    vertices.push(CylinderVertex {
        position: [0.0, 1.0, 0.0],
        normal: [0.0, 1.0, 0.0],
    });
    let top_ring_start = vertices.len() as u32;
    for i in 0..segments {
        let angle = (i as f32) * std::f32::consts::TAU / (segments as f32);
        let x = angle.cos() * radius;
        let z = angle.sin() * radius;
        vertices.push(CylinderVertex {
            position: [x, 1.0, z],
            normal: [0.0, 1.0, 0.0],
        });
    }
    for i in 0..segments {
        let i0 = top_center_idx;
        let i1 = top_ring_start + i as u32;
        let i2 = top_ring_start + ((i + 1) % segments) as u32;
        indices.push(i0);
        indices.push(i1);
        indices.push(i2);
    }

    (vertices, indices)
}
