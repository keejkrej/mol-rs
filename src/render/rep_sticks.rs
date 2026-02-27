use bytemuck::{Pod, Zeroable};
use glam::Vec3;
use wgpu::util::DeviceExt;

use crate::core::atom::REP_STICKS;
use crate::core::molecule::Molecule;

const CYLINDER_SEGMENTS: usize = 8;
const STICK_RADIUS: f32 = 0.15;

/// Vertex for cylinder mesh.
#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable)]
pub struct CylinderVertex {
    pub position: [f32; 3],
    pub normal: [f32; 3],
    pub color: [f32; 3],
}

/// Holds GPU resources for stick rendering.
pub struct StickRep {
    pub vertex_buffer: Option<wgpu::Buffer>,
    pub index_buffer: Option<wgpu::Buffer>,
    pub index_count: u32,
}

impl StickRep {
    pub fn new() -> Self {
        Self {
            vertex_buffer: None,
            index_buffer: None,
            index_count: 0,
        }
    }

    /// Rebuild geometry from molecules.
    pub fn update(&mut self, device: &wgpu::Device, molecules: &[Molecule]) {
        let mut vertices: Vec<CylinderVertex> = Vec::new();
        let mut indices: Vec<u32> = Vec::new();

        for mol in molecules {
            if !mol.visible {
                continue;
            }
            for bond in &mol.bonds {
                let a = &mol.atoms[bond.atom_a];
                let b = &mol.atoms[bond.atom_b];

                if (a.vis_rep & REP_STICKS) == 0 || (b.vis_rep & REP_STICKS) == 0 {
                    continue;
                }

                let pa = Vec3::from(mol.coords[bond.atom_a]);
                let pb = Vec3::from(mol.coords[bond.atom_b]);
                let mid = (pa + pb) * 0.5;

                // First half: atom A color
                build_cylinder_segment(
                    &pa, &mid, STICK_RADIUS, a.color,
                    &mut vertices, &mut indices,
                );

                // Second half: atom B color
                build_cylinder_segment(
                    &mid, &pb, STICK_RADIUS, b.color,
                    &mut vertices, &mut indices,
                );
            }
        }

        self.index_count = indices.len() as u32;

        if vertices.is_empty() {
            self.vertex_buffer = None;
            self.index_buffer = None;
            return;
        }

        self.vertex_buffer = Some(device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("stick_vertices"),
            contents: bytemuck::cast_slice(&vertices),
            usage: wgpu::BufferUsages::VERTEX,
        }));

        self.index_buffer = Some(device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("stick_indices"),
            contents: bytemuck::cast_slice(&indices),
            usage: wgpu::BufferUsages::INDEX,
        }));
    }
}

/// Generate a cylinder mesh segment from `start` to `end`.
fn build_cylinder_segment(
    start: &Vec3,
    end: &Vec3,
    radius: f32,
    color: [f32; 3],
    vertices: &mut Vec<CylinderVertex>,
    indices: &mut Vec<u32>,
) {
    let axis = *end - *start;
    let length = axis.length();
    if length < 1e-6 {
        return;
    }
    let dir = axis / length;

    // Build an orthonormal basis
    let up = if dir.y.abs() < 0.99 {
        Vec3::Y
    } else {
        Vec3::X
    };
    let right = dir.cross(up).normalize();
    let actual_up = right.cross(dir).normalize();

    let base_idx = vertices.len() as u32;
    let seg = CYLINDER_SEGMENTS;

    // Generate ring vertices at start and end
    for ring in 0..2 {
        let center = if ring == 0 { *start } else { *end };
        for i in 0..seg {
            let angle = (i as f32) * std::f32::consts::TAU / (seg as f32);
            let cos_a = angle.cos();
            let sin_a = angle.sin();

            let normal = right * cos_a + actual_up * sin_a;
            let pos = center + normal * radius;

            vertices.push(CylinderVertex {
                position: pos.to_array(),
                normal: normal.to_array(),
                color,
            });
        }
    }

    // Generate indices for the cylinder wall (triangle strip as triangles)
    for i in 0..seg {
        let i0 = base_idx + i as u32;
        let i1 = base_idx + ((i + 1) % seg) as u32;
        let i2 = base_idx + (seg + i) as u32;
        let i3 = base_idx + (seg + (i + 1) % seg) as u32;

        // Two triangles per quad
        indices.push(i0);
        indices.push(i2);
        indices.push(i1);

        indices.push(i1);
        indices.push(i2);
        indices.push(i3);
    }

    // Cap at start
    let start_center_idx = vertices.len() as u32;
    vertices.push(CylinderVertex {
        position: start.to_array(),
        normal: (-dir).to_array(),
        color,
    });
    for i in 0..seg {
        let i0 = base_idx + i as u32;
        let i1 = base_idx + ((i + 1) % seg) as u32;
        indices.push(start_center_idx);
        indices.push(i1);
        indices.push(i0);
    }

    // Cap at end
    let end_center_idx = vertices.len() as u32;
    vertices.push(CylinderVertex {
        position: end.to_array(),
        normal: dir.to_array(),
        color,
    });
    for i in 0..seg {
        let i0 = base_idx + (seg + i) as u32;
        let i1 = base_idx + (seg + (i + 1) % seg) as u32;
        indices.push(end_center_idx);
        indices.push(i0);
        indices.push(i1);
    }
}
