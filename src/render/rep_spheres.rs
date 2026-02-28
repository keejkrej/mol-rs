use bytemuck::{Pod, Zeroable};
use wgpu::util::DeviceExt;

use crate::core::atom::REP_SPHERES;
use crate::core::molecule::Molecule;

/// Per-instance data for sphere impostors.
#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable)]
pub struct SphereInstance {
    pub center: [f32; 3],
    pub radius: f32,
    pub color: [f32; 3],
    pub _pad: f32,
}

/// Holds the GPU resources for sphere rendering.
pub struct SphereRep {
    pub instance_buffer: Option<wgpu::Buffer>,
    pub instance_count: u32,
}

impl SphereRep {
    pub fn new() -> Self {
        Self {
            instance_buffer: None,
            instance_count: 0,
        }
    }

    /// Rebuild the instance buffer from molecules.
    pub fn update(
        &mut self,
        device: &wgpu::Device,
        molecules: &[Molecule],
        current_state: usize,
        all_states: bool,
    ) {
        let mut instances: Vec<SphereInstance> = Vec::new();

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
                for (i, atom) in mol.atoms.iter().enumerate() {
                    if (atom.vis_rep & REP_SPHERES) == 0 {
                        continue;
                    }
                    let pos = coords[i];
                    instances.push(SphereInstance {
                        center: pos,
                        radius: atom.vdw * 0.25, // Scale down for ball-and-stick style
                        color: atom.color,
                        _pad: 0.0,
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
            label: Some("sphere_instances"),
            contents: bytemuck::cast_slice(&instances),
            usage: wgpu::BufferUsages::VERTEX,
        }));
    }
}
