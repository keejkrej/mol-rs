use bytemuck::{Pod, Zeroable};
use glam::Vec3;
use wgpu::util::DeviceExt;

use crate::core::atom::REP_CARTOON;
use crate::core::molecule::Molecule;
use crate::core::secondary_structure::SSType;

const SPLINE_SUBDIVISIONS: usize = 4; // interpolated points between each CA pair
const TUBE_SEGMENTS: usize = 8;       // radial segments for tubes
const TUBE_RADIUS: f32 = 0.25;        // loop tube radius
const HELIX_WIDTH: f32 = 1.4;         // helix ribbon half-width
const HELIX_THICKNESS: f32 = 0.25;    // helix ribbon half-thickness
const SHEET_WIDTH: f32 = 1.2;         // sheet ribbon half-width
const SHEET_THICKNESS: f32 = 0.15;    // sheet ribbon half-thickness
const ARROW_HEAD_WIDTH: f32 = 1.8;    // arrow head half-width at start of taper

/// Vertex for cartoon mesh — same layout as CylinderVertex for pipeline reuse.
#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable)]
pub struct CartoonVertex {
    pub position: [f32; 3],
    pub normal: [f32; 3],
    pub color: [f32; 3],
}

/// Holds GPU resources for cartoon rendering.
pub struct CartoonRep {
    pub vertex_buffer: Option<wgpu::Buffer>,
    pub index_buffer: Option<wgpu::Buffer>,
    pub index_count: u32,
}

/// A point along the spline with associated metadata.
struct SplinePoint {
    pos: Vec3,
    /// Tangent direction (forward along backbone).
    tangent: Vec3,
    /// Normal vector (used for ribbon orientation).
    normal: Vec3,
    /// Secondary structure type.
    ss: SSType,
    /// Color at this point.
    color: [f32; 3],
    /// Fractional position within the current SS segment (0..1), used for arrow tapering.
    ss_frac: f32,
}

impl CartoonRep {
    pub fn new() -> Self {
        Self {
            vertex_buffer: None,
            index_buffer: None,
            index_count: 0,
        }
    }

    /// Rebuild geometry from molecules.
    pub fn update(&mut self, device: &wgpu::Device, molecules: &[Molecule]) {
        let mut vertices: Vec<CartoonVertex> = Vec::new();
        let mut indices: Vec<u32> = Vec::new();

        for mol in molecules {
            if !mol.visible {
                continue;
            }
            // Extract backbone traces per chain
            let chains = extract_backbone_chains(mol);
            for chain in &chains {
                if chain.len() < 2 {
                    continue;
                }
                let spline = build_spline(chain);
                if spline.len() < 2 {
                    continue;
                }
                generate_cartoon_mesh(&spline, &mut vertices, &mut indices);
            }
        }

        self.index_count = indices.len() as u32;

        if vertices.is_empty() {
            self.vertex_buffer = None;
            self.index_buffer = None;
            return;
        }

        self.vertex_buffer = Some(device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("cartoon_vertices"),
            contents: bytemuck::cast_slice(&vertices),
            usage: wgpu::BufferUsages::VERTEX,
        }));

        self.index_buffer = Some(device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("cartoon_indices"),
            contents: bytemuck::cast_slice(&indices),
            usage: wgpu::BufferUsages::INDEX,
        }));
    }
}

/// A backbone control point extracted from the molecule.
struct ControlPoint {
    pos: Vec3,
    ss: SSType,
    color: [f32; 3],
}

/// Extract per-chain backbone CA traces from a molecule.
/// Only includes atoms with REP_CARTOON enabled.
fn extract_backbone_chains(mol: &Molecule) -> Vec<Vec<ControlPoint>> {
    let mut chains: Vec<Vec<ControlPoint>> = Vec::new();
    let mut current_chain: Vec<ControlPoint> = Vec::new();
    let mut last_chain_id: Option<char> = None;

    for res in &mol.residues {
        // Use CA atom if available
        let ca_idx = match res.ca_index {
            Some(idx) => idx,
            None => continue,
        };

        let atom = &mol.atoms[ca_idx];
        if (atom.vis_rep & REP_CARTOON) == 0 {
            continue;
        }

        // Chain break detection
        if let Some(prev_chain) = last_chain_id {
            if atom.chain != prev_chain {
                if current_chain.len() >= 2 {
                    chains.push(std::mem::take(&mut current_chain));
                } else {
                    current_chain.clear();
                }
            }
        }
        last_chain_id = Some(atom.chain);

        current_chain.push(ControlPoint {
            pos: Vec3::from(mol.coords[ca_idx]),
            ss: res.ss_type,
            color: atom.color,
        });
    }

    if current_chain.len() >= 2 {
        chains.push(current_chain);
    }

    chains
}

/// Build a smooth spline from backbone control points using Catmull-Rom interpolation.
fn build_spline(controls: &[ControlPoint]) -> Vec<SplinePoint> {
    let n = controls.len();
    if n < 2 {
        return Vec::new();
    }

    let mut spline_points: Vec<SplinePoint> = Vec::new();
    let total_segments = n - 1;

    for seg in 0..total_segments {
        let p0 = if seg > 0 { controls[seg - 1].pos } else { controls[seg].pos * 2.0 - controls[seg + 1].pos };
        let p1 = controls[seg].pos;
        let p2 = controls[seg + 1].pos;
        let p3 = if seg + 2 < n { controls[seg + 2].pos } else { controls[seg + 1].pos * 2.0 - controls[seg].pos };

        let subdivs = if seg == total_segments - 1 { SPLINE_SUBDIVISIONS + 1 } else { SPLINE_SUBDIVISIONS };

        for sub in 0..subdivs {
            let t = sub as f32 / SPLINE_SUBDIVISIONS as f32;
            let pos = catmull_rom(p0, p1, p2, p3, t);
            let tangent = catmull_rom_tangent(p0, p1, p2, p3, t).normalize_or_zero();

            // Interpolate SS type and color from the control points
            let (ss, color) = if t < 0.5 {
                (controls[seg].ss, controls[seg].color)
            } else {
                (controls[seg + 1].ss, controls[seg + 1].color)
            };

            spline_points.push(SplinePoint {
                pos,
                tangent,
                normal: Vec3::ZERO, // filled in next pass
                ss,
                color,
                ss_frac: 0.0, // filled in next pass
            });
        }
    }

    // Compute consistent normals using parallel transport
    compute_normals(&mut spline_points);

    // Compute ss_frac for arrow tapering on sheet segments
    compute_ss_fractions(&mut spline_points);

    spline_points
}

/// Catmull-Rom interpolation.
fn catmull_rom(p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3, t: f32) -> Vec3 {
    let t2 = t * t;
    let t3 = t2 * t;
    0.5 * ((2.0 * p1)
        + (-p0 + p2) * t
        + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t2
        + (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t3)
}

/// Catmull-Rom tangent (derivative).
fn catmull_rom_tangent(p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3, t: f32) -> Vec3 {
    let t2 = t * t;
    0.5 * ((-p0 + p2)
        + (4.0 * p0 - 10.0 * p1 + 8.0 * p2 - 2.0 * p3) * t
        + (-3.0 * p0 + 9.0 * p1 - 9.0 * p2 + 3.0 * p3) * t2)
}

/// Compute normals along the spline using parallel transport.
fn compute_normals(points: &mut [SplinePoint]) {
    if points.is_empty() {
        return;
    }

    // Initial normal: perpendicular to first tangent
    let t0 = points[0].tangent;
    let initial_up = if t0.y.abs() < 0.99 { Vec3::Y } else { Vec3::X };
    points[0].normal = t0.cross(initial_up).normalize_or_zero().cross(t0).normalize_or_zero();

    // Parallel transport
    for i in 1..points.len() {
        let prev_t = points[i - 1].tangent;
        let curr_t = points[i].tangent;
        let prev_n = points[i - 1].normal;

        // Rotate the previous normal to align with the new tangent
        let b = prev_t.cross(curr_t);
        if b.length_squared() < 1e-10 {
            // Tangents nearly parallel
            points[i].normal = prev_n;
        } else {
            let b = b.normalize();
            let angle = prev_t.dot(curr_t).clamp(-1.0, 1.0).acos();
            let rot = glam::Quat::from_axis_angle(b, angle);
            points[i].normal = (rot * prev_n).normalize();
        }
    }
}

/// Compute fractional position within consecutive SS segments (for sheet arrow tapering).
fn compute_ss_fractions(points: &mut [SplinePoint]) {
    let n = points.len();
    if n == 0 {
        return;
    }

    // Find runs of same SS type
    let mut run_start = 0;
    for i in 1..=n {
        let new_run = i == n || points[i].ss != points[run_start].ss;
        if new_run {
            let run_len = i - run_start;
            for j in run_start..i {
                points[j].ss_frac = (j - run_start) as f32 / (run_len as f32).max(1.0);
            }
            if i < n {
                run_start = i;
            }
        }
    }
}

/// Generate cartoon mesh from spline points.
fn generate_cartoon_mesh(
    spline: &[SplinePoint],
    vertices: &mut Vec<CartoonVertex>,
    indices: &mut Vec<u32>,
) {
    if spline.len() < 2 {
        return;
    }

    // For each spline point, generate a cross-section ring
    // The ring shape depends on the SS type:
    //   Loop  → circular tube
    //   Helix → wide flat ribbon (elliptical)
    //   Sheet → flat ribbon with arrow tapering at end

    let ring_size = TUBE_SEGMENTS;

    // Generate rings
    let mut rings: Vec<Vec<(Vec3, Vec3)>> = Vec::new(); // (position, normal) per ring vertex

    for pt in spline {
        let t = pt.tangent;
        let n = pt.normal;
        let b = t.cross(n).normalize_or_zero(); // binormal

        let ring = match pt.ss {
            SSType::Loop => {
                // Circular cross-section
                generate_circle_ring(pt.pos, n, b, TUBE_RADIUS, ring_size)
            }
            SSType::Helix => {
                // Elliptical cross-section (wide ribbon)
                generate_ellipse_ring(pt.pos, n, b, HELIX_WIDTH, HELIX_THICKNESS, ring_size)
            }
            SSType::Sheet => {
                // Flat ribbon with arrow taper
                let width = if pt.ss_frac > 0.7 {
                    // Arrow head: widen then taper to point
                    let arrow_t = (pt.ss_frac - 0.7) / 0.3;
                    let w = ARROW_HEAD_WIDTH * (1.0 - arrow_t);
                    w.max(0.05)
                } else {
                    SHEET_WIDTH
                };
                generate_ellipse_ring(pt.pos, n, b, width, SHEET_THICKNESS, ring_size)
            }
        };
        rings.push(ring);
    }

    // Emit vertices and stitch adjacent rings with triangles
    let base = vertices.len() as u32;

    for ring in &rings {
        for &(pos, normal) in ring {
            // Color will be set per-ring below
            vertices.push(CartoonVertex {
                position: pos.to_array(),
                normal: normal.to_array(),
                color: [0.0; 3], // placeholder
            });
        }
    }

    // Set colors per ring (each ring corresponds to one spline point)
    for (ring_idx, pt) in spline.iter().enumerate() {
        let start = base as usize + ring_idx * ring_size;
        for v in &mut vertices[start..start + ring_size] {
            v.color = pt.color;
        }
    }

    // Stitch rings
    let num_rings = rings.len();
    for r in 0..num_rings - 1 {
        let r0 = base + (r * ring_size) as u32;
        let r1 = base + ((r + 1) * ring_size) as u32;

        for i in 0..ring_size {
            let i0 = r0 + i as u32;
            let i1 = r0 + ((i + 1) % ring_size) as u32;
            let i2 = r1 + i as u32;
            let i3 = r1 + ((i + 1) % ring_size) as u32;

            indices.push(i0);
            indices.push(i2);
            indices.push(i1);

            indices.push(i1);
            indices.push(i2);
            indices.push(i3);
        }
    }

    // Cap at start
    let cap_center = vertices.len() as u32;
    let pt = &spline[0];
    vertices.push(CartoonVertex {
        position: pt.pos.to_array(),
        normal: (-pt.tangent).to_array(),
        color: pt.color,
    });
    for i in 0..ring_size {
        let i0 = base + i as u32;
        let i1 = base + ((i + 1) % ring_size) as u32;
        indices.push(cap_center);
        indices.push(i1);
        indices.push(i0);
    }

    // Cap at end
    let end_cap_center = vertices.len() as u32;
    let pt = &spline[spline.len() - 1];
    vertices.push(CartoonVertex {
        position: pt.pos.to_array(),
        normal: pt.tangent.to_array(),
        color: pt.color,
    });
    let last_ring_base = base + ((num_rings - 1) * ring_size) as u32;
    for i in 0..ring_size {
        let i0 = last_ring_base + i as u32;
        let i1 = last_ring_base + ((i + 1) % ring_size) as u32;
        indices.push(end_cap_center);
        indices.push(i0);
        indices.push(i1);
    }
}

/// Generate a circular cross-section ring.
fn generate_circle_ring(
    center: Vec3,
    normal: Vec3,
    binormal: Vec3,
    radius: f32,
    segments: usize,
) -> Vec<(Vec3, Vec3)> {
    let mut ring = Vec::with_capacity(segments);
    for i in 0..segments {
        let angle = (i as f32) * std::f32::consts::TAU / (segments as f32);
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        let dir = normal * cos_a + binormal * sin_a;
        let pos = center + dir * radius;
        ring.push((pos, dir));
    }
    ring
}

/// Generate an elliptical cross-section ring (for ribbons).
/// `half_width` is along the normal, `half_thickness` along the binormal.
fn generate_ellipse_ring(
    center: Vec3,
    normal: Vec3,
    binormal: Vec3,
    half_width: f32,
    half_thickness: f32,
    segments: usize,
) -> Vec<(Vec3, Vec3)> {
    let mut ring = Vec::with_capacity(segments);
    for i in 0..segments {
        let angle = (i as f32) * std::f32::consts::TAU / (segments as f32);
        let cos_a = angle.cos();
        let sin_a = angle.sin();

        // Position on ellipse
        let pos = center + normal * (cos_a * half_width) + binormal * (sin_a * half_thickness);

        // Normal to ellipse surface (not the same as the radial direction for an ellipse)
        // For an ellipse with semi-axes a,b: outward normal at angle θ is (cos θ / a, sin θ / b) normalized
        let nx = cos_a / half_width;
        let ny = sin_a / half_thickness;
        let n_len = (nx * nx + ny * ny).sqrt();
        let out_normal = (normal * (nx / n_len) + binormal * (ny / n_len)).normalize_or_zero();

        ring.push((pos, out_normal));
    }
    ring
}
