// Cylinder (stick) shader with Blinn-Phong lighting
// Uses pre-built mesh geometry: vertices have position, normal, and color.

struct Uniforms {
    view_proj: mat4x4<f32>,
    view: mat4x4<f32>,
    proj: mat4x4<f32>,
    eye_pos: vec4<f32>,
    light_dir: vec4<f32>,
    viewport_size: vec4<f32>,
};

@group(0) @binding(0) var<uniform> u: Uniforms;

struct VertexInput {
    @location(0) position: vec3<f32>,
    @location(1) normal: vec3<f32>,
    // Instance attributes
    @location(2) start: vec3<f32>,
    @location(3) end: vec3<f32>,
    @location(4) color: vec3<f32>,
    @location(5) radius: f32,
};

struct VertexOutput {
    @builtin(position) clip_pos: vec4<f32>,
    @location(0) color: vec3<f32>,
    @location(1) world_normal: vec3<f32>,
    @location(2) world_pos: vec3<f32>,
};

@vertex
fn vs_main(in: VertexInput) -> VertexOutput {
    var out: VertexOutput;

    let axis_vec = in.end - in.start;
    let len = length(axis_vec);
    var dir = vec3<f32>(0.0, 1.0, 0.0);
    if (len > 0.0001) {
        dir = axis_vec / len;
    }

    // Build orthonormal basis
    var temp_up = vec3<f32>(0.0, 1.0, 0.0);
    if (abs(dir.y) > 0.99) {
        temp_up = vec3<f32>(1.0, 0.0, 0.0);
    }
    let right = normalize(cross(dir, temp_up));
    let fwd = normalize(cross(right, dir));

    // Transform position
    // Local Y is along the stick axis (0 to 1)
    // Local X and Z are radial
    let p = in.position;
    let pos_world = in.start 
        + right * (p.x * in.radius) 
        + dir * (p.y * len) 
        + fwd * (p.z * in.radius);

    // Transform normal
    // Just rotate
    let norm_world = normalize(
          right * in.normal.x 
        + dir * in.normal.y 
        + fwd * in.normal.z
    );

    out.clip_pos = u.view_proj * vec4<f32>(pos_world, 1.0);
    out.color = in.color;
    out.world_normal = norm_world;
    out.world_pos = pos_world;
    return out;
}

@fragment
fn fs_main(in: VertexOutput) -> @location(0) vec4<f32> {
    let normal = normalize(in.world_normal);
    let light_dir = normalize(u.light_dir.xyz);
    let view_dir = normalize(u.eye_pos.xyz - in.world_pos);

    // Diffuse
    let n_dot_l = max(dot(normal, -light_dir), 0.0);

    // Specular (Blinn-Phong)
    let half_vec = normalize(view_dir - light_dir);
    let spec = pow(max(dot(normal, half_vec), 0.0), 32.0);

    let ambient = vec3<f32>(0.15, 0.15, 0.15);
    let diffuse = in.color * n_dot_l * 0.75;
    let specular = vec3<f32>(0.3, 0.3, 0.3) * spec;

    let final_color = ambient + diffuse + specular;
    return vec4<f32>(clamp(final_color, vec3<f32>(0.0), vec3<f32>(1.0)), 1.0);
}
