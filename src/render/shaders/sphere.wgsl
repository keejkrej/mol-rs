// Sphere impostor shader: renders pixel-perfect spheres via ray-sphere intersection
// Each atom is a billboard quad (2 triangles from 6 vertices).
// The vertex shader expands the quad; the fragment shader does ray-sphere math.

struct Uniforms {
    view_proj: mat4x4<f32>,
    view: mat4x4<f32>,
    proj: mat4x4<f32>,
    eye_pos: vec4<f32>,
    light_dir: vec4<f32>,
    viewport_size: vec4<f32>,   // (width, height, 0, 0)
};

@group(0) @binding(0) var<uniform> u: Uniforms;

// Per-instance data: center + radius + color
struct InstanceInput {
    @location(0) center: vec3<f32>,
    @location(1) radius: f32,
    @location(2) color: vec3<f32>,
};

struct VertexOutput {
    @builtin(position) clip_pos: vec4<f32>,
    @location(0) color: vec3<f32>,
    @location(1) view_center: vec3<f32>,
    @location(2) view_pos: vec3<f32>,
    @location(3) radius: f32,
};

// We emit 6 vertices per instance (2 triangles forming a quad)
// vertex_index: 0-5

@vertex
fn vs_main(
    @builtin(vertex_index) vid: u32,
    instance: InstanceInput,
) -> VertexOutput {
    var out: VertexOutput;

    // Transform center to view space
    let view_center = (u.view * vec4<f32>(instance.center, 1.0)).xyz;

    // Billboard offsets: expand quad in view space
    // Quad vertices: (-1,-1), (1,-1), (1,1), (-1,-1), (1,1), (-1,1)
    var offsets = array<vec2<f32>, 6>(
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 1.0, -1.0),
        vec2<f32>( 1.0,  1.0),
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 1.0,  1.0),
        vec2<f32>(-1.0,  1.0),
    );

    let offset = offsets[vid];
    // Scale by radius with a small margin for edge antialiasing
    let scale = instance.radius * 1.2;
    let view_pos = view_center + vec3<f32>(offset.x * scale, offset.y * scale, 0.0);

    out.clip_pos = u.proj * vec4<f32>(view_pos, 1.0);
    out.color = instance.color;
    out.view_center = view_center;
    out.view_pos = view_pos;
    out.radius = instance.radius;
    return out;
}

@fragment
fn fs_main(in: VertexOutput) -> @location(0) vec4<f32> {
    // Ray-sphere intersection in view space
    // Ray origin = in.view_pos (on the billboard), direction = (0,0,-1) in view space
    let oc = in.view_pos - in.view_center;

    // For orthographic-like close approximation, ray dir is (0,0,-1)
    // But for perspective, ray goes from eye through fragment:
    let ray_dir = normalize(in.view_pos);

    let a = dot(ray_dir, ray_dir);
    let b = 2.0 * dot(oc, ray_dir);
    let c = dot(oc, oc) - in.radius * in.radius;
    let discriminant = b * b - 4.0 * a * c;

    if (discriminant < 0.0) {
        discard;
    }

    let t = (-b - sqrt(discriminant)) / (2.0 * a);
    let hit_pos = in.view_pos + ray_dir * t;
    let normal = normalize(hit_pos - in.view_center);

    // Blinn-Phong lighting in view space
    let light_dir_view = normalize((u.view * vec4<f32>(normalize(u.light_dir.xyz), 0.0)).xyz);
    let n_dot_l = max(dot(normal, -light_dir_view), 0.0);

    let view_dir = normalize(-hit_pos);
    let half_vec = normalize(view_dir - light_dir_view);
    let spec = pow(max(dot(normal, half_vec), 0.0), 32.0);

    let ambient = vec3<f32>(0.15, 0.15, 0.15);
    let diffuse = in.color * n_dot_l * 0.75;
    let specular = vec3<f32>(0.4, 0.4, 0.4) * spec;

    let final_color = ambient + diffuse + specular;

    // Write correct depth
    let clip = u.proj * vec4<f32>(hit_pos, 1.0);
    let ndc_depth = clip.z / clip.w;

    return vec4<f32>(clamp(final_color, vec3<f32>(0.0), vec3<f32>(1.0)), 1.0);
}
