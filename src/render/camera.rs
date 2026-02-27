use glam::{Mat4, Quat, Vec3};

/// Arcball-style camera for 3D molecular viewing.
pub struct Camera {
    /// Camera target point (center of rotation).
    pub target: Vec3,
    /// Distance from target.
    pub distance: f32,
    /// Rotation quaternion.
    pub rotation: Quat,
    /// Field of view in radians.
    pub fov: f32,
    /// Near clip plane.
    pub near: f32,
    /// Far clip plane.
    pub far: f32,
}

impl Default for Camera {
    fn default() -> Self {
        Self {
            target: Vec3::ZERO,
            distance: 50.0,
            rotation: Quat::IDENTITY,
            fov: 45.0_f32.to_radians(),
            near: 0.1,
            far: 1000.0,
        }
    }
}

impl Camera {
    /// Compute the view matrix.
    pub fn view_matrix(&self) -> Mat4 {
        let eye = self.eye_position();
        let up = self.rotation * Vec3::Y;
        Mat4::look_at_rh(eye, self.target, up)
    }

    /// Compute the projection matrix for a given aspect ratio.
    pub fn projection_matrix(&self, aspect: f32) -> Mat4 {
        Mat4::perspective_rh(self.fov, aspect, self.near, self.far)
    }

    /// Get the eye (camera) position in world space.
    pub fn eye_position(&self) -> Vec3 {
        let forward = self.rotation * Vec3::NEG_Z;
        self.target - forward * self.distance
    }

    /// Rotate the camera by a delta in screen-space pixels.
    /// `dx` and `dy` are normalized mouse deltas.
    pub fn rotate(&mut self, dx: f32, dy: f32) {
        let sensitivity = 0.005;
        let angle_x = -dx * sensitivity;
        let angle_y = -dy * sensitivity;

        // Rotate around world Y axis for horizontal movement
        let rot_y = Quat::from_rotation_y(angle_x);
        // Rotate around camera-local X axis for vertical movement
        let right = self.rotation * Vec3::X;
        let rot_x = Quat::from_axis_angle(right, angle_y);

        self.rotation = (rot_y * rot_x * self.rotation).normalize();
    }

    /// Zoom by changing distance. Positive = zoom in, negative = zoom out.
    pub fn zoom(&mut self, delta: f32) {
        let factor = 1.0 - delta * 0.1;
        self.distance = (self.distance * factor).clamp(1.0, 500.0);
    }

    /// Pan the camera in the screen plane.
    pub fn pan(&mut self, dx: f32, dy: f32) {
        let scale = self.distance * 0.002;
        let right = self.rotation * Vec3::X;
        let up = self.rotation * Vec3::Y;
        self.target += right * (-dx * scale) + up * (dy * scale);
    }

    /// Reset to view a sphere centered at `center` with given `radius`.
    pub fn reset_to_fit(&mut self, center: [f32; 3], radius: f32) {
        self.target = Vec3::from_array(center);
        self.distance = radius * 2.5;
        self.rotation = Quat::IDENTITY;
        self.near = 0.1;
        self.far = (radius * 10.0).max(100.0);
    }
}
