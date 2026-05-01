use std::collections::BTreeSet;

use crate::core::bond::BondInfo;
use crate::core::molecule::Molecule;
use crate::render::camera::Camera;
use crate::scene::color::{apply_color_scheme, ColorScheme};
use crate::selection::{evaluate_with_coords, parser::Selector};

#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct Measurement {
    pub p1: [f32; 3],
    pub p2: [f32; 3],
    pub distance: f32,
    pub label: String,
}

pub struct Scene {
    pub molecules: Vec<Molecule>,
    pub measurements: Vec<Measurement>,
    pub camera: Camera,
    pub color_scheme: ColorScheme,
    pub current_state: usize,
    pub all_states: bool,
    /// True when geometry buffers need rebuilding.
    pub geometry_dirty: bool,
    /// Background color [r, g, b].
    pub bg_color: [f32; 3],
}

impl Default for Scene {
    fn default() -> Self {
        Self {
            molecules: Vec::new(),
            measurements: Vec::new(),
            camera: Camera::default(),
            color_scheme: ColorScheme::ByElement,
            current_state: 1,
            all_states: false,
            geometry_dirty: false,
            bg_color: [0.0, 0.0, 0.0],
        }
    }
}

impl Scene {
    pub fn max_state_count(&self) -> usize {
        self.molecules
            .iter()
            .map(|m| m.state_count())
            .max()
            .unwrap_or(1)
            .max(1)
    }

    pub fn count_selected_atoms(&self, sel: &Selector, state: usize) -> usize {
        self.molecules
            .iter()
            .map(|mol| {
                let coords = mol.coords_for_state(state);
                evaluate_with_coords(sel, mol, coords)
                    .into_iter()
                    .enumerate()
                    .filter(|(idx, selected)| *selected && *idx < coords.len())
                    .count()
            })
            .sum()
    }

    pub fn count_selection_states(&self, sel: &Selector) -> usize {
        (1..=self.max_state_count())
            .filter(|state| {
                self.molecules.iter().any(|mol| {
                    let coords = mol.coords_for_state(*state);
                    !coords.is_empty()
                        && evaluate_with_coords(sel, mol, coords)
                            .into_iter()
                            .enumerate()
                            .any(|(idx, selected)| selected && idx < coords.len())
                })
            })
            .count()
    }

    pub fn set_state_clamped(&mut self, state: usize) {
        let max_state = self.max_state_count();
        let clamped = state.clamp(1, max_state);
        if self.current_state != clamped {
            self.current_state = clamped;
            self.geometry_dirty = true;
        }
    }

    pub fn next_state(&mut self) {
        self.set_state_clamped(self.current_state + 1);
    }

    pub fn prev_state(&mut self) {
        self.set_state_clamped(self.current_state.saturating_sub(1));
    }

    /// Add a molecule to the scene and adjust camera to fit.
    pub fn add_molecule(&mut self, mut mol: Molecule) {
        let requested_state = self.current_state;
        apply_color_scheme(&mut mol, self.color_scheme);
        let center = mol.centroid_for_state(requested_state);
        let radius = mol.radius_for_state(requested_state);
        self.molecules.push(mol);
        self.set_state_clamped(requested_state);

        // Fit camera to the newly loaded molecule
        self.camera.reset_to_fit(center, radius);
        self.geometry_dirty = true;
    }

    /// Copy an existing molecular object by exact name.
    pub fn copy_object(&mut self, target: &str, source: &str) -> bool {
        let Some(mut copy) = self
            .molecules
            .iter()
            .find(|mol| mol.name == source)
            .cloned()
        else {
            return false;
        };

        copy.name = target.to_string();
        self.molecules.push(copy);
        self.set_state_clamped(self.current_state);
        self.geometry_dirty = true;
        true
    }

    /// Rename a molecular object by exact name.
    pub fn rename_object(
        &mut self,
        old_name: &str,
        new_name: &str,
    ) -> Result<(), RenameObjectError> {
        if old_name == new_name {
            return Ok(());
        }
        if self.molecules.iter().any(|mol| mol.name == new_name) {
            return Err(RenameObjectError::NameExists);
        }

        let Some(mol) = self.molecules.iter_mut().find(|mol| mol.name == old_name) else {
            return Err(RenameObjectError::NotFound);
        };
        mol.name = new_name.to_string();
        Ok(())
    }

    /// Create a molecule object from selected atoms.
    ///
    /// `selection_state` chooses which coordinates are used for state-dependent
    /// selectors like `present` or distance selectors. When `source_state` is
    /// `None`, all scene states are copied. Otherwise, only that source state is
    /// copied into state 1 of the new object. With `extract`, selected atoms are
    /// also removed from their source molecules.
    pub fn create_object_from_selection(
        &mut self,
        name: &str,
        sel: &Selector,
        selection_state: usize,
        source_state: Option<usize>,
        extract: bool,
    ) -> usize {
        let original_len = self.molecules.len();
        if original_len == 0 {
            return 0;
        }

        let masks: Vec<Vec<bool>> = self
            .molecules
            .iter()
            .take(original_len)
            .map(|mol| evaluate_with_coords(sel, mol, mol.coords_for_state(selection_state)))
            .collect();

        let state_count = source_state
            .map(|_| 1)
            .unwrap_or_else(|| self.max_state_count());
        let mut out = Molecule::new(name.to_string());
        out.coord_sets = vec![Vec::new(); state_count];

        for (mol, mask) in self.molecules.iter().take(original_len).zip(masks.iter()) {
            let selection_coords = mol.coords_for_state(selection_state);
            let mut old_to_new = vec![None; mol.atoms.len()];

            for (idx, atom) in mol.atoms.iter().cloned().enumerate() {
                if mask.get(idx).copied().unwrap_or(false) && idx < selection_coords.len() {
                    old_to_new[idx] = Some(out.atoms.len());
                    out.atoms.push(atom);
                }
            }

            for output_state in 0..state_count {
                let state = source_state.unwrap_or(output_state + 1);
                let coords = mol.coords_for_state(state);
                for (idx, mapped) in old_to_new.iter().enumerate() {
                    if mapped.is_some() {
                        if let Some(coord) = coords.get(idx).copied() {
                            out.coord_sets[output_state].push(coord);
                        }
                    }
                }
            }

            for bond in &mol.bonds {
                let Some(atom_a) = old_to_new.get(bond.atom_a).copied().flatten() else {
                    continue;
                };
                let Some(atom_b) = old_to_new.get(bond.atom_b).copied().flatten() else {
                    continue;
                };
                out.bonds.push(BondInfo {
                    atom_a,
                    atom_b,
                    order: bond.order,
                });
            }
        }

        let created = out.atoms.len();
        if created == 0 {
            return 0;
        }

        out.build_residues();

        if extract {
            for (mol, mask) in self
                .molecules
                .iter_mut()
                .take(original_len)
                .zip(masks.iter())
            {
                mol.remove_atoms(mask);
            }
        }

        self.molecules.push(out);
        self.set_state_clamped(self.current_state);
        self.geometry_dirty = true;
        created
    }

    /// Recolor all molecules with the given scheme.
    pub fn set_color_scheme(&mut self, scheme: ColorScheme) {
        self.color_scheme = scheme;
        for mol in &mut self.molecules {
            apply_color_scheme(mol, scheme);
        }
        self.geometry_dirty = true;
    }

    /// Set object visibility for names matching a PyMOL-like object pattern.
    pub fn set_object_visibility(&mut self, pattern: &str, visible: bool) -> usize {
        let pattern = pattern.trim();
        let pattern = if pattern.is_empty() { "all" } else { pattern };
        let mut changed = 0usize;

        for mol in &mut self.molecules {
            if object_name_matches(pattern, &mol.name) && mol.visible != visible {
                mol.visible = visible;
                changed += 1;
            }
        }

        if changed > 0 {
            self.geometry_dirty = true;
        }
        changed
    }

    /// Delete molecular objects matching a PyMOL-like object pattern.
    pub fn delete_objects(&mut self, pattern: &str) -> usize {
        let pattern = pattern.trim();
        if pattern.is_empty() {
            return 0;
        }

        let old_len = self.molecules.len();
        self.molecules
            .retain(|mol| !object_name_matches(pattern, &mol.name));
        let deleted = old_len - self.molecules.len();

        if deleted > 0 {
            if is_all_object_pattern(pattern) || self.molecules.is_empty() {
                self.measurements.clear();
            }
            self.set_state_clamped(self.current_state);
            self.geometry_dirty = true;
        }

        deleted
    }

    /// Sort atoms in molecular objects matching a PyMOL-like object pattern.
    pub fn sort_objects(&mut self, pattern: &str) -> usize {
        let pattern = pattern.trim();
        let pattern = if pattern.is_empty() { "all" } else { pattern };
        let mut matched = 0usize;
        let mut changed = false;

        for mol in &mut self.molecules {
            if object_name_matches(pattern, &mol.name) {
                matched += 1;
                changed |= mol.sort_atoms();
            }
        }

        if changed {
            self.geometry_dirty = true;
        }

        matched
    }

    /// Reorder molecular objects by PyMOL-like name patterns.
    pub fn order_objects(
        &mut self,
        names: &str,
        sort_selected: bool,
        location: OrderLocation,
    ) -> usize {
        let patterns: Vec<&str> = names.split_whitespace().collect();
        if patterns.is_empty() {
            return 0;
        }

        let matched_indices: Vec<usize> = self
            .molecules
            .iter()
            .enumerate()
            .filter_map(|(idx, mol)| {
                patterns
                    .iter()
                    .any(|pattern| object_name_matches(pattern, &mol.name))
                    .then_some(idx)
            })
            .collect();

        let matched_count = matched_indices.len();
        if matched_count == 0 {
            return 0;
        }

        let original_names: Vec<String> =
            self.molecules.iter().map(|mol| mol.name.clone()).collect();
        let first_match = matched_indices[0];
        let mut selected = Vec::with_capacity(matched_count);
        let mut remaining = Vec::with_capacity(self.molecules.len() - matched_count);

        for (idx, mol) in self.molecules.drain(..).enumerate() {
            if matched_indices.contains(&idx) {
                selected.push(mol);
            } else {
                remaining.push(mol);
            }
        }

        if sort_selected {
            selected.sort_by_key(|mol| mol.name.to_ascii_lowercase());
        }

        let insert_at = match location {
            OrderLocation::Top => 0,
            OrderLocation::Current => first_match.min(remaining.len()),
            OrderLocation::Bottom => remaining.len(),
        };

        remaining.splice(insert_at..insert_at, selected);
        let changed = original_names
            != remaining
                .iter()
                .map(|mol| mol.name.clone())
                .collect::<Vec<_>>();
        self.molecules = remaining;

        if changed {
            self.geometry_dirty = true;
        }

        matched_count
    }

    pub fn clip_camera(
        &mut self,
        mode: ClipMode,
        distance: f32,
        selection: Option<&Selector>,
        state: usize,
    ) -> Result<(), ClipError> {
        let mut near = self.camera.near;
        let mut far = self.camera.far;

        match mode {
            ClipMode::Near => near -= distance,
            ClipMode::Far => far -= distance,
            ClipMode::Move => {
                near -= distance;
                far -= distance;
            }
            ClipMode::Slab => {
                let center = if let Some(sel) = selection {
                    self.selection_clip_depth(sel, state)
                        .map(|(min_depth, max_depth)| (min_depth + max_depth) * 0.5)
                } else {
                    None
                }
                .unwrap_or((near + far) * 0.5);
                let half_width = (distance.max(0.0)) * 0.5;
                near = center - half_width;
                far = center + half_width;
            }
            ClipMode::Atoms => {
                let sel = selection.unwrap_or(&Selector::All);
                let Some((min_depth, max_depth)) = self.selection_clip_depth(sel, state) else {
                    return Err(ClipError::EmptySelection);
                };
                let buffer = distance.max(0.0);
                near = min_depth - buffer;
                far = max_depth + buffer;
            }
            ClipMode::NearSet => near = distance,
            ClipMode::FarSet => far = distance,
        }

        set_camera_clip(&mut self.camera.near, &mut self.camera.far, near, far);
        self.geometry_dirty = true;
        Ok(())
    }

    fn selection_clip_depth(&self, sel: &Selector, state: usize) -> Option<(f32, f32)> {
        let eye = self.camera.eye_position();
        let forward = self.camera.rotation * glam::Vec3::NEG_Z;
        let mut min_depth = f32::INFINITY;
        let mut max_depth = f32::NEG_INFINITY;
        let mut found = false;

        for mol in &self.molecules {
            let coords = mol.coords_for_state(state);
            let mask = evaluate_with_coords(sel, mol, coords);
            for (idx, selected) in mask.iter().enumerate() {
                if *selected {
                    let Some(coord) = coords.get(idx).copied() else {
                        continue;
                    };
                    let depth = (glam::Vec3::from_array(coord) - eye).dot(forward);
                    min_depth = min_depth.min(depth);
                    max_depth = max_depth.max(depth);
                    found = true;
                }
            }
        }

        found.then_some((min_depth, max_depth))
    }

    pub fn selection_bounds(&self, sel: &Selector, state: usize) -> Option<SelectionBounds> {
        let mut min = [f32::INFINITY; 3];
        let mut max = [f32::NEG_INFINITY; 3];
        let mut atom_count = 0usize;

        for mol in &self.molecules {
            let coords = mol.coords_for_state(state);
            let mask = evaluate_with_coords(sel, mol, coords);

            for (idx, selected) in mask.iter().enumerate() {
                if *selected {
                    let Some(p) = coords.get(idx).copied() else {
                        continue;
                    };
                    min[0] = min[0].min(p[0]);
                    min[1] = min[1].min(p[1]);
                    min[2] = min[2].min(p[2]);
                    max[0] = max[0].max(p[0]);
                    max[1] = max[1].max(p[1]);
                    max[2] = max[2].max(p[2]);
                    atom_count += 1;
                }
            }
        }

        (atom_count > 0).then_some(SelectionBounds {
            min,
            max,
            atom_count,
        })
    }

    pub fn selection_coords(&self, sel: &Selector, state: usize) -> Vec<[f32; 3]> {
        let mut points = Vec::new();

        for mol in &self.molecules {
            let coords = mol.coords_for_state(state);
            let mask = evaluate_with_coords(sel, mol, coords);

            for (idx, selected) in mask.iter().enumerate() {
                if *selected {
                    if let Some(p) = coords.get(idx).copied() {
                        points.push(p);
                    }
                }
            }
        }

        points
    }

    pub fn selection_chains(&self, sel: &Selector, state: usize) -> Vec<char> {
        let mut chains = BTreeSet::new();

        for mol in &self.molecules {
            let coords = mol.coords_for_state(state);
            let mask = evaluate_with_coords(sel, mol, coords);

            for (idx, selected) in mask.iter().enumerate() {
                if *selected && idx < coords.len() {
                    if let Some(atom) = mol.atoms.get(idx) {
                        chains.insert(atom.chain);
                    }
                }
            }
        }

        chains.into_iter().collect()
    }

    pub fn selection_indices(&self, sel: &Selector, state: usize) -> Vec<AtomIndex> {
        let mut indices = Vec::new();

        for mol in &self.molecules {
            let coords = mol.coords_for_state(state);
            let mask = evaluate_with_coords(sel, mol, coords);

            for (idx, selected) in mask.iter().enumerate() {
                if *selected && idx < coords.len() {
                    indices.push(AtomIndex {
                        object: mol.name.clone(),
                        index: idx + 1,
                    });
                }
            }
        }

        indices
    }

    pub fn selection_ids(&self, sel: &Selector, state: usize) -> Vec<AtomId> {
        let mut ids = Vec::new();

        for mol in &self.molecules {
            let coords = mol.coords_for_state(state);
            let mask = evaluate_with_coords(sel, mol, coords);

            for (idx, selected) in mask.iter().enumerate() {
                if *selected && idx < coords.len() {
                    if let Some(atom) = mol.atoms.get(idx) {
                        ids.push(AtomId {
                            object: mol.name.clone(),
                            id: if atom.serial == 0 {
                                (idx + 1) as u32
                            } else {
                                atom.serial
                            },
                        });
                    }
                }
            }
        }

        ids
    }

    pub fn object_names(
        &self,
        enabled_only: bool,
        selection: Option<&Selector>,
        state: usize,
    ) -> Vec<String> {
        self.molecules
            .iter()
            .filter(|mol| !enabled_only || mol.visible)
            .filter(|mol| {
                let Some(sel) = selection else {
                    return true;
                };
                let coords = mol.coords_for_state(state);
                evaluate_with_coords(sel, mol, coords)
                    .into_iter()
                    .enumerate()
                    .any(|(idx, selected)| selected && idx < coords.len())
            })
            .map(|mol| mol.name.clone())
            .collect()
    }

    pub fn object_type(&self, name: &str) -> Option<&'static str> {
        self.molecules
            .iter()
            .any(|mol| mol.name == name)
            .then_some("object:molecule")
    }

    pub fn object_names_of_type(&self, object_type: &str) -> Vec<String> {
        if object_type == "object:molecule" {
            self.object_names(false, None, self.current_state)
        } else {
            Vec::new()
        }
    }

    pub fn closest_distance(
        &self,
        sel1: &Selector,
        sel2: &Selector,
        state: usize,
    ) -> Option<SelectionDistance> {
        type SelPoint = (usize, usize, [f32; 3]);
        let mut pts1: Vec<SelPoint> = Vec::new();
        let mut pts2: Vec<SelPoint> = Vec::new();

        for (mol_idx, mol) in self.molecules.iter().enumerate() {
            let coords = mol.coords_for_state(state);
            let mask1 = evaluate_with_coords(sel1, mol, coords);
            let mask2 = evaluate_with_coords(sel2, mol, coords);
            for (idx, coord) in coords.iter().enumerate() {
                if mask1.get(idx).copied().unwrap_or(false) {
                    pts1.push((mol_idx, idx, *coord));
                }
                if mask2.get(idx).copied().unwrap_or(false) {
                    pts2.push((mol_idx, idx, *coord));
                }
            }
        }

        let mut best: Option<SelectionDistance> = None;
        for (mol1, idx1, p1) in &pts1 {
            for (mol2, idx2, p2) in &pts2 {
                if mol1 == mol2 && idx1 == idx2 {
                    continue;
                }
                let dx = p1[0] - p2[0];
                let dy = p1[1] - p2[1];
                let dz = p1[2] - p2[2];
                let distance2 = dx * dx + dy * dy + dz * dz;

                match best {
                    Some(existing) if distance2 >= existing.distance * existing.distance => {}
                    _ => {
                        best = Some(SelectionDistance {
                            p1: *p1,
                            p2: *p2,
                            distance: distance2.sqrt(),
                        });
                    }
                }
            }
        }

        best
    }

    pub fn selection_angle(
        &self,
        sel1: &Selector,
        sel2: &Selector,
        sel3: &Selector,
        state: usize,
    ) -> Result<f32, SelectionPointError> {
        let p1 = self.single_selection_point(sel1, state)?;
        let p2 = self.single_selection_point(sel2, state)?;
        let p3 = self.single_selection_point(sel3, state)?;

        let v1 = [p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]];
        let v2 = [p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2]];
        let len1 = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt();
        let len2 = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]).sqrt();
        if len1 == 0.0 || len2 == 0.0 {
            return Err(SelectionPointError::Degenerate);
        }

        let dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
        let cosine = (dot / (len1 * len2)).clamp(-1.0, 1.0);
        Ok(cosine.acos().to_degrees())
    }

    pub fn selection_dihedral(
        &self,
        sel1: &Selector,
        sel2: &Selector,
        sel3: &Selector,
        sel4: &Selector,
        state: usize,
    ) -> Result<f32, SelectionPointError> {
        let p1 = self.single_selection_point(sel1, state)?;
        let p2 = self.single_selection_point(sel2, state)?;
        let p3 = self.single_selection_point(sel3, state)?;
        let p4 = self.single_selection_point(sel4, state)?;

        let b0 = [p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]];
        let b1 = [p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2]];
        let b2 = [p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2]];

        let len_b1 = vec_len(b1);
        if len_b1 == 0.0 {
            return Err(SelectionPointError::Degenerate);
        }
        let b1n = [b1[0] / len_b1, b1[1] / len_b1, b1[2] / len_b1];

        let v = vec_sub(b0, vec_scale(b1n, dot3(b0, b1n)));
        let w = vec_sub(b2, vec_scale(b1n, dot3(b2, b1n)));
        if vec_len(v) == 0.0 || vec_len(w) == 0.0 {
            return Err(SelectionPointError::Degenerate);
        }

        let x = dot3(v, w);
        let y = dot3(cross3(b1n, v), w);
        Ok(y.atan2(x).to_degrees())
    }

    fn single_selection_point(
        &self,
        sel: &Selector,
        state: usize,
    ) -> Result<[f32; 3], SelectionPointError> {
        match self.selection_coords(sel, state).as_slice() {
            [] => Err(SelectionPointError::Empty),
            [point] => Ok(*point),
            _ => Err(SelectionPointError::Multiple),
        }
    }

    /// Compute bounding-sphere data for a selection in the current state.
    pub fn selection_extent(&self, sel: &Selector) -> Option<SelectionExtent> {
        let bounds = self.selection_bounds(sel, self.current_state)?;
        let points = self.selection_coords(sel, self.current_state);

        let center = [
            (bounds.min[0] + bounds.max[0]) * 0.5,
            (bounds.min[1] + bounds.max[1]) * 0.5,
            (bounds.min[2] + bounds.max[2]) * 0.5,
        ];
        let radius = points
            .iter()
            .map(|p| {
                let dx = p[0] - center[0];
                let dy = p[1] - center[1];
                let dz = p[2] - center[2];
                dx * dx + dy * dy + dz * dz
            })
            .fold(0.0f32, f32::max)
            .sqrt()
            .max(1.0);

        Some(SelectionExtent {
            center,
            radius,
            atom_count: bounds.atom_count,
        })
    }
}

fn dot3(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross3(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn vec_len(v: [f32; 3]) -> f32 {
    dot3(v, v).sqrt()
}

fn vec_scale(v: [f32; 3], scale: f32) -> [f32; 3] {
    [v[0] * scale, v[1] * scale, v[2] * scale]
}

fn vec_sub(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn object_name_matches(pattern: &str, name: &str) -> bool {
    let pattern = pattern.trim();
    if is_all_object_pattern(pattern) {
        return true;
    }

    wildcard_match_ci(pattern, name)
}

fn is_all_object_pattern(pattern: &str) -> bool {
    pattern.eq_ignore_ascii_case("all") || pattern == "*"
}

fn set_camera_clip(camera_near: &mut f32, camera_far: &mut f32, near: f32, far: f32) {
    let near = near.max(0.01);
    let far = far.max(near + 0.01);
    *camera_near = near;
    *camera_far = far;
}

fn wildcard_match_ci(pattern: &str, text: &str) -> bool {
    let pattern = pattern.to_ascii_lowercase();
    let text = text.to_ascii_lowercase();
    let pattern = pattern.as_bytes();
    let text = text.as_bytes();

    let (mut p, mut t) = (0usize, 0usize);
    let mut star = None;
    let mut star_text = 0usize;

    while t < text.len() {
        if p < pattern.len() && pattern[p] == text[t] {
            p += 1;
            t += 1;
        } else if p < pattern.len() && pattern[p] == b'*' {
            star = Some(p);
            p += 1;
            star_text = t;
        } else if let Some(star_idx) = star {
            p = star_idx + 1;
            star_text += 1;
            t = star_text;
        } else {
            return false;
        }
    }

    while p < pattern.len() && pattern[p] == b'*' {
        p += 1;
    }

    p == pattern.len()
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SelectionBounds {
    pub min: [f32; 3],
    pub max: [f32; 3],
    pub atom_count: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RenameObjectError {
    NotFound,
    NameExists,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OrderLocation {
    Top,
    Current,
    Bottom,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClipMode {
    Near,
    Far,
    Move,
    Slab,
    Atoms,
    NearSet,
    FarSet,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClipError {
    EmptySelection,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomIndex {
    pub object: String,
    pub index: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomId {
    pub object: String,
    pub id: u32,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SelectionDistance {
    pub p1: [f32; 3],
    pub p2: [f32; 3],
    pub distance: f32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SelectionPointError {
    Empty,
    Multiple,
    Degenerate,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SelectionExtent {
    pub center: [f32; 3],
    pub radius: f32,
    pub atom_count: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::atom::AtomInfo;
    use crate::selection::parser::Selector;

    fn atom(chain: char, serial: u32) -> AtomInfo {
        AtomInfo {
            chain,
            serial,
            ..AtomInfo::default()
        }
    }

    #[test]
    fn selection_extent_uses_current_state_and_selection() {
        let mut mol = Molecule::new("extent".to_string());
        mol.atoms = vec![atom('A', 1), atom('A', 2), atom('B', 3)];
        mol.coord_sets = vec![
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [100.0, 0.0, 0.0]],
            vec![[10.0, 0.0, 0.0], [14.0, 0.0, 0.0], [100.0, 0.0, 0.0]],
        ];

        let mut scene = Scene::default();
        scene.add_molecule(mol);
        scene.set_state_clamped(2);

        let extent = scene.selection_extent(&Selector::Chain('A')).unwrap();

        assert_eq!(extent.center, [12.0, 0.0, 0.0]);
        assert_eq!(extent.radius, 2.0);
        assert_eq!(extent.atom_count, 2);
    }

    #[test]
    fn selection_bounds_uses_requested_state_and_selection() {
        let mut mol = Molecule::new("extent".to_string());
        mol.atoms = vec![atom('A', 1), atom('A', 2), atom('B', 3)];
        mol.coord_sets = vec![
            vec![[-1.0, 2.0, 0.0], [2.0, -3.0, 4.0], [100.0, 0.0, 0.0]],
            vec![[10.0, 0.0, 0.0], [14.0, 1.0, -2.0], [100.0, 0.0, 0.0]],
        ];

        let mut scene = Scene::default();
        scene.add_molecule(mol);

        let bounds = scene.selection_bounds(&Selector::Chain('A'), 1).unwrap();

        assert_eq!(bounds.min, [-1.0, -3.0, 0.0]);
        assert_eq!(bounds.max, [2.0, 2.0, 4.0]);
        assert_eq!(bounds.atom_count, 2);
    }

    #[test]
    fn selection_coords_uses_requested_state_and_present_atoms() {
        let mut mol = Molecule::new("coords".to_string());
        mol.atoms = vec![atom('A', 1), atom('A', 2), atom('B', 3)];
        mol.coord_sets = vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            Vec::new(),
        ];

        let mut scene = Scene::default();
        scene.add_molecule(mol);

        assert_eq!(
            scene.selection_coords(&Selector::Chain('A'), 1),
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
        );
        assert!(scene.selection_coords(&Selector::Chain('A'), 2).is_empty());
    }

    #[test]
    fn selection_chains_returns_sorted_present_chains() {
        let mut mol = Molecule::new("chains".to_string());
        mol.atoms = vec![atom('B', 1), atom('A', 2), atom('B', 3)];
        mol.coord_sets = vec![vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]];

        let mut scene = Scene::default();
        scene.add_molecule(mol);

        assert_eq!(scene.selection_chains(&Selector::All, 1), vec!['A', 'B']);
        assert_eq!(scene.selection_chains(&Selector::Chain('B'), 1), vec!['B']);
    }

    #[test]
    fn selection_indices_returns_present_atom_indices() {
        let mut mol = Molecule::new("idx".to_string());
        mol.atoms = vec![atom('A', 1), atom('B', 2), atom('B', 3)];
        mol.coord_sets = vec![vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]];

        let mut scene = Scene::default();
        scene.add_molecule(mol);

        assert_eq!(
            scene.selection_indices(&Selector::All, 1),
            vec![
                AtomIndex {
                    object: "idx".to_string(),
                    index: 1,
                },
                AtomIndex {
                    object: "idx".to_string(),
                    index: 2,
                },
            ]
        );
        assert_eq!(
            scene.selection_indices(&Selector::Chain('B'), 1),
            vec![AtomIndex {
                object: "idx".to_string(),
                index: 2,
            }]
        );
    }

    #[test]
    fn selection_ids_returns_present_atom_serials() {
        let mut mol = Molecule::new("ids".to_string());
        mol.atoms = vec![atom('A', 10), atom('B', 20), atom('B', 30)];
        mol.coord_sets = vec![vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]];

        let mut scene = Scene::default();
        scene.add_molecule(mol);

        assert_eq!(
            scene.selection_ids(&Selector::All, 1),
            vec![
                AtomId {
                    object: "ids".to_string(),
                    id: 10,
                },
                AtomId {
                    object: "ids".to_string(),
                    id: 20,
                },
            ]
        );
        assert_eq!(
            scene.selection_ids(&Selector::Chain('B'), 1),
            vec![AtomId {
                object: "ids".to_string(),
                id: 20,
            }]
        );

        scene.molecules[0].atoms[1].serial = 0;
        assert_eq!(
            scene.selection_ids(&Selector::Chain('B'), 1),
            vec![AtomId {
                object: "ids".to_string(),
                id: 2,
            }]
        );
    }

    #[test]
    fn copy_object_duplicates_exact_source() {
        let mut scene = Scene::default();
        let mut mol = Molecule::new("source".to_string());
        mol.atoms = vec![atom('A', 1)];
        mol.coord_sets = vec![vec![[1.0, 2.0, 3.0]]];
        scene.add_molecule(mol);
        scene.geometry_dirty = false;

        assert!(scene.copy_object("copy", "source"));
        assert_eq!(scene.molecules.len(), 2);
        assert_eq!(scene.molecules[1].name, "copy");
        assert_eq!(scene.molecules[1].atoms.len(), 1);
        assert_eq!(scene.molecules[1].coord_sets[0], vec![[1.0, 2.0, 3.0]]);
        assert!(scene.geometry_dirty);

        assert!(!scene.copy_object("missing_copy", "missing"));
        assert_eq!(scene.molecules.len(), 2);
    }

    #[test]
    fn rename_object_updates_exact_source_and_rejects_conflicts() {
        let mut scene = Scene::default();
        scene.molecules.push(Molecule::new("source".to_string()));
        scene.molecules.push(Molecule::new("existing".to_string()));

        assert_eq!(
            scene.rename_object("missing", "renamed"),
            Err(RenameObjectError::NotFound)
        );
        assert_eq!(
            scene.rename_object("source", "existing"),
            Err(RenameObjectError::NameExists)
        );

        assert_eq!(scene.rename_object("source", "renamed"), Ok(()));
        assert_eq!(scene.molecules[0].name, "renamed");
        assert_eq!(scene.molecules[1].name, "existing");
        assert_eq!(scene.object_type("source"), None);
        assert_eq!(scene.object_type("renamed"), Some("object:molecule"));
    }

    #[test]
    fn create_object_from_selection_copies_atoms_states_and_bonds() {
        let mut mol = Molecule::new("source".to_string());
        mol.atoms = vec![atom('A', 10), atom('A', 20), atom('B', 30)];
        mol.coord_sets = vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec![[0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [2.0, 1.0, 0.0]],
        ];
        mol.bonds = vec![
            BondInfo {
                atom_a: 0,
                atom_b: 1,
                order: 1,
            },
            BondInfo {
                atom_a: 1,
                atom_b: 2,
                order: 1,
            },
        ];
        mol.build_residues();

        let mut scene = Scene::default();
        scene.add_molecule(mol);
        scene.geometry_dirty = false;

        let created =
            scene.create_object_from_selection("subset", &Selector::Chain('A'), 1, None, false);

        assert_eq!(created, 2);
        assert_eq!(scene.molecules.len(), 2);
        let subset = &scene.molecules[1];
        assert_eq!(subset.name, "subset");
        assert_eq!(subset.atoms.len(), 2);
        assert_eq!(subset.coord_sets.len(), 2);
        assert_eq!(
            subset.coord_sets,
            vec![
                vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
                vec![[0.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
            ]
        );
        assert_eq!(subset.bonds.len(), 1);
        assert_eq!(subset.bonds[0].atom_a, 0);
        assert_eq!(subset.bonds[0].atom_b, 1);
        assert_eq!(subset.residues.len(), 1);
        assert_eq!(scene.molecules[0].atoms.len(), 3);
        assert!(scene.geometry_dirty);
    }

    #[test]
    fn extract_object_from_selection_removes_source_atoms() {
        let mut mol = Molecule::new("source".to_string());
        mol.atoms = vec![atom('A', 10), atom('B', 20), atom('B', 30)];
        mol.coord_sets = vec![vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]];
        mol.bonds = vec![
            BondInfo {
                atom_a: 0,
                atom_b: 1,
                order: 1,
            },
            BondInfo {
                atom_a: 1,
                atom_b: 2,
                order: 1,
            },
        ];
        mol.build_residues();

        let mut scene = Scene::default();
        scene.add_molecule(mol);

        let created =
            scene.create_object_from_selection("picked", &Selector::Chain('B'), 1, Some(1), true);

        assert_eq!(created, 2);
        assert_eq!(scene.molecules.len(), 2);
        assert_eq!(scene.molecules[0].atoms.len(), 1);
        assert_eq!(scene.molecules[0].coord_sets[0], vec![[0.0, 0.0, 0.0]]);
        assert!(scene.molecules[0].bonds.is_empty());
        assert_eq!(scene.molecules[1].name, "picked");
        assert_eq!(scene.molecules[1].atoms.len(), 2);
        assert_eq!(
            scene.molecules[1].coord_sets,
            vec![vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]]
        );
        assert_eq!(scene.molecules[1].bonds.len(), 1);
    }

    #[test]
    fn object_names_filters_enabled_and_selection() {
        let mut first = Molecule::new("first".to_string());
        first.atoms = vec![atom('A', 1)];
        first.coord_sets = vec![vec![[0.0, 0.0, 0.0]]];

        let mut second = Molecule::new("second".to_string());
        second.atoms = vec![atom('B', 2)];
        second.coord_sets = vec![vec![[1.0, 0.0, 0.0]]];
        second.visible = false;

        let mut scene = Scene::default();
        scene.molecules.push(first);
        scene.molecules.push(second);

        assert_eq!(
            scene.object_names(false, None, 1),
            vec!["first".to_string(), "second".to_string()]
        );
        assert_eq!(scene.object_names(true, None, 1), vec!["first".to_string()]);
        assert_eq!(
            scene.object_names(false, Some(&Selector::Chain('B')), 1),
            vec!["second".to_string()]
        );
    }

    #[test]
    fn object_type_reports_molecule_objects() {
        let mut scene = Scene::default();
        scene.molecules.push(Molecule::new("mol1".to_string()));

        assert_eq!(scene.object_type("mol1"), Some("object:molecule"));
        assert_eq!(scene.object_type("missing"), None);
    }

    #[test]
    fn object_names_of_type_reports_molecule_objects() {
        let mut scene = Scene::default();
        scene.molecules.push(Molecule::new("mol1".to_string()));
        scene.molecules.push(Molecule::new("mol2".to_string()));

        assert_eq!(
            scene.object_names_of_type("object:molecule"),
            vec!["mol1".to_string(), "mol2".to_string()]
        );
        assert!(scene.object_names_of_type("object:map").is_empty());
    }

    #[test]
    fn closest_distance_uses_requested_state_and_excludes_same_atom() {
        let mut mol = Molecule::new("distance".to_string());
        mol.atoms = vec![atom('A', 1), atom('B', 2), atom('B', 3)];
        mol.coord_sets = vec![
            vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [10.0, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [1.5, 0.0, 0.0]],
        ];

        let mut scene = Scene::default();
        scene.add_molecule(mol);

        let first_state = scene
            .closest_distance(&Selector::Chain('A'), &Selector::Chain('B'), 1)
            .unwrap();
        assert_eq!(first_state.distance, 3.0);

        let second_state = scene
            .closest_distance(&Selector::Chain('A'), &Selector::Chain('B'), 2)
            .unwrap();
        assert_eq!(second_state.distance, 1.5);

        let all = scene
            .closest_distance(&Selector::All, &Selector::All, 1)
            .unwrap();
        assert_eq!(all.distance, 3.0);
        assert!(scene
            .closest_distance(&Selector::Chain('A'), &Selector::None, 1)
            .is_none());
    }

    #[test]
    fn selection_angle_uses_single_points_and_state() {
        let mut mol = Molecule::new("angle".to_string());
        mol.atoms = vec![atom('A', 1), atom('B', 2), atom('C', 3)];
        mol.coord_sets = vec![
            vec![[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [-1.0, 0.0, 0.0]],
        ];

        let mut scene = Scene::default();
        scene.add_molecule(mol);

        let first = scene
            .selection_angle(
                &Selector::Chain('A'),
                &Selector::Chain('B'),
                &Selector::Chain('C'),
                1,
            )
            .unwrap();
        assert!((first - 90.0).abs() < 1e-5);

        let second = scene
            .selection_angle(
                &Selector::Chain('A'),
                &Selector::Chain('B'),
                &Selector::Chain('C'),
                2,
            )
            .unwrap();
        assert!((second - 180.0).abs() < 1e-5);

        assert_eq!(
            scene.selection_angle(
                &Selector::All,
                &Selector::Chain('B'),
                &Selector::Chain('C'),
                1
            ),
            Err(SelectionPointError::Multiple)
        );
        assert_eq!(
            scene.selection_angle(
                &Selector::None,
                &Selector::Chain('B'),
                &Selector::Chain('C'),
                1
            ),
            Err(SelectionPointError::Empty)
        );
    }

    #[test]
    fn selection_dihedral_uses_single_points_and_state() {
        let mut mol = Molecule::new("dihedral".to_string());
        mol.atoms = vec![atom('A', 1), atom('B', 2), atom('C', 3), atom('D', 4)];
        mol.coord_sets = vec![
            vec![
                [1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 1.0, 1.0],
            ],
            vec![
                [1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 1.0, -1.0],
            ],
        ];

        let mut scene = Scene::default();
        scene.add_molecule(mol);

        let first = scene
            .selection_dihedral(
                &Selector::Chain('A'),
                &Selector::Chain('B'),
                &Selector::Chain('C'),
                &Selector::Chain('D'),
                1,
            )
            .unwrap();
        assert!((first - -90.0).abs() < 1e-5);

        let second = scene
            .selection_dihedral(
                &Selector::Chain('A'),
                &Selector::Chain('B'),
                &Selector::Chain('C'),
                &Selector::Chain('D'),
                2,
            )
            .unwrap();
        assert!((second - 90.0).abs() < 1e-5);

        assert_eq!(
            scene.selection_dihedral(
                &Selector::All,
                &Selector::Chain('B'),
                &Selector::Chain('C'),
                &Selector::Chain('D'),
                1,
            ),
            Err(SelectionPointError::Multiple)
        );
        assert_eq!(
            scene.selection_dihedral(
                &Selector::None,
                &Selector::Chain('B'),
                &Selector::Chain('C'),
                &Selector::Chain('D'),
                1,
            ),
            Err(SelectionPointError::Empty)
        );
    }

    #[test]
    fn selection_extent_returns_none_for_empty_selection() {
        let mut mol = Molecule::new("extent".to_string());
        mol.atoms = vec![atom('A', 1)];
        mol.coord_sets = vec![vec![[0.0, 0.0, 0.0]]];

        let mut scene = Scene::default();
        scene.add_molecule(mol);

        assert!(scene.selection_extent(&Selector::None).is_none());
    }

    #[test]
    fn set_object_visibility_matches_names_and_wildcards() {
        let mut scene = Scene::default();
        scene.molecules.push(Molecule::new("protein".to_string()));
        scene.molecules.push(Molecule::new("ligand".to_string()));
        scene
            .molecules
            .push(Molecule::new("protein_copy".to_string()));

        let changed = scene.set_object_visibility("prot*", false);

        assert_eq!(changed, 2);
        assert!(!scene.molecules[0].visible);
        assert!(scene.molecules[1].visible);
        assert!(!scene.molecules[2].visible);
        assert!(scene.geometry_dirty);

        scene.geometry_dirty = false;
        assert_eq!(scene.set_object_visibility("PROTEIN", true), 1);
        assert!(scene.molecules[0].visible);
        assert!(!scene.molecules[2].visible);
        assert!(scene.geometry_dirty);

        scene.geometry_dirty = false;
        assert_eq!(scene.set_object_visibility("all", false), 2);
        assert!(scene.molecules.iter().all(|mol| !mol.visible));
        assert!(scene.geometry_dirty);
    }

    #[test]
    fn set_object_visibility_noop_does_not_dirty_scene() {
        let mut scene = Scene::default();
        scene.molecules.push(Molecule::new("protein".to_string()));

        assert_eq!(scene.set_object_visibility("missing*", false), 0);
        assert!(scene.molecules[0].visible);
        assert!(!scene.geometry_dirty);
    }

    #[test]
    fn delete_objects_matches_names_and_wildcards() {
        let mut scene = Scene::default();
        scene.molecules.push(Molecule::new("protein".to_string()));
        scene.molecules.push(Molecule::new("ligand".to_string()));
        scene
            .molecules
            .push(Molecule::new("protein_copy".to_string()));

        let deleted = scene.delete_objects("prot*");

        assert_eq!(deleted, 2);
        assert_eq!(scene.molecules.len(), 1);
        assert_eq!(scene.molecules[0].name, "ligand");
        assert!(scene.geometry_dirty);
    }

    #[test]
    fn delete_all_clears_molecules_measurements_and_clamps_state() {
        let mut scene = Scene::default();
        let mut mol = Molecule::new("protein".to_string());
        mol.coord_sets = vec![vec![[0.0, 0.0, 0.0]], vec![[1.0, 0.0, 0.0]]];
        scene.molecules.push(mol);
        scene.molecules.push(Molecule::new("ligand".to_string()));
        scene.measurements.push(Measurement {
            p1: [0.0, 0.0, 0.0],
            p2: [1.0, 0.0, 0.0],
            distance: 1.0,
            label: "1.00 A".to_string(),
        });
        scene.current_state = 2;

        let deleted = scene.delete_objects("all");

        assert_eq!(deleted, 2);
        assert!(scene.molecules.is_empty());
        assert!(scene.measurements.is_empty());
        assert_eq!(scene.current_state, 1);
        assert!(scene.geometry_dirty);
    }

    #[test]
    fn delete_objects_noop_does_not_dirty_scene() {
        let mut scene = Scene::default();
        scene.molecules.push(Molecule::new("protein".to_string()));

        assert_eq!(scene.delete_objects("missing*"), 0);
        assert_eq!(scene.molecules.len(), 1);
        assert!(!scene.geometry_dirty);
    }

    #[test]
    fn sort_objects_matches_patterns_and_marks_dirty_on_change() {
        let mut first = Molecule::new("first".to_string());
        first.atoms = vec![atom('B', 2), atom('A', 1)];
        first.coord_sets = vec![vec![[2.0, 0.0, 0.0], [1.0, 0.0, 0.0]]];

        let mut second = Molecule::new("second".to_string());
        second.atoms = vec![atom('B', 2), atom('A', 1)];
        second.coord_sets = vec![vec![[4.0, 0.0, 0.0], [3.0, 0.0, 0.0]]];

        let mut scene = Scene::default();
        scene.add_molecule(first);
        scene.add_molecule(second);
        scene.geometry_dirty = false;

        assert_eq!(scene.sort_objects("first"), 1);
        assert_eq!(
            scene.molecules[0]
                .atoms
                .iter()
                .map(|atom| atom.chain)
                .collect::<Vec<_>>(),
            vec!['A', 'B']
        );
        assert_eq!(
            scene.molecules[1]
                .atoms
                .iter()
                .map(|atom| atom.chain)
                .collect::<Vec<_>>(),
            vec!['B', 'A']
        );
        assert!(scene.geometry_dirty);

        scene.geometry_dirty = false;
        assert_eq!(scene.sort_objects("missing"), 0);
        assert!(!scene.geometry_dirty);
    }

    #[test]
    fn order_objects_moves_matching_names_to_requested_location() {
        let mut scene = Scene::default();
        scene.molecules.push(Molecule::new("gamma".to_string()));
        scene.molecules.push(Molecule::new("alpha".to_string()));
        scene.molecules.push(Molecule::new("beta".to_string()));
        scene.molecules.push(Molecule::new("water".to_string()));

        assert_eq!(
            scene.order_objects("alpha beta", true, OrderLocation::Top),
            2
        );
        assert_eq!(
            scene
                .molecules
                .iter()
                .map(|mol| mol.name.as_str())
                .collect::<Vec<_>>(),
            vec!["alpha", "beta", "gamma", "water"]
        );
        assert!(scene.geometry_dirty);

        scene.geometry_dirty = false;
        assert_eq!(scene.order_objects("w*", false, OrderLocation::Top), 1);
        assert_eq!(
            scene
                .molecules
                .iter()
                .map(|mol| mol.name.as_str())
                .collect::<Vec<_>>(),
            vec!["water", "alpha", "beta", "gamma"]
        );
        assert!(scene.geometry_dirty);

        scene.geometry_dirty = false;
        assert_eq!(
            scene.order_objects("missing", false, OrderLocation::Bottom),
            0
        );
        assert!(!scene.geometry_dirty);
    }

    #[test]
    fn clip_camera_adjusts_planes_and_selection_depth() {
        let mut scene = Scene::default();
        scene.camera.near = 10.0;
        scene.camera.far = 110.0;

        scene.clip_camera(ClipMode::Near, -5.0, None, 1).unwrap();
        assert_eq!(scene.camera.near, 15.0);
        assert_eq!(scene.camera.far, 110.0);

        scene.clip_camera(ClipMode::Far, 10.0, None, 1).unwrap();
        assert_eq!(scene.camera.near, 15.0);
        assert_eq!(scene.camera.far, 100.0);

        scene.clip_camera(ClipMode::Slab, 20.0, None, 1).unwrap();
        assert_eq!(scene.camera.near, 47.5);
        assert_eq!(scene.camera.far, 67.5);

        let mut mol = Molecule::new("clip".to_string());
        mol.atoms = vec![atom('A', 1), atom('A', 2), atom('B', 3)];
        mol.coord_sets = vec![vec![[0.0, 0.0, 20.0], [0.0, 0.0, 10.0], [0.0, 0.0, 0.0]]];
        scene.molecules.push(mol);
        scene
            .clip_camera(ClipMode::Atoms, 5.0, Some(&Selector::Chain('A')), 1)
            .unwrap();
        assert_eq!(scene.camera.near, 25.0);
        assert_eq!(scene.camera.far, 45.0);

        assert_eq!(
            scene.clip_camera(ClipMode::Atoms, 5.0, Some(&Selector::None), 1),
            Err(ClipError::EmptySelection)
        );
    }
}
