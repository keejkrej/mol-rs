use std::path::PathBuf;

use glam::{Quat, Vec3};

use crate::core::atom::{REP_ALL, REP_CARTOON, REP_LINES, REP_SPHERES, REP_STICKS};
use crate::io::pdb::{write_pdb, PdbWriteSource};
use crate::scene::scene::{
    AtomId, AtomIndex, ClipError, ClipMode, Measurement, OrderLocation, RenameObjectError,
    SelectionPointError,
};
use crate::selection::{
    evaluate_with_coords, evaluator::count_selected, parse_selection, parser::Selector,
};

use super::MolApp;

impl MolApp {
    /// Parse "rep_name, selection_expr" from a comma-separated argument string.
    /// If no comma, the entire string is the rep name and selection defaults to "all".
    fn parse_rep_selection(args: &str) -> (String, String) {
        if let Some(comma) = args.find(',') {
            let rep = args[..comma].trim().to_lowercase();
            let sel = args[comma + 1..].trim().to_string();
            (rep, sel)
        } else {
            (args.trim().to_lowercase(), String::new())
        }
    }

    fn parse_file_selection(args: &str) -> (String, String) {
        if let Some(comma) = args.find(',') {
            (
                args[..comma].trim().to_string(),
                args[comma + 1..].trim().to_string(),
            )
        } else {
            (args.trim().to_string(), String::new())
        }
    }

    fn parse_copy_args(args: &str) -> Result<(String, String), String> {
        let parts: Vec<&str> = args.split(',').map(str::trim).collect();
        if parts.len() < 2 || parts[0].is_empty() || parts[1].is_empty() {
            return Err("Usage: copy <target>, <source>".to_string());
        }

        Ok((parts[0].to_string(), parts[1].to_string()))
    }

    fn parse_set_name_args(args: &str) -> Result<(String, String), String> {
        let parts: Vec<&str> = args.split(',').map(str::trim).collect();
        if parts.len() < 2 || parts[0].is_empty() || parts[1].is_empty() {
            return Err("Usage: set_name <old_name>, <new_name>".to_string());
        }

        Ok((parts[0].to_string(), parts[1].to_string()))
    }

    fn parse_order_args(args: &str) -> Result<(String, bool, OrderLocation), String> {
        let usage = "Usage: order <names> [, sort [, location]]";
        let parts: Vec<&str> = args.split(',').map(str::trim).collect();
        let names = parts.first().copied().unwrap_or("");
        if names.is_empty() {
            return Err(usage.to_string());
        }

        let mut sort_selected = false;
        let mut sort_seen = false;
        let mut location = OrderLocation::Current;

        for part in parts.iter().skip(1).filter(|part| !part.is_empty()) {
            let (key, value) = if let Some((key, value)) = part.split_once('=') {
                (Some(key.trim().to_ascii_lowercase()), value.trim())
            } else {
                (None, *part)
            };

            match key.as_deref() {
                Some("sort") => {
                    sort_selected = Self::parse_bool_arg(value, usage)?;
                    sort_seen = true;
                }
                Some("location") => {
                    location = Self::parse_order_location(value, usage)?;
                }
                Some(_) => return Err(usage.to_string()),
                None if !sort_seen => match Self::parse_bool_arg(value, usage) {
                    Ok(parsed) => {
                        sort_selected = parsed;
                        sort_seen = true;
                    }
                    Err(_) => {
                        location = Self::parse_order_location(value, usage)?;
                    }
                },
                None => {
                    location = Self::parse_order_location(value, usage)?;
                }
            }
        }

        Ok((names.to_string(), sort_selected, location))
    }

    fn parse_bool_arg(value: &str, usage: &str) -> Result<bool, String> {
        match value.to_ascii_lowercase().as_str() {
            "1" | "on" | "true" | "yes" => Ok(true),
            "0" | "off" | "false" | "no" => Ok(false),
            _ => Err(usage.to_string()),
        }
    }

    fn parse_order_location(value: &str, usage: &str) -> Result<OrderLocation, String> {
        match value.to_ascii_lowercase().as_str() {
            "top" => Ok(OrderLocation::Top),
            "current" => Ok(OrderLocation::Current),
            "bottom" => Ok(OrderLocation::Bottom),
            _ => Err(usage.to_string()),
        }
    }

    fn parse_clip_args(
        args: &str,
        usage: &str,
    ) -> Result<(ClipMode, f32, Option<Selector>, Option<usize>), String> {
        let parts: Vec<&str> = args.split(',').map(str::trim).collect();
        if parts.len() < 2 || parts[0].is_empty() || parts[1].is_empty() {
            return Err(usage.to_string());
        }

        let mode = match parts[0].to_ascii_lowercase().as_str() {
            "near" => ClipMode::Near,
            "far" => ClipMode::Far,
            "move" => ClipMode::Move,
            "slab" => ClipMode::Slab,
            "atoms" => ClipMode::Atoms,
            "near_set" => ClipMode::NearSet,
            "far_set" => ClipMode::FarSet,
            _ => return Err(usage.to_string()),
        };
        let distance = parts[1].parse::<f32>().map_err(|_| usage.to_string())?;
        let selection = parts
            .get(2)
            .filter(|part| !part.is_empty())
            .map(|part| parse_selection(part))
            .transpose()?;
        let state = parts
            .get(3)
            .filter(|part| !part.is_empty())
            .map(|part| part.parse::<usize>().map_err(|_| usage.to_string()))
            .transpose()?;

        Ok((mode, distance, selection, state))
    }

    fn parse_create_args(
        args: &str,
        usage: &str,
    ) -> Result<(String, Selector, Option<usize>), String> {
        let parts: Vec<&str> = args.split(',').map(str::trim).collect();
        if parts.len() < 2 || parts[0].is_empty() || parts[1].is_empty() {
            return Err(usage.to_string());
        }

        let sel = parse_selection(parts[1])?;
        let source_state = parts
            .get(2)
            .filter(|part| !part.is_empty())
            .map(|part| {
                part.parse::<isize>()
                    .map_err(|_| usage.to_string())
                    .and_then(|state| {
                        if state <= 0 {
                            Ok(None)
                        } else {
                            Ok(Some(state as usize))
                        }
                    })
            })
            .transpose()?
            .flatten();

        Ok((parts[0].to_string(), sel, source_state))
    }

    fn parse_count_atoms_args(args: &str) -> Result<(Selector, Option<usize>), String> {
        let parts: Vec<&str> = args.split(',').map(str::trim).collect();
        let sel = parse_selection(parts.first().copied().unwrap_or(""))?;
        let state = parts
            .get(2)
            .filter(|part| !part.is_empty())
            .map(|part| {
                part.parse::<usize>()
                    .map_err(|_| "Usage: count_atoms [selection [, quiet [, state]]]".to_string())
            })
            .transpose()?;

        Ok((sel, state))
    }

    fn parse_selection_state_args(
        args: &str,
        usage: &str,
    ) -> Result<(Selector, Option<usize>), String> {
        let parts: Vec<&str> = args.split(',').map(str::trim).collect();
        let sel = parse_selection(parts.first().copied().unwrap_or(""))?;
        let state = parts
            .get(1)
            .filter(|part| !part.is_empty())
            .map(|part| part.parse::<usize>().map_err(|_| usage.to_string()))
            .transpose()?;

        Ok((sel, state))
    }

    fn parse_selection_quiet_args(args: &str, usage: &str) -> Result<Selector, String> {
        let parts: Vec<&str> = args.split(',').map(str::trim).collect();
        let sel = parse_selection(parts.first().copied().unwrap_or(""))?;
        if let Some(quiet) = parts.get(1).filter(|part| !part.is_empty()) {
            Self::parse_bool_arg(quiet, usage)?;
        }

        Ok(sel)
    }

    fn parse_selection_mode_args(args: &str, usage: &str) -> Result<(Selector, u8), String> {
        let parts: Vec<&str> = args.split(',').map(str::trim).collect();
        let sel = parse_selection(parts.first().copied().unwrap_or(""))?;
        let mode = parts
            .get(1)
            .filter(|part| !part.is_empty())
            .map(|part| {
                part.parse::<u8>()
                    .map_err(|_| usage.to_string())
                    .and_then(|mode| {
                        if mode <= 1 {
                            Ok(mode)
                        } else {
                            Err(usage.to_string())
                        }
                    })
            })
            .transpose()?
            .unwrap_or(0);

        Ok((sel, mode))
    }

    fn parse_two_selection_state_args(
        args: &str,
        usage: &str,
    ) -> Result<(Selector, Selector, Option<usize>), String> {
        let parts: Vec<&str> = args.split(',').map(str::trim).collect();
        if parts.len() < 2 || parts[0].is_empty() || parts[1].is_empty() {
            return Err(usage.to_string());
        }

        let sel1 = parse_selection(parts[0]).map_err(|e| format!("Selection 1 error: {}", e))?;
        let sel2 = parse_selection(parts[1]).map_err(|e| format!("Selection 2 error: {}", e))?;
        let state = parts
            .get(2)
            .filter(|part| !part.is_empty())
            .map(|part| part.parse::<usize>().map_err(|_| usage.to_string()))
            .transpose()?;

        Ok((sel1, sel2, state))
    }

    fn parse_three_selection_state_args(
        args: &str,
        usage: &str,
    ) -> Result<(Selector, Selector, Selector, Option<usize>), String> {
        let parts: Vec<&str> = args.split(',').map(str::trim).collect();
        if parts.len() < 3 || parts[0].is_empty() || parts[1].is_empty() || parts[2].is_empty() {
            return Err(usage.to_string());
        }

        let sel1 = parse_selection(parts[0]).map_err(|e| format!("Selection 1 error: {}", e))?;
        let sel2 = parse_selection(parts[1]).map_err(|e| format!("Selection 2 error: {}", e))?;
        let sel3 = parse_selection(parts[2]).map_err(|e| format!("Selection 3 error: {}", e))?;
        let state = parts
            .get(3)
            .filter(|part| !part.is_empty())
            .map(|part| part.parse::<usize>().map_err(|_| usage.to_string()))
            .transpose()?;

        Ok((sel1, sel2, sel3, state))
    }

    fn parse_four_selection_state_args(
        args: &str,
        usage: &str,
    ) -> Result<(Selector, Selector, Selector, Selector, Option<usize>), String> {
        let parts: Vec<&str> = args.split(',').map(str::trim).collect();
        if parts.len() < 4
            || parts[0].is_empty()
            || parts[1].is_empty()
            || parts[2].is_empty()
            || parts[3].is_empty()
        {
            return Err(usage.to_string());
        }

        let sel1 = parse_selection(parts[0]).map_err(|e| format!("Selection 1 error: {}", e))?;
        let sel2 = parse_selection(parts[1]).map_err(|e| format!("Selection 2 error: {}", e))?;
        let sel3 = parse_selection(parts[2]).map_err(|e| format!("Selection 3 error: {}", e))?;
        let sel4 = parse_selection(parts[3]).map_err(|e| format!("Selection 4 error: {}", e))?;
        let state = parts
            .get(4)
            .filter(|part| !part.is_empty())
            .map(|part| part.parse::<usize>().map_err(|_| usage.to_string()))
            .transpose()?;

        Ok((sel1, sel2, sel3, sel4, state))
    }

    fn parse_get_names_args(args: &str) -> Result<(String, bool, Option<Selector>), String> {
        let parts: Vec<&str> = args.split(',').map(str::trim).collect();
        let name_type = parts
            .first()
            .filter(|part| !part.is_empty())
            .copied()
            .unwrap_or("public_objects")
            .to_ascii_lowercase();
        let enabled_only = parts
            .get(1)
            .filter(|part| !part.is_empty())
            .map(|part| match part.to_ascii_lowercase().as_str() {
                "1" | "on" | "true" => Ok(true),
                "0" | "off" | "false" => Ok(false),
                _ => Err("Usage: get_names [type [, enabled_only [, selection]]]".to_string()),
            })
            .transpose()?
            .unwrap_or(false);
        let selection = parts
            .get(2)
            .filter(|part| !part.is_empty())
            .map(|part| parse_selection(part))
            .transpose()?;

        Ok((name_type, enabled_only, selection))
    }

    fn get_names_includes_objects(name_type: &str) -> Result<bool, String> {
        match name_type {
            "objects"
            | "all"
            | "public"
            | "public_objects"
            | "public_nongroup_objects"
            | "nongroup_objects" => Ok(true),
            "selections" | "public_selections" | "public_group_objects" | "group_objects" => {
                Ok(false)
            }
            _ => Err(format!("Unknown name type: '{}'", name_type)),
        }
    }

    fn parse_type_public_args(args: &str, usage: &str) -> Result<(String, bool), String> {
        let parts: Vec<&str> = args.split(',').map(str::trim).collect();
        let object_type = parts.first().copied().unwrap_or("");
        if object_type.is_empty() {
            return Err(usage.to_string());
        }
        let public = parts
            .get(1)
            .filter(|part| !part.is_empty())
            .map(|part| match part.to_ascii_lowercase().as_str() {
                "1" | "on" | "true" => Ok(true),
                "0" | "off" | "false" => Ok(false),
                _ => Err(usage.to_string()),
            })
            .transpose()?
            .unwrap_or(true);

        Ok((object_type.to_string(), public))
    }

    fn representation_mask(name: &str) -> Option<u32> {
        match name {
            "lines" | "line" | "wire" | "wires" => Some(REP_LINES),
            "sticks" | "stick" => Some(REP_STICKS),
            "spheres" | "sphere" => Some(REP_SPHERES),
            "cartoon" | "ribbon" => Some(REP_CARTOON),
            "everything" | "all" => Some(REP_ALL),
            _ => None,
        }
    }

    fn parse_showhide_selection(args: &str, default_rep: &str) -> (String, String) {
        if args.is_empty() {
            return (default_rep.to_string(), String::new());
        }

        let (rep, sel) = Self::parse_rep_selection(args);
        if sel.is_empty() && rep == "all" {
            (default_rep.to_string(), "all".to_string())
        } else {
            (rep, sel)
        }
    }

    fn parse_color(name: &str) -> Option<[f32; 3]> {
        match name {
            "red" => Some([1.0, 0.2, 0.2]),
            "green" => Some([0.2, 1.0, 0.2]),
            "blue" => Some([0.2, 0.2, 1.0]),
            "yellow" => Some([1.0, 1.0, 0.2]),
            "cyan" => Some([0.2, 1.0, 1.0]),
            "magenta" => Some([1.0, 0.2, 1.0]),
            "orange" => Some([1.0, 0.6, 0.2]),
            "white" => Some([1.0, 1.0, 1.0]),
            "gray" | "grey" => Some([0.5, 0.5, 0.5]),
            "pink" => Some([1.0, 0.65, 0.85]),
            "salmon" => Some([1.0, 0.6, 0.5]),
            "purple" => Some([0.6, 0.2, 0.8]),
            _ => None,
        }
    }

    fn format_chain_list(chains: &[char]) -> String {
        let items: Vec<String> = chains
            .iter()
            .map(|chain| {
                if *chain == '\0' || *chain == ' ' {
                    "''".to_string()
                } else {
                    format!("'{}'", chain)
                }
            })
            .collect();
        format!("[{}]", items.join(", "))
    }

    fn format_string_list(items: &[String]) -> String {
        format!(
            "[{}]",
            items
                .iter()
                .map(|item| format!("'{}'", item))
                .collect::<Vec<_>>()
                .join(", ")
        )
    }

    fn format_atom_index_list(items: &[AtomIndex]) -> String {
        format!(
            "[{}]",
            items
                .iter()
                .map(|item| format!("('{}', {})", item.object, item.index))
                .collect::<Vec<_>>()
                .join(", ")
        )
    }

    fn format_atom_id_list(items: &[AtomId], mode: u8) -> String {
        format!(
            "[{}]",
            items
                .iter()
                .map(|item| {
                    if mode == 1 {
                        format!("('{}', {})", item.object, item.id)
                    } else {
                        item.id.to_string()
                    }
                })
                .collect::<Vec<_>>()
                .join(", ")
        )
    }

    fn format_selection_point_error(error: SelectionPointError) -> &'static str {
        match error {
            SelectionPointError::Empty => "One or more selections are empty.",
            SelectionPointError::Multiple => "Selections must match exactly one atom.",
            SelectionPointError::Degenerate => "Angle is undefined for coincident atoms.",
        }
    }

    pub(super) fn handle_command(&mut self, cmd: &str) {
        let parts: Vec<&str> = cmd.splitn(2, char::is_whitespace).collect();
        let verb = parts[0].to_lowercase();
        let args = if parts.len() > 1 { parts[1].trim() } else { "" };

        match verb.as_str() {
            "load" => {
                if args.is_empty() {
                    self.command_line.log("Usage: load <filename>");
                } else {
                    let path = PathBuf::from(args);
                    self.load_file(path);
                }
            }
            "save" => {
                let (filename, sel_str) = Self::parse_file_selection(args);
                if filename.is_empty() {
                    self.command_line.log("Usage: save <file>[, <selection>]");
                    return;
                }

                let sel = match parse_selection(&sel_str) {
                    Ok(s) => s,
                    Err(e) => {
                        self.command_line.log(format!("Selection error: {}", e));
                        return;
                    }
                };

                let current_state = self.scene.current_state;
                let masks: Vec<Vec<bool>> = self
                    .scene
                    .molecules
                    .iter()
                    .map(|mol| evaluate_with_coords(&sel, mol, mol.coords_for_state(current_state)))
                    .collect();
                let sources: Vec<PdbWriteSource<'_>> = self
                    .scene
                    .molecules
                    .iter()
                    .zip(masks.iter())
                    .map(|(mol, mask)| PdbWriteSource {
                        molecule: mol,
                        coords: mol.coords_for_state(current_state),
                        mask: Some(mask.as_slice()),
                    })
                    .collect();

                match write_pdb(&PathBuf::from(&filename), &sources) {
                    Ok(atom_count) => self
                        .command_line
                        .log(format!("Saved {} atoms to {}", atom_count, filename)),
                    Err(e) => self.command_line.log(e),
                }
            }
            "copy" => {
                let (target, source) = match Self::parse_copy_args(args) {
                    Ok(parsed) => parsed,
                    Err(e) => {
                        self.command_line.log(e);
                        return;
                    }
                };

                if self.scene.copy_object(&target, &source) {
                    self.command_line
                        .log(format!("copy {}, {}: 1 object(s)", target, source));
                } else {
                    self.command_line
                        .log(format!("copy: unknown object '{}'", source));
                }
            }
            "set_name" => {
                let (old_name, new_name) = match Self::parse_set_name_args(args) {
                    Ok(parsed) => parsed,
                    Err(e) => {
                        self.command_line.log(e);
                        return;
                    }
                };

                match self.scene.rename_object(&old_name, &new_name) {
                    Ok(()) => self
                        .command_line
                        .log(format!("set_name {}, {}", old_name, new_name)),
                    Err(RenameObjectError::NotFound) => self
                        .command_line
                        .log(format!("set_name: unknown object '{}'", old_name)),
                    Err(RenameObjectError::NameExists) => self
                        .command_line
                        .log(format!("set_name: object '{}' already exists", new_name)),
                }
            }
            "create" | "extract" => {
                let usage = if verb == "create" {
                    "Usage: create <name>, <selection> [, source_state]"
                } else {
                    "Usage: extract <name>, <selection> [, source_state]"
                };
                let (name, sel, source_state) = match Self::parse_create_args(args, usage) {
                    Ok(parsed) => parsed,
                    Err(e) => {
                        self.command_line.log(e);
                        return;
                    }
                };
                let selection_state = source_state.unwrap_or(self.scene.current_state);
                let extract = verb == "extract";
                let created = self.scene.create_object_from_selection(
                    &name,
                    &sel,
                    selection_state,
                    source_state,
                    extract,
                );
                self.command_line
                    .log(format!("{} {}: {} atoms", verb, name, created));
            }
            "show" => {
                let (rep_name, sel_str) = Self::parse_showhide_selection(args, "wire");
                let flag = match Self::representation_mask(&rep_name) {
                    Some(f) => f,
                    None => {
                        self.command_line.log(format!(
                            "Unknown rep: '{}'. Use wire/lines/sticks/spheres/cartoon/everything",
                            rep_name
                        ));
                        return;
                    }
                };
                let sel = match parse_selection(&sel_str) {
                    Ok(s) => s,
                    Err(e) => {
                        self.command_line.log(format!("Selection error: {}", e));
                        return;
                    }
                };
                let mut total = 0usize;
                let current_state = self.scene.current_state;
                for mol in &mut self.scene.molecules {
                    let mask = evaluate_with_coords(&sel, mol, mol.coords_for_state(current_state));
                    for (i, atom) in mol.atoms.iter_mut().enumerate() {
                        if mask[i] {
                            atom.vis_rep |= flag;
                            total += 1;
                        }
                    }
                }
                self.scene.geometry_dirty = true;
                self.command_line
                    .log(format!("show {}: {} atoms", rep_name, total));
            }
            "hide" => {
                let (rep_name, sel_str) = Self::parse_showhide_selection(args, "everything");
                let flag = match Self::representation_mask(&rep_name) {
                    Some(f) => f,
                    None => {
                        self.command_line.log(format!(
                            "Unknown rep: '{}'. Use wire/lines/sticks/spheres/cartoon/everything",
                            rep_name
                        ));
                        return;
                    }
                };
                let sel = match parse_selection(&sel_str) {
                    Ok(s) => s,
                    Err(e) => {
                        self.command_line.log(format!("Selection error: {}", e));
                        return;
                    }
                };
                let mut total = 0usize;
                let current_state = self.scene.current_state;
                for mol in &mut self.scene.molecules {
                    let mask = evaluate_with_coords(&sel, mol, mol.coords_for_state(current_state));
                    for (i, atom) in mol.atoms.iter_mut().enumerate() {
                        if mask[i] {
                            atom.vis_rep &= !flag;
                            total += 1;
                        }
                    }
                }
                self.scene.geometry_dirty = true;
                self.command_line
                    .log(format!("hide {}: {} atoms", rep_name, total));
            }
            "as" | "show_as" => {
                let (rep_name, sel_str) = Self::parse_showhide_selection(args, "wire");
                let flag = match Self::representation_mask(&rep_name) {
                    Some(f) => f,
                    None => {
                        self.command_line.log(format!(
                            "Unknown rep: '{}'. Use wire/lines/sticks/spheres/cartoon/everything",
                            rep_name
                        ));
                        return;
                    }
                };
                let sel = match parse_selection(&sel_str) {
                    Ok(s) => s,
                    Err(e) => {
                        self.command_line.log(format!("Selection error: {}", e));
                        return;
                    }
                };
                let mut total = 0usize;
                let current_state = self.scene.current_state;
                for mol in &mut self.scene.molecules {
                    let mask = evaluate_with_coords(&sel, mol, mol.coords_for_state(current_state));
                    for (i, atom) in mol.atoms.iter_mut().enumerate() {
                        if mask[i] {
                            atom.vis_rep = flag;
                            total += 1;
                        }
                    }
                }
                self.scene.geometry_dirty = true;
                self.command_line
                    .log(format!("as {}: {} atoms", rep_name, total));
            }
            "color" => {
                // color <color_name>, <selection>
                let (color_name, sel_str) = Self::parse_rep_selection(args);
                let rgb = match Self::parse_color(&color_name) {
                    Some(c) => c,
                    None => {
                        self.command_line.log(format!(
                            "Unknown color: '{}'. Try: red green blue yellow cyan magenta orange white gray pink salmon purple",
                            color_name
                        ));
                        return;
                    }
                };
                let sel = match parse_selection(&sel_str) {
                    Ok(s) => s,
                    Err(e) => {
                        self.command_line.log(format!("Selection error: {}", e));
                        return;
                    }
                };
                let mut total = 0usize;
                let current_state = self.scene.current_state;
                for mol in &mut self.scene.molecules {
                    let mask = evaluate_with_coords(&sel, mol, mol.coords_for_state(current_state));
                    for (i, atom) in mol.atoms.iter_mut().enumerate() {
                        if mask[i] {
                            atom.color = rgb;
                            total += 1;
                        }
                    }
                }
                self.scene.geometry_dirty = true;
                self.command_line
                    .log(format!("color {}: {} atoms", color_name, total));
            }
            "select" => {
                // select <selection> — just counts matching atoms
                let sel = match parse_selection(args) {
                    Ok(s) => s,
                    Err(e) => {
                        self.command_line.log(format!("Selection error: {}", e));
                        return;
                    }
                };
                let mut total = 0usize;
                let current_state = self.scene.current_state;
                for mol in &self.scene.molecules {
                    let mask = evaluate_with_coords(&sel, mol, mol.coords_for_state(current_state));
                    total += count_selected(&mask);
                }
                self.command_line.log(format!("Selected {} atoms", total));
            }
            "count_atoms" => {
                let (sel, state) = match Self::parse_count_atoms_args(args) {
                    Ok(parsed) => parsed,
                    Err(e) => {
                        self.command_line.log(e);
                        return;
                    }
                };
                let state = state.unwrap_or(self.scene.current_state);
                let total = self.scene.count_selected_atoms(&sel, state);
                self.command_line
                    .log(format!("count_atoms: {} atoms", total));
            }
            "count_states" => {
                let sel = match parse_selection(args) {
                    Ok(s) => s,
                    Err(e) => {
                        self.command_line.log(format!("Selection error: {}", e));
                        return;
                    }
                };
                let total = self.scene.count_selection_states(&sel);
                self.command_line
                    .log(format!("count_states: {} states", total));
            }
            "get_extent" => {
                let (sel, state) = match Self::parse_selection_state_args(
                    args,
                    "Usage: get_extent [selection [, state]]",
                ) {
                    Ok(parsed) => parsed,
                    Err(e) => {
                        self.command_line.log(e);
                        return;
                    }
                };
                let state = state.unwrap_or(self.scene.current_state);
                let Some(bounds) = self.scene.selection_bounds(&sel, state) else {
                    self.command_line.log("No atoms in selection.");
                    return;
                };

                self.command_line.log(format!(
                    "cmd.extent: min: [{:8.3},{:8.3},{:8.3}]",
                    bounds.min[0], bounds.min[1], bounds.min[2]
                ));
                self.command_line.log(format!(
                    "cmd.extent: max: [{:8.3},{:8.3},{:8.3}]",
                    bounds.max[0], bounds.max[1], bounds.max[2]
                ));
            }
            "get_atom_coords" => {
                let (sel, state) = match Self::parse_selection_state_args(
                    args,
                    "Usage: get_atom_coords <selection> [, state]",
                ) {
                    Ok(parsed) => parsed,
                    Err(e) => {
                        self.command_line.log(e);
                        return;
                    }
                };
                let state = state.unwrap_or(self.scene.current_state);
                let coords = self.scene.selection_coords(&sel, state);
                match coords.as_slice() {
                    [coord] => {
                        self.command_line.log(format!(
                            "cmd.get_atom_coords: [{:8.3},{:8.3},{:8.3}]",
                            coord[0], coord[1], coord[2]
                        ));
                    }
                    [] => self.command_line.log("No atoms in selection."),
                    _ => self
                        .command_line
                        .log("get_atom_coords: selection must match exactly one atom."),
                }
            }
            "get_chains" => {
                let (sel, state) = match Self::parse_selection_state_args(
                    args,
                    "Usage: get_chains [selection [, state]]",
                ) {
                    Ok(parsed) => parsed,
                    Err(e) => {
                        self.command_line.log(e);
                        return;
                    }
                };
                let state = state.unwrap_or(self.scene.current_state);
                let chains = self.scene.selection_chains(&sel, state);
                self.command_line.log(format!(
                    "cmd.get_chains: {}",
                    Self::format_chain_list(&chains)
                ));
            }
            "index" => {
                let sel = match parse_selection(args) {
                    Ok(s) => s,
                    Err(e) => {
                        self.command_line.log(format!("Selection error: {}", e));
                        return;
                    }
                };
                let indices = self.scene.selection_indices(&sel, self.scene.current_state);
                self.command_line.log(format!(
                    "cmd.index: {}",
                    Self::format_atom_index_list(&indices)
                ));
            }
            "identify" => {
                let (sel, mode) = match Self::parse_selection_mode_args(
                    args,
                    "Usage: identify [selection [, mode]]",
                ) {
                    Ok(parsed) => parsed,
                    Err(e) => {
                        self.command_line.log(e);
                        return;
                    }
                };
                let ids = self.scene.selection_ids(&sel, self.scene.current_state);
                self.command_line.log(format!(
                    "cmd.identify: {}",
                    Self::format_atom_id_list(&ids, mode)
                ));
            }
            "get_names" => {
                let (name_type, enabled_only, selection) = match Self::parse_get_names_args(args) {
                    Ok(parsed) => parsed,
                    Err(e) => {
                        self.command_line.log(e);
                        return;
                    }
                };
                let include_objects = match Self::get_names_includes_objects(&name_type) {
                    Ok(include_objects) => include_objects,
                    Err(e) => {
                        self.command_line.log(e);
                        return;
                    }
                };
                let names = if include_objects {
                    self.scene.object_names(
                        enabled_only,
                        selection.as_ref(),
                        self.scene.current_state,
                    )
                } else {
                    Vec::new()
                };
                self.command_line.log(format!(
                    "cmd.get_names: {}",
                    Self::format_string_list(&names)
                ));
            }
            "get_object_list" => {
                let sel = match Self::parse_selection_quiet_args(
                    args,
                    "Usage: get_object_list [selection [, quiet]]",
                ) {
                    Ok(sel) => sel,
                    Err(e) => {
                        self.command_line.log(e);
                        return;
                    }
                };
                let names = self
                    .scene
                    .object_names(false, Some(&sel), self.scene.current_state);
                self.command_line.log(format!(
                    "cmd.get_object_list: {}",
                    Self::format_string_list(&names)
                ));
            }
            "get_type" => {
                if args.is_empty() {
                    self.command_line.log("Usage: get_type <object>");
                    return;
                }
                match self.scene.object_type(args) {
                    Some(object_type) => self
                        .command_line
                        .log(format!("cmd.get_type: {}", object_type)),
                    None => self
                        .command_line
                        .log(format!("get_type: unknown object '{}'", args)),
                }
            }
            "get_names_of_type" => {
                let (object_type, _public) = match Self::parse_type_public_args(
                    args,
                    "Usage: get_names_of_type <type> [, public]",
                ) {
                    Ok(parsed) => parsed,
                    Err(e) => {
                        self.command_line.log(e);
                        return;
                    }
                };
                let names = self.scene.object_names_of_type(&object_type);
                self.command_line.log(format!(
                    "cmd.get_names_of_type: {}",
                    Self::format_string_list(&names)
                ));
            }
            "get_distance" => {
                let (sel1, sel2, state) = match Self::parse_two_selection_state_args(
                    args,
                    "Usage: get_distance <sel1>, <sel2> [, state]",
                ) {
                    Ok(parsed) => parsed,
                    Err(e) => {
                        self.command_line.log(e);
                        return;
                    }
                };
                let state = state.unwrap_or(self.scene.current_state);
                let Some(result) = self.scene.closest_distance(&sel1, &sel2, state) else {
                    self.command_line.log("One or both selections are empty.");
                    return;
                };
                self.command_line.log(format!(
                    "cmd.get_distance: {:.3} Angstroms.",
                    result.distance
                ));
            }
            "get_angle" => {
                let (sel1, sel2, sel3, state) = match Self::parse_three_selection_state_args(
                    args,
                    "Usage: get_angle <sel1>, <sel2>, <sel3> [, state]",
                ) {
                    Ok(parsed) => parsed,
                    Err(e) => {
                        self.command_line.log(e);
                        return;
                    }
                };
                let state = state.unwrap_or(self.scene.current_state);
                match self.scene.selection_angle(&sel1, &sel2, &sel3, state) {
                    Ok(angle) => self
                        .command_line
                        .log(format!("cmd.get_angle: {:.3} degrees.", angle)),
                    Err(e) => self.command_line.log(Self::format_selection_point_error(e)),
                }
            }
            "get_dihedral" => {
                let (sel1, sel2, sel3, sel4, state) = match Self::parse_four_selection_state_args(
                    args,
                    "Usage: get_dihedral <sel1>, <sel2>, <sel3>, <sel4> [, state]",
                ) {
                    Ok(parsed) => parsed,
                    Err(e) => {
                        self.command_line.log(e);
                        return;
                    }
                };
                let state = state.unwrap_or(self.scene.current_state);
                match self
                    .scene
                    .selection_dihedral(&sel1, &sel2, &sel3, &sel4, state)
                {
                    Ok(angle) => self
                        .command_line
                        .log(format!("cmd.get_dihedral: {:.3} degrees.", angle)),
                    Err(e) => self.command_line.log(Self::format_selection_point_error(e)),
                }
            }
            "get_position" => {
                let target = self.scene.camera.target;
                self.command_line.log(format!(
                    "cmd.get_position: [{:8.3},{:8.3},{:8.3}]",
                    target.x, target.y, target.z
                ));
            }
            "get_clip" => {
                self.command_line.log(format!(
                    "cmd.get_clip: [{:.3}, {:.3}]",
                    self.scene.camera.near, self.scene.camera.far
                ));
            }
            "remove" | "rm" => {
                if args.is_empty() {
                    self.command_line.log("Usage: remove <selection>");
                    return;
                }

                let sel = match parse_selection(args) {
                    Ok(s) => s,
                    Err(e) => {
                        self.command_line.log(format!("Selection error: {}", e));
                        return;
                    }
                };

                let current_state = self.scene.current_state;
                let mut total = 0usize;
                for mol in &mut self.scene.molecules {
                    let mask = evaluate_with_coords(&sel, mol, mol.coords_for_state(current_state));
                    total += mol.remove_atoms(&mask);
                }

                if total > 0 {
                    self.scene.geometry_dirty = true;
                }
                self.command_line.log(format!("Removed {} atoms", total));
            }
            "zoom" | "orient" => {
                let sel = match parse_selection(args) {
                    Ok(s) => s,
                    Err(e) => {
                        self.command_line.log(format!("Selection error: {}", e));
                        return;
                    }
                };

                let Some(extent) = self.scene.selection_extent(&sel) else {
                    self.command_line.log("No atoms in selection.");
                    return;
                };

                self.scene.camera.reset_to_fit(extent.center, extent.radius);
                if verb == "orient" {
                    self.scene.camera.rotation = Quat::IDENTITY;
                }
                self.command_line
                    .log(format!("{}: {} atoms", verb, extent.atom_count));
            }
            "center" => {
                let sel = match parse_selection(args) {
                    Ok(s) => s,
                    Err(e) => {
                        self.command_line.log(format!("Selection error: {}", e));
                        return;
                    }
                };

                let Some(extent) = self.scene.selection_extent(&sel) else {
                    self.command_line.log("No atoms in selection.");
                    return;
                };

                self.scene.camera.target = Vec3::from_array(extent.center);
                self.command_line
                    .log(format!("center: {} atoms", extent.atom_count));
            }
            "enable" | "disable" => {
                let visible = verb == "enable";
                let pattern = if args.is_empty() { "all" } else { args };
                let changed = self.scene.set_object_visibility(pattern, visible);
                self.command_line
                    .log(format!("{} {}: {} object(s)", verb, pattern, changed));
            }
            "delete" | "del" => {
                if args.is_empty() {
                    self.command_line.log("Usage: delete <object|all>");
                    return;
                }

                let deleted = self.scene.delete_objects(args);
                self.command_line
                    .log(format!("delete {}: {} object(s)", args, deleted));
            }
            "sort" => {
                let pattern = if args.is_empty() { "all" } else { args };
                let sorted = self.scene.sort_objects(pattern);
                self.command_line
                    .log(format!("sort {}: {} object(s)", pattern, sorted));
            }
            "order" => {
                let (names, sort_selected, location) = match Self::parse_order_args(args) {
                    Ok(parsed) => parsed,
                    Err(e) => {
                        self.command_line.log(e);
                        return;
                    }
                };
                let ordered = self.scene.order_objects(&names, sort_selected, location);
                self.command_line
                    .log(format!("order {}: {} object(s)", names, ordered));
            }
            "clip" => {
                let usage = "Usage: clip <near|far|move|slab|atoms|near_set|far_set>, <distance> [, selection [, state]]";
                let (mode, distance, selection, state) = match Self::parse_clip_args(args, usage) {
                    Ok(parsed) => parsed,
                    Err(e) => {
                        self.command_line.log(e);
                        return;
                    }
                };
                let state = state.unwrap_or(self.scene.current_state);
                match self
                    .scene
                    .clip_camera(mode, distance, selection.as_ref(), state)
                {
                    Ok(()) => self.command_line.log(format!(
                        "clip: near {:.3}, far {:.3}",
                        self.scene.camera.near, self.scene.camera.far
                    )),
                    Err(ClipError::EmptySelection) => {
                        self.command_line.log("No atoms in selection.")
                    }
                }
            }
            "distance" | "dist" => {
                let (sel1, sel2, state) = match Self::parse_two_selection_state_args(
                    args,
                    "Usage: distance <sel1>, <sel2>",
                ) {
                    Ok(parsed) => parsed,
                    Err(e) => {
                        self.command_line.log(e);
                        return;
                    }
                };
                let state = state.unwrap_or(self.scene.current_state);
                let Some(result) = self.scene.closest_distance(&sel1, &sel2, state) else {
                    self.command_line.log("One or both selections are empty.");
                    return;
                };

                self.scene.measurements.push(Measurement {
                    p1: result.p1,
                    p2: result.p2,
                    distance: result.distance,
                    label: format!("{:.2} Å", result.distance),
                });
                self.command_line
                    .log(format!("Distance: {:.2} Å", result.distance));
            }
            "png" => {
                if args.is_empty() {
                    self.command_line.log("Usage: png <filename>");
                } else {
                    self.screenshot_requested = Some(PathBuf::from(args));
                    self.command_line.log("Screenshot requested...");
                }
            }
            "reset" => {
                if let Some(mol) = self.scene.molecules.first() {
                    let c = mol.centroid_for_state(self.scene.current_state);
                    let r = mol.radius_for_state(self.scene.current_state);
                    self.scene.camera.reset_to_fit(c, r);
                }
                self.scene.measurements.clear();
                self.command_line.log("View reset.");
            }
            "state" | "frame" => {
                if args.is_empty() {
                    self.command_line.log(format!(
                        "State: {}/{}",
                        self.scene.current_state,
                        self.scene.max_state_count()
                    ));
                    return;
                }
                match args.parse::<usize>() {
                    Ok(n) => {
                        self.scene.set_state_clamped(n);
                        self.command_line.log(format!(
                            "State: {}/{}",
                            self.scene.current_state,
                            self.scene.max_state_count()
                        ));
                    }
                    Err(_) => {
                        self.command_line.log("Usage: state <n>");
                    }
                }
            }
            "next" => {
                self.scene.next_state();
                self.command_line.log(format!(
                    "State: {}/{}",
                    self.scene.current_state,
                    self.scene.max_state_count()
                ));
            }
            "prev" => {
                self.scene.prev_state();
                self.command_line.log(format!(
                    "State: {}/{}",
                    self.scene.current_state,
                    self.scene.max_state_count()
                ));
            }
            "all_states" => {
                if args.is_empty() {
                    self.command_line.log(format!(
                        "all_states is {}",
                        if self.scene.all_states { "on" } else { "off" }
                    ));
                    return;
                }
                let v = args.to_lowercase();
                let parsed = match v.as_str() {
                    "on" | "1" | "true" => Some(true),
                    "off" | "0" | "false" => Some(false),
                    _ => None,
                };
                if let Some(on) = parsed {
                    if self.scene.all_states != on {
                        self.scene.all_states = on;
                        self.scene.geometry_dirty = true;
                    }
                    self.command_line.log(format!(
                        "all_states {}",
                        if self.scene.all_states { "on" } else { "off" }
                    ));
                } else {
                    self.command_line.log("Usage: all_states <on|off>");
                }
            }
            "bg_color" | "bg" | "bgcolor" => {
                let color_name = args.trim().to_lowercase();
                if let Some(rgb) = Self::parse_color(&color_name) {
                    self.scene.bg_color = rgb;
                    self.command_line
                        .log(format!("Background set to {}", color_name));
                } else {
                    self.command_line
                        .log(format!("Unknown color: '{}'", color_name));
                }
            }
            "help" => {
                self.command_line.log("Commands:");
                self.command_line
                    .log("  load <file>             — Load a PDB/CIF file");
                self.command_line
                    .log("  save <file>[, <sel>]    — Save current state as PDB");
                self.command_line
                    .log("  copy <target>, <source> — Duplicate a molecule object");
                self.command_line
                    .log("  set_name <old>, <new>  — Rename a molecule object");
                self.command_line
                    .log("  create <name>, <sel>    — Create object from selection");
                self.command_line
                    .log("  extract <name>, <sel>   — Move selection into new object");
                self.command_line.log(
                    "  show <rep>[, <sel>]     — Show representation (wire/lines/sticks/spheres/cartoon/everything)",
                );
                self.command_line
                    .log("  hide <rep>[, <sel>]     — Hide representation");
                self.command_line
                    .log("  as <rep>[, <sel>]       — Replace visible representation");
                self.command_line
                    .log("  color <color>[, <sel>]  — Color atoms");
                self.command_line
                    .log("  select <sel>            — Count matching atoms");
                self.command_line
                    .log("  count_atoms [sel]       — Count matching atoms");
                self.command_line
                    .log("  count_states [sel]      — Count states with matching atoms");
                self.command_line
                    .log("  get_extent [sel]        — Print selection min/max coordinates");
                self.command_line
                    .log("  get_atom_coords <sel>   — Print one atom coordinate");
                self.command_line
                    .log("  get_chains [sel]        — Print selected chain IDs");
                self.command_line
                    .log("  index [sel]             — Print object atom indices");
                self.command_line
                    .log("  identify [sel]          — Print selected atom IDs");
                self.command_line
                    .log("  get_names [type]        — Print object names");
                self.command_line
                    .log("  get_object_list [sel]   — Print objects covered by selection");
                self.command_line
                    .log("  get_type <object>       — Print object type");
                self.command_line
                    .log("  get_names_of_type <type> — Print objects of a type");
                self.command_line
                    .log("  get_distance <s1>, <s2> — Query closest distance");
                self.command_line
                    .log("  get_angle <s1>, <s2>, <s3> — Query atom angle");
                self.command_line
                    .log("  get_dihedral <s1>, <s2>, <s3>, <s4> — Query torsion");
                self.command_line
                    .log("  get_position          — Print viewer center");
                self.command_line
                    .log("  get_clip              — Print clipping planes");
                self.command_line
                    .log("  remove <sel>            — Delete matching atoms");
                self.command_line
                    .log("  zoom|center|orient [sel] — Move view to selection");
                self.command_line
                    .log("  enable|disable [obj]    — Toggle object visibility");
                self.command_line
                    .log("  delete <obj|all>        — Delete matching objects");
                self.command_line
                    .log("  sort [obj]              — Sort atoms in molecule objects");
                self.command_line
                    .log("  order <names>[, sort]   — Reorder molecule objects");
                self.command_line
                    .log("  clip <mode>, <dist>     — Adjust clipping planes");
                self.command_line
                    .log("  distance <s1>, <s2>     — Measure distance");
                self.command_line
                    .log("  state <n> / frame <n>   — Set current state");
                self.command_line
                    .log("  next / prev             — Step through states");
                self.command_line
                    .log("  all_states <on|off>     — Render all states or current state");
                self.command_line
                    .log("  png <file>              — Save screenshot");
                self.command_line
                    .log("  bg_color | bgcolor | bg <color> — Set background color");
                self.command_line
                    .log("  reset                   — Reset camera view");
                self.command_line.log(
                    "Selections: chain A+B, model/object obj*, rep wire/sticks/spheres/cartoon/everything, color red, resi 1-50 (resi/residue), resn/resname ALA, name/ca C*, bb/backbone, sidechain, elem/element/symbol C, alt A, ss H/S/L, serial/index/ID 5, b/q/x/y/z <value>, byres/byca/bychain/byobject/bymolecule chain A, first/last chain A, within/around/expand/beyond/near_to 5 chain A, extend 2 chain A, neighbor chain A, bonded, guide, het/hetatm, hydro/hydrogens, solvent, polymer/protein/nucleic, metals, organic, inorganic, visible/vis, present/pr., enabled, all/*, !/not, &/and, |/or, ()",
                );
            }
            _ => {
                self.command_line.log(format!(
                    "Unknown command: '{}'. Type 'help' for usage.",
                    verb
                ));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use glam::Vec3;

    use crate::core::atom::{AtomInfo, REP_ALL, REP_CARTOON, REP_LINES, REP_SPHERES, REP_STICKS};
    use crate::core::bond::BondInfo;
    use crate::core::molecule::Molecule;
    use crate::scene::scene::Scene;
    use crate::ui::command_line::CommandLine;
    use crate::ui::control_panel::ControlPanelState;

    use super::MolApp;

    fn test_app() -> MolApp {
        let mut mol = Molecule::new("obj".to_string());
        mol.atoms.push(AtomInfo {
            name: "CA".to_string(),
            elem_symbol: "C".to_string(),
            chain: 'A',
            vis_rep: 0,
            ..AtomInfo::default()
        });
        mol.atoms.push(AtomInfo {
            name: "O".to_string(),
            elem_symbol: "O".to_string(),
            chain: 'B',
            vis_rep: 0,
            ..AtomInfo::default()
        });
        mol.coord_sets = vec![vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]];

        let mut scene = Scene::default();
        scene.molecules.push(mol);

        MolApp {
            scene,
            command_line: CommandLine::default(),
            control_panel_state: ControlPanelState::default(),
            render_state: None,
            renderer: None,
            offscreen: None,
            open_file_requested: false,
            pending_file: None,
            screenshot_requested: None,
        }
    }

    #[test]
    fn show_hide_accept_pymol_representation_aliases() {
        let mut app = test_app();

        app.handle_command("show");
        assert_eq!(app.scene.molecules[0].atoms[0].vis_rep, REP_LINES);
        assert_eq!(app.scene.molecules[0].atoms[1].vis_rep, REP_LINES);

        app.handle_command("hide");
        assert_eq!(app.scene.molecules[0].atoms[0].vis_rep, 0);
        assert_eq!(app.scene.molecules[0].atoms[1].vis_rep, 0);

        app.handle_command("show everything, all");
        assert_eq!(app.scene.molecules[0].atoms[0].vis_rep, REP_ALL);
        assert_eq!(app.scene.molecules[0].atoms[1].vis_rep, REP_ALL);

        app.handle_command("hide wire, name CA");
        assert_eq!(
            app.scene.molecules[0].atoms[0].vis_rep,
            REP_STICKS | REP_SPHERES | REP_CARTOON
        );
        assert_eq!(app.scene.molecules[0].atoms[1].vis_rep, REP_ALL);
    }

    #[test]
    fn show_as_replaces_representations_for_selection() {
        let mut app = test_app();
        app.scene.molecules[0].atoms[0].vis_rep = REP_LINES | REP_SPHERES;
        app.scene.molecules[0].atoms[1].vis_rep = REP_STICKS;

        app.handle_command("as sticks, name CA");
        assert_eq!(app.scene.molecules[0].atoms[0].vis_rep, REP_STICKS);
        assert_eq!(app.scene.molecules[0].atoms[1].vis_rep, REP_STICKS);

        app.handle_command("show_as ribbon, name O");
        assert_eq!(app.scene.molecules[0].atoms[0].vis_rep, REP_STICKS);
        assert_eq!(app.scene.molecules[0].atoms[1].vis_rep, REP_CARTOON);
    }

    #[test]
    fn count_commands_report_atoms_and_states() {
        let mut app = test_app();
        app.scene.molecules[0].coord_sets = vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0]],
            Vec::new(),
        ];
        app.scene.current_state = 2;

        app.handle_command("count_atoms present");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "count_atoms: 1 atoms"
        );

        app.handle_command("count_atoms present, 1, 1");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "count_atoms: 2 atoms"
        );

        app.handle_command("count_states present");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "count_states: 2 states"
        );

        app.handle_command("count_states name O");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "count_states: 1 states"
        );
    }

    #[test]
    fn get_extent_reports_selection_bounds() {
        let mut app = test_app();
        app.scene.molecules[0].coord_sets = vec![
            vec![[-1.0, 2.0, 0.5], [4.0, -3.0, 2.0]],
            vec![[10.0, 0.0, 0.0], [12.0, 3.0, -2.0]],
        ];
        app.scene.current_state = 2;

        app.handle_command("get_extent all");
        let len = app.command_line.output.len();
        assert_eq!(
            app.command_line.output[len - 2],
            "cmd.extent: min: [  10.000,   0.000,  -2.000]"
        );
        assert_eq!(
            app.command_line.output[len - 1],
            "cmd.extent: max: [  12.000,   3.000,   0.000]"
        );

        app.handle_command("get_extent name CA, 1");
        let len = app.command_line.output.len();
        assert_eq!(
            app.command_line.output[len - 2],
            "cmd.extent: min: [  -1.000,   2.000,   0.500]"
        );
        assert_eq!(
            app.command_line.output[len - 1],
            "cmd.extent: max: [  -1.000,   2.000,   0.500]"
        );

        app.handle_command("get_extent name ZZ");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "No atoms in selection."
        );
    }

    #[test]
    fn get_atom_coords_reports_single_selected_coordinate() {
        let mut app = test_app();
        app.scene.molecules[0].coord_sets = vec![
            vec![[-1.0, 2.0, 0.5], [4.0, -3.0, 2.0]],
            vec![[10.0, 0.0, 0.0], [12.0, 3.0, -2.0]],
        ];
        app.scene.current_state = 2;

        app.handle_command("get_atom_coords name O");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.get_atom_coords: [  12.000,   3.000,  -2.000]"
        );

        app.handle_command("get_atom_coords name CA, 1");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.get_atom_coords: [  -1.000,   2.000,   0.500]"
        );

        app.handle_command("get_atom_coords all");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "get_atom_coords: selection must match exactly one atom."
        );

        app.handle_command("get_atom_coords name ZZ");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "No atoms in selection."
        );
    }

    #[test]
    fn get_chains_reports_selected_chain_ids() {
        let mut app = test_app();
        app.scene.molecules[0].coord_sets = vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec![[2.0, 0.0, 0.0]],
        ];
        app.scene.current_state = 2;

        app.handle_command("get_chains all");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.get_chains: ['A']"
        );

        app.handle_command("get_chains all, 1");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.get_chains: ['A', 'B']"
        );

        app.handle_command("get_chains name ZZ");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.get_chains: []"
        );
    }

    #[test]
    fn index_reports_object_atom_indices() {
        let mut app = test_app();
        app.scene.molecules[0].name = "obj1".to_string();
        app.scene.molecules[0].coord_sets = vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec![[2.0, 0.0, 0.0]],
        ];
        app.scene.current_state = 2;

        app.handle_command("index all");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.index: [('obj1', 1)]"
        );

        app.handle_command("index chain B");
        assert_eq!(app.command_line.output.last().unwrap(), "cmd.index: []");

        app.scene.current_state = 1;
        app.handle_command("index chain B");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.index: [('obj1', 2)]"
        );
    }

    #[test]
    fn identify_reports_atom_ids() {
        let mut app = test_app();
        app.scene.molecules[0].name = "obj1".to_string();
        app.scene.molecules[0].atoms[0].serial = 10;
        app.scene.molecules[0].atoms[1].serial = 20;
        app.scene.molecules[0].coord_sets = vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec![[2.0, 0.0, 0.0]],
        ];
        app.scene.current_state = 2;

        app.handle_command("identify all");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.identify: [10]"
        );

        app.handle_command("identify all, 1");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.identify: [('obj1', 10)]"
        );

        app.scene.current_state = 1;
        app.handle_command("identify chain B");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.identify: [20]"
        );

        app.handle_command("identify all, 2");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "Usage: identify [selection [, mode]]"
        );
    }

    #[test]
    fn copy_create_and_extract_manage_molecule_objects() {
        let mut app = test_app();
        app.scene.molecules[0].name = "source".to_string();
        app.scene.molecules[0].bonds = vec![BondInfo {
            atom_a: 0,
            atom_b: 1,
            order: 1,
        }];
        app.scene.molecules[0].build_residues();

        app.handle_command("copy duplicate, source");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "copy duplicate, source: 1 object(s)"
        );
        assert_eq!(app.scene.molecules.len(), 2);
        assert_eq!(app.scene.molecules[1].name, "duplicate");
        assert_eq!(app.scene.molecules[1].atoms.len(), 2);

        app.handle_command("create just_a, model source and chain A");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "create just_a: 1 atoms"
        );
        assert_eq!(app.scene.molecules.len(), 3);
        assert_eq!(app.scene.molecules[2].name, "just_a");
        assert_eq!(app.scene.molecules[2].atoms[0].chain, 'A');

        app.handle_command("extract moved_b, model source and chain B, 1");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "extract moved_b: 1 atoms"
        );
        assert_eq!(app.scene.molecules.len(), 4);
        assert_eq!(app.scene.molecules[0].atoms.len(), 1);
        assert_eq!(app.scene.molecules[0].atoms[0].chain, 'A');
        assert_eq!(app.scene.molecules[3].name, "moved_b");
        assert_eq!(app.scene.molecules[3].atoms[0].chain, 'B');
    }

    #[test]
    fn set_name_renames_molecule_objects() {
        let mut app = test_app();
        app.scene.molecules[0].name = "source".to_string();

        app.handle_command("set_name source, renamed");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "set_name source, renamed"
        );
        assert_eq!(app.scene.molecules[0].name, "renamed");

        app.handle_command("set_name source, other");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "set_name: unknown object 'source'"
        );
    }

    #[test]
    fn sort_command_reorders_atoms_for_matching_objects() {
        let mut app = test_app();
        app.scene.molecules[0].name = "source".to_string();
        app.scene.molecules[0].atoms[0].chain = 'B';
        app.scene.molecules[0].atoms[0].resi = 2;
        app.scene.molecules[0].atoms[1].chain = 'A';
        app.scene.molecules[0].atoms[1].resi = 1;
        app.scene.molecules[0].coord_sets = vec![vec![[2.0, 0.0, 0.0], [1.0, 0.0, 0.0]]];

        app.handle_command("sort source");

        assert_eq!(
            app.command_line.output.last().unwrap(),
            "sort source: 1 object(s)"
        );
        assert_eq!(
            app.scene.molecules[0]
                .atoms
                .iter()
                .map(|atom| atom.chain)
                .collect::<Vec<_>>(),
            vec!['A', 'B']
        );
        assert_eq!(
            app.scene.molecules[0].coord_sets[0],
            vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
        );
    }

    #[test]
    fn order_command_reorders_molecule_objects() {
        let mut app = test_app();
        app.scene.molecules[0].name = "gamma".to_string();
        app.scene.molecules.push(Molecule::new("alpha".to_string()));
        app.scene.molecules.push(Molecule::new("beta".to_string()));

        app.handle_command("order alpha beta, yes, top");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "order alpha beta: 2 object(s)"
        );
        assert_eq!(
            app.scene
                .molecules
                .iter()
                .map(|mol| mol.name.as_str())
                .collect::<Vec<_>>(),
            vec!["alpha", "beta", "gamma"]
        );

        app.handle_command("order gamma, location=bottom");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "order gamma: 1 object(s)"
        );
        assert_eq!(
            app.scene
                .molecules
                .iter()
                .map(|mol| mol.name.as_str())
                .collect::<Vec<_>>(),
            vec!["alpha", "beta", "gamma"]
        );

        app.handle_command("order");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "Usage: order <names> [, sort [, location]]"
        );
    }

    #[test]
    fn get_names_reports_objects_with_filters() {
        let mut app = test_app();
        app.scene.molecules[0].name = "first".to_string();

        let mut second = Molecule::new("second".to_string());
        second.atoms.push(AtomInfo {
            name: "N".to_string(),
            elem_symbol: "N".to_string(),
            chain: 'C',
            ..AtomInfo::default()
        });
        second.coord_sets = vec![vec![[2.0, 0.0, 0.0]]];
        second.visible = false;
        app.scene.molecules.push(second);

        app.handle_command("get_names");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.get_names: ['first', 'second']"
        );

        app.handle_command("get_names objects, 1");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.get_names: ['first']"
        );

        app.handle_command("get_names objects, 0, chain C");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.get_names: ['second']"
        );

        app.handle_command("get_names selections");
        assert_eq!(app.command_line.output.last().unwrap(), "cmd.get_names: []");
    }

    #[test]
    fn get_object_list_reports_objects_covered_by_selection() {
        let mut app = test_app();
        app.scene.molecules[0].name = "first".to_string();
        app.scene.molecules[0].coord_sets = vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0]],
        ];

        let mut second = Molecule::new("second".to_string());
        second.atoms.push(AtomInfo {
            name: "N".to_string(),
            elem_symbol: "N".to_string(),
            chain: 'B',
            ..AtomInfo::default()
        });
        second.coord_sets = vec![vec![[2.0, 0.0, 0.0]], Vec::new()];
        app.scene.molecules.push(second);
        app.scene.current_state = 2;

        app.handle_command("get_object_list all");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.get_object_list: ['first']"
        );

        app.handle_command("get_object_list chain B, 0");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.get_object_list: []"
        );

        app.scene.current_state = 1;
        app.handle_command("get_object_list chain B");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.get_object_list: ['first', 'second']"
        );

        app.handle_command("get_object_list all, maybe");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "Usage: get_object_list [selection [, quiet]]"
        );
    }

    #[test]
    fn get_type_reports_molecule_object_type() {
        let mut app = test_app();
        app.scene.molecules[0].name = "mol1".to_string();

        app.handle_command("get_type mol1");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.get_type: object:molecule"
        );

        app.handle_command("get_type missing");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "get_type: unknown object 'missing'"
        );

        app.handle_command("get_type");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "Usage: get_type <object>"
        );
    }

    #[test]
    fn get_names_of_type_reports_matching_objects() {
        let mut app = test_app();
        app.scene.molecules[0].name = "mol1".to_string();
        app.scene.molecules.push(Molecule::new("mol2".to_string()));

        app.handle_command("get_names_of_type object:molecule");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.get_names_of_type: ['mol1', 'mol2']"
        );

        app.handle_command("get_names_of_type object:map");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.get_names_of_type: []"
        );

        app.handle_command("get_names_of_type");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "Usage: get_names_of_type <type> [, public]"
        );
    }

    #[test]
    fn get_distance_reports_closest_distance_without_measurement() {
        let mut app = test_app();
        app.scene.molecules[0].coord_sets = vec![
            vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]],
        ];
        app.scene.current_state = 2;

        app.handle_command("get_distance name CA, name O");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.get_distance: 5.000 Angstroms."
        );
        assert!(app.scene.measurements.is_empty());

        app.handle_command("get_distance name CA, name O, 1");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.get_distance: 3.000 Angstroms."
        );
        assert!(app.scene.measurements.is_empty());

        app.handle_command("distance name CA, name O, 1");
        assert_eq!(app.command_line.output.last().unwrap(), "Distance: 3.00 Å");
        assert_eq!(app.scene.measurements.len(), 1);

        app.handle_command("get_distance name ZZ, name O");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "One or both selections are empty."
        );
    }

    #[test]
    fn get_angle_reports_atom_angle() {
        let mut app = test_app();
        app.scene.molecules[0].atoms.push(AtomInfo {
            name: "N".to_string(),
            elem_symbol: "N".to_string(),
            chain: 'C',
            ..AtomInfo::default()
        });
        app.scene.molecules[0].coord_sets = vec![
            vec![[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [-1.0, 0.0, 0.0]],
        ];
        app.scene.current_state = 2;

        app.handle_command("get_angle name CA, name O, name N");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.get_angle: 180.000 degrees."
        );

        app.handle_command("get_angle name CA, name O, name N, 1");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.get_angle: 90.000 degrees."
        );

        app.handle_command("get_angle all, name O, name N");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "Selections must match exactly one atom."
        );

        app.handle_command("get_angle name ZZ, name O, name N");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "One or more selections are empty."
        );
    }

    #[test]
    fn get_dihedral_reports_atom_torsion() {
        let mut app = test_app();
        app.scene.molecules[0].atoms.push(AtomInfo {
            name: "N".to_string(),
            elem_symbol: "N".to_string(),
            chain: 'C',
            ..AtomInfo::default()
        });
        app.scene.molecules[0].atoms.push(AtomInfo {
            name: "CB".to_string(),
            elem_symbol: "C".to_string(),
            chain: 'D',
            ..AtomInfo::default()
        });
        app.scene.molecules[0].coord_sets = vec![
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
        app.scene.current_state = 2;

        app.handle_command("get_dihedral name CA, name O, name N, name CB");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.get_dihedral: 90.000 degrees."
        );

        app.handle_command("get_dihedral name CA, name O, name N, name CB, 1");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.get_dihedral: -90.000 degrees."
        );

        app.handle_command("get_dihedral all, name O, name N, name CB");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "Selections must match exactly one atom."
        );

        app.handle_command("get_dihedral name ZZ, name O, name N, name CB");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "One or more selections are empty."
        );
    }

    #[test]
    fn get_position_reports_camera_target() {
        let mut app = test_app();
        app.scene.camera.target = Vec3::new(1.25, -2.5, 3.75);

        app.handle_command("get_position");

        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.get_position: [   1.250,  -2.500,   3.750]"
        );
    }

    #[test]
    fn clip_commands_adjust_and_report_camera_clip() {
        let mut app = test_app();
        app.scene.camera.near = 10.0;
        app.scene.camera.far = 110.0;

        app.handle_command("clip near, -5");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "clip: near 15.000, far 110.000"
        );

        app.handle_command("clip slab, 20");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "clip: near 52.500, far 72.500"
        );

        app.scene.molecules[0].coord_sets = vec![vec![[0.0, 0.0, 20.0], [0.0, 0.0, 10.0]]];
        app.handle_command("clip atoms, 5, all");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "clip: near 25.000, far 45.000"
        );

        app.handle_command("get_clip");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "cmd.get_clip: [25.000, 45.000]"
        );

        app.handle_command("clip atoms, 5, none");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "No atoms in selection."
        );

        app.handle_command("clip bogus, 1");
        assert_eq!(
            app.command_line.output.last().unwrap(),
            "Usage: clip <near|far|move|slab|atoms|near_set|far_set>, <distance> [, selection [, state]]"
        );
    }
}
