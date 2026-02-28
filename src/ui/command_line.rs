use egui::Ui;

/// State for the command line input.
pub struct CommandLine {
    pub input: String,
    pub history: Vec<String>,
    pub output: Vec<String>,
}

impl Default for CommandLine {
    fn default() -> Self {
        Self {
            input: String::new(),
            history: Vec::new(),
            output: vec!["mol ready. Type 'help' for commands.".into()],
        }
    }
}

impl CommandLine {
    /// Draw the command line bar. Returns Some(command) if the user pressed Enter.
    pub fn draw(&mut self, ui: &mut Ui) -> Option<String> {
        let mut submitted = None;

        ui.horizontal(|ui| {
            ui.label("PyMOL>");
            let response = ui.text_edit_singleline(&mut self.input);
            if response.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                let cmd = self.input.trim().to_string();
                if !cmd.is_empty() {
                    self.history.push(cmd.clone());
                    submitted = Some(cmd);
                }
                self.input.clear();
                response.request_focus();
            }
        });

        submitted
    }

    pub fn log(&mut self, msg: impl Into<String>) {
        self.output.push(msg.into());
    }

    /// Draw the output log.
    pub fn draw_output(&self, ui: &mut Ui) {
        egui::ScrollArea::vertical()
            .max_height(100.0)
            .stick_to_bottom(true)
            .show(ui, |ui| {
                for line in &self.output {
                    ui.label(line);
                }
            });
    }
}
