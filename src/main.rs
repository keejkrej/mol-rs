mod app;
mod core;
mod io;
mod render;
mod scene;
mod ui;

use app::MolApp;

fn main() -> eframe::Result {
    env_logger::init();

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1280.0, 800.0])
            .with_title("mol-rs — Molecular Viewer"),
        renderer: eframe::Renderer::Wgpu,
        ..Default::default()
    };

    eframe::run_native(
        "mol-rs",
        options,
        Box::new(|cc| Ok(Box::new(MolApp::new(cc)))),
    )
}
