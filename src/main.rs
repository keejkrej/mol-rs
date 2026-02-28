mod app;
mod core;
mod io;
mod render;
mod scene;
mod selection;
mod ui;

use app::MolApp;

#[cfg(target_os = "windows")]
fn configure_wgpu_backend() {
    // On some hybrid AMD/NVIDIA setups, probing all backends can crash in driver layers.
    // Default to DX12 unless the user explicitly selected a backend.
    if std::env::var_os("WGPU_BACKEND").is_none() {
        std::env::set_var("WGPU_BACKEND", "dx12");
    }
}

#[cfg(not(target_os = "windows"))]
fn configure_wgpu_backend() {}

fn main() -> eframe::Result {
    env_logger::init();
    configure_wgpu_backend();

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1280.0, 800.0])
            .with_title("mol — Molecular Viewer"),
        renderer: eframe::Renderer::Wgpu,
        ..Default::default()
    };

    eframe::run_native(
        "mol",
        options,
        Box::new(|cc| Ok(Box::new(MolApp::new(cc)))),
    )
}
