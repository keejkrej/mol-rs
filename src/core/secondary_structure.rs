#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SSType {
    #[default]
    Loop,
    Helix,
    Sheet,
}
