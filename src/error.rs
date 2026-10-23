#[derive(thiserror::Error, Debug)]
pub enum DfntransError {
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error("Invalid input file: {0}")]
    InvalidInput(String),
    #[error("Invalid mesh or flow solution: {0}")]
    InvalidNetwork(String),
    #[error(
        "The mesh has no point data array boundary__{0}. Use the name of a boundary__<name> array of the mesh."
    )]
    UnknownBoundary(String),
    #[error("No particle was placed: {0}")]
    NoParticles(String),
    #[error(transparent)]
    Flow(#[from] fehm::FehmError),
}
