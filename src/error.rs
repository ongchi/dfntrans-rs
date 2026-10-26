#[derive(thiserror::Error, Debug)]
pub enum DfntransError {
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error("Invalid input file: cannot read {file}: {source}")]
    ReadInput {
        file: String,
        source: std::io::Error,
    },
    /// The TOML text: syntax, an unknown or a missing key, a wrong type
    #[error("Invalid input file: {0}")]
    InvalidToml(String),
    #[error("Invalid input file: porosity must be above 0 and at most 1")]
    Porosity,
    #[error("Invalid input file: max_time_steps must be above 0")]
    MaxTimeSteps,
    #[error("Invalid input file: inflow and outflow must each name at least one boundary")]
    InflowOutflowMissing,
    #[error("Invalid input file: the boundary {boundary} is in inflow and in outflow")]
    BoundaryInInflowAndOutflow { boundary: String },
    #[error("Invalid input file: the number of particles of seeding must be above 0")]
    ParticleCount,
    #[error("Invalid input file: tdrw.matrix_porosity must be above 0 and at most 1")]
    MatrixPorosity,
    #[error("Invalid input file: tdrw.matrix_diffusivity must be above 0")]
    MatrixDiffusivity,
    #[error("Invalid input file: control_planes.spacing must be above 0")]
    ControlPlaneSpacing,
    #[error(
        "Invalid input file: control_planes.spacing and length_to_well must be above 0, and well_thickness at least 0"
    )]
    WellControlSurface,
    /// A defect of the mesh
    #[error("Invalid mesh or flow solution: {0}")]
    InvalidNetwork(String),
    #[error(
        "Invalid mesh or flow solution: the in-flow and the out-flow boundary must each have at least one point (a point in both is in the in-flow boundary)"
    )]
    EmptyInflowOutflow,
    #[error("Invalid mesh or flow solution: the flow solution is not a solution of this mesh")]
    FlowNotOfMesh,
    #[error(
        "The mesh has no point data array boundary__{boundary}. Use the name of a boundary__<name> array of the mesh."
    )]
    UnknownBoundary { boundary: String },
    #[error("No particle was placed: the seeding gives no particle; check the in-flow boundary")]
    NoSeededParticle,
    #[error(
        "No particle was placed: all cells (but the last one) have a node on the out-flow boundary"
    )]
    NoCellForSeeding,
    #[error("No particle was placed: the in-flow volume flux is {volume_flux:e} m³/s, not above 0")]
    NoInflowFlux { volume_flux: f64 },
    #[error(transparent)]
    Flow(#[from] fehm::FehmError),
}
