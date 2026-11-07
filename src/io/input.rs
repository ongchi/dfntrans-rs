//! The input file (TOML) of a particle tracking run. All validation of the
//! input file is here; the names of the boundaries are checked against the
//! mesh in `network`.

use std::path::{Path, PathBuf};

use serde::Deserialize;

use crate::DfntransError;

/// Options of a particle tracking run
#[derive(Clone, Debug, PartialEq)]
pub struct Input {
    /// The fehm-rs input file, as written in the input file (relative to its folder)
    pub flow_file: PathBuf,
    /// Seed of the random numbers (`srand48`). The binary takes a seed from
    /// the clock for 0.
    pub seed: u32,
    pub time_unit: TimeUnit,
    /// Fracture porosity
    pub porosity: f64,
    /// Largest number of time steps of a particle
    pub max_time_steps: usize,
    /// Names of the `boundary__<name>` point arrays of the in-flow boundary
    pub inflow: Vec<String>,
    /// Names of the `boundary__<name>` point arrays of the out-flow boundary
    pub outflow: Vec<String>,
    pub routing: RoutingRule,
    pub seeding: Seeding,
    /// Matrix diffusion (time domain random walk)
    pub matrix_diffusion: Option<MatrixDiffusion>,
    pub control_surfaces: Option<ControlSurfaces>,
    pub output: OutputOptions,
}

/// The unit of the times and velocities of the output
#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum TimeUnit {
    Seconds,
    Minutes,
    Hours,
    Days,
    Years,
}

impl TimeUnit {
    /// Seconds in one unit (a year has 365 days)
    pub fn seconds(self) -> f64 {
        match self {
            TimeUnit::Seconds => 1.0,
            TimeUnit::Minutes => 60.0,
            TimeUnit::Hours => 60.0 * 60.0,
            TimeUnit::Days => 24.0 * 60.0 * 60.0,
            TimeUnit::Years => 365.0 * 24.0 * 60.0 * 60.0,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            TimeUnit::Seconds => "s",
            TimeUnit::Minutes => "min",
            TimeUnit::Hours => "h",
            TimeUnit::Days => "d",
            TimeUnit::Years => "a",
        }
    }
}

/// The rule which selects the out-going fracture at an intersection
#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum RoutingRule {
    CompleteMixing,
    StreamlineRouting,
}

/// The initial positions of the particles
#[derive(Clone, Debug, Deserialize, PartialEq)]
#[serde(tag = "method", rename_all = "snake_case", deny_unknown_fields)]
pub enum Seeding {
    /// The same number of particles on each in-flow line, at equal distances
    /// (DFNTrans `init_nf`)
    UniformPerLine {
        particles_per_line: usize,
        #[serde(default = "default_true")]
        weight_by_flux: bool,
    },
    /// Particles at one distance on all in-flow lines; the distance is the
    /// total length / (particles_per_line · number of lines) (`init_eqd`)
    Equidistant {
        particles_per_line: usize,
        #[serde(default = "default_true")]
        weight_by_flux: bool,
    },
    /// Particles at the centers of random triangles, weighted by aperture (`init_random`)
    Random { particle_count: usize },
    /// Particles on the in-flow lines in proportion to the in-flow flux of
    /// each point; all particles have the same weight (`init_fluxw`)
    FluxWeighted { particle_count: usize },
    /// Particles at each in-flow point, weighted by aperture (`init_well`)
    AtNodes { particles_per_node: usize },
    /// Particles at given points (x, y, z in the coordinates of the mesh),
    /// each in its nearest triangle, weighted by aperture (not in DFNTrans)
    Points {
        points: Vec<[f64; 3]>,
        particles_per_point: usize,
    },
}

fn default_true() -> bool {
    true
}

/// Matrix diffusion with the time domain random walk of DFNTrans (`tdrw`)
#[derive(Clone, Debug, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct MatrixDiffusion {
    /// Matrix porosity
    pub matrix_porosity: f64,
    /// Matrix diffusion coefficient (m²/s)
    pub matrix_diffusivity: f64,
}

/// Control planes or cylinders: the state of a particle is recorded where
/// it crosses one (DFNTrans `ControlPlane`, `ControlCylinder`)
#[derive(Clone, Debug, Deserialize, PartialEq)]
#[serde(tag = "shape", rename_all = "snake_case", deny_unknown_fields)]
pub enum ControlSurfaces {
    /// Planes normal to `direction`, from the in-flow boundary to the
    /// out-flow boundary, at a distance of `spacing`
    Plane { direction: Axis, spacing: f64 },
    /// Cylinders about the line through the origin along `axis`, from the
    /// radius `length_to_well` inward, at a distance of `spacing`, down to
    /// the radius `well_thickness` / 2
    Cylinder {
        axis: Axis,
        spacing: f64,
        length_to_well: f64,
        well_thickness: f64,
    },
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Axis {
    X,
    Y,
    Z,
}

impl Axis {
    pub fn index(self) -> usize {
        match self {
            Axis::X => 0,
            Axis::Y => 1,
            Axis::Z => 2,
        }
    }
}

/// Optional output files
#[derive(Clone, Debug, Default, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct OutputOptions {
    /// `trajectories.vtp`: the path of each particle
    #[serde(default)]
    pub trajectories: bool,
    /// Also the particles which do not reach the out-flow boundary
    #[serde(default)]
    pub all_particles: bool,
    /// `velocities.vtp`: the reconstructed velocity at each point
    #[serde(default)]
    pub velocities: bool,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct InputFile {
    flow: PathBuf,
    seed: u32,
    time_unit: TimeUnit,
    porosity: f64,
    max_time_steps: usize,
    inflow: Vec<String>,
    outflow: Vec<String>,
    routing: RoutingRule,
    seeding: Seeding,
    tdrw: Option<MatrixDiffusion>,
    control_planes: Option<ControlSurfaces>,
    #[serde(default)]
    output: OutputOptions,
}

/// Reads the input file. The flow file is relative to the folder of the input file.
pub fn read_input(path: &Path) -> Result<Input, DfntransError> {
    let text = std::fs::read_to_string(path).map_err(|source| DfntransError::ReadInput {
        file: path.display().to_string(),
        source,
    })?;
    let mut input = build_input_from_text(&text)?;
    if let Some(folder) = path.parent() {
        input.flow_file = folder.join(&input.flow_file);
    }
    Ok(input)
}

/// Builds the input from the text of the input file.
pub fn build_input_from_text(text: &str) -> Result<Input, DfntransError> {
    let file: InputFile =
        toml::from_str(text).map_err(|error| DfntransError::InvalidToml(error.to_string()))?;
    let input = Input {
        flow_file: file.flow,
        seed: file.seed,
        time_unit: file.time_unit,
        porosity: file.porosity,
        max_time_steps: file.max_time_steps,
        inflow: file.inflow,
        outflow: file.outflow,
        routing: file.routing,
        seeding: file.seeding,
        matrix_diffusion: file.tdrw,
        control_surfaces: file.control_planes,
        output: file.output,
    };
    validate(&input)?;
    Ok(input)
}

fn validate(input: &Input) -> Result<(), DfntransError> {
    if !(input.porosity > 0.0 && input.porosity <= 1.0) {
        return Err(DfntransError::Porosity);
    }
    if input.max_time_steps == 0 {
        return Err(DfntransError::MaxTimeSteps);
    }
    if input.inflow.is_empty() || input.outflow.is_empty() {
        return Err(DfntransError::InflowOutflowMissing);
    }
    if let Some(name) = input
        .inflow
        .iter()
        .find(|name| input.outflow.contains(name))
    {
        return Err(DfntransError::BoundaryInInflowAndOutflow {
            boundary: name.clone(),
        });
    }

    let count = match &input.seeding {
        Seeding::UniformPerLine {
            particles_per_line, ..
        }
        | Seeding::Equidistant {
            particles_per_line, ..
        } => particles_per_line,
        Seeding::Random { particle_count } | Seeding::FluxWeighted { particle_count } => {
            particle_count
        }
        Seeding::AtNodes { particles_per_node } => particles_per_node,
        Seeding::Points { points, .. } if points.is_empty() => {
            return Err(DfntransError::NoSeedingPoints);
        }
        Seeding::Points {
            particles_per_point,
            ..
        } => particles_per_point,
    };
    if *count == 0 {
        return Err(DfntransError::ParticleCount);
    }

    if let Some(diffusion) = &input.matrix_diffusion {
        if !(diffusion.matrix_porosity > 0.0 && diffusion.matrix_porosity <= 1.0) {
            return Err(DfntransError::MatrixPorosity);
        }
        if diffusion.matrix_diffusivity.is_nan() || diffusion.matrix_diffusivity <= 0.0 {
            return Err(DfntransError::MatrixDiffusivity);
        }
    }

    match input.control_surfaces {
        Some(ControlSurfaces::Plane { spacing, .. }) if spacing.is_nan() || spacing <= 0.0 => {
            Err(DfntransError::ControlPlaneSpacing)
        }
        Some(ControlSurfaces::Cylinder {
            spacing,
            length_to_well,
            well_thickness,
            ..
        }) if !(spacing > 0.0 && length_to_well > 0.0 && well_thickness >= 0.0) => {
            Err(DfntransError::WellControlSurface)
        }
        _ => Ok(()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const VALID: &str = r#"
flow = "flow.toml"
seed = 1
time_unit = "days"
porosity = 1.0
max_time_steps = 1000
inflow = ["west"]
outflow = ["east"]
routing = "complete_mixing"

[seeding]
method = "flux_weighted"
particle_count = 100
"#;

    #[test]
    fn reads_a_valid_file() {
        let input = build_input_from_text(VALID).unwrap();
        assert_eq!(input.time_unit, TimeUnit::Days);
        assert_eq!(
            input.seeding,
            Seeding::FluxWeighted {
                particle_count: 100
            }
        );
        assert_eq!(input.output, OutputOptions::default());
        assert!(input.matrix_diffusion.is_none());
    }

    #[test]
    fn reads_the_optional_sections() {
        let text = format!(
            "{VALID}\n[tdrw]\nmatrix_porosity = 0.02\nmatrix_diffusivity = 1e-13\n\
             [control_planes]\nshape = \"plane\"\ndirection = \"x\"\nspacing = 1.0\n\
             [output]\ntrajectories = true\n"
        );
        let input = build_input_from_text(&text).unwrap();
        assert_eq!(
            input.control_surfaces,
            Some(ControlSurfaces::Plane {
                direction: Axis::X,
                spacing: 1.0
            })
        );
        assert!(input.output.trajectories);
        assert!(input.matrix_diffusion.is_some());
    }

    fn input_error(text: &str) -> DfntransError {
        build_input_from_text(text).unwrap_err()
    }

    #[test]
    fn rejects_an_unknown_key() {
        let text = VALID.replace("seed = 1", "seed = 1\ndensity = 1000.0");
        let error = input_error(&text);
        assert!(matches!(error, DfntransError::InvalidToml(_)), "{error}");
    }

    #[test]
    fn rejects_an_unknown_key_of_the_seeding() {
        let text = VALID.replace("particle_count = 100", "particle_count = 100\nextra = 1");
        let error = input_error(&text);
        assert!(matches!(error, DfntransError::InvalidToml(_)), "{error}");
    }

    #[test]
    fn rejects_a_missing_key() {
        let text = VALID.replace("porosity = 1.0", "");
        let error = input_error(&text);
        assert!(matches!(error, DfntransError::InvalidToml(_)), "{error}");
    }

    #[test]
    fn rejects_a_boundary_in_inflow_and_outflow() {
        let text = VALID.replace(r#"["east"]"#, r#"["west"]"#);
        let error = input_error(&text);
        assert!(
            matches!(&error, DfntransError::BoundaryInInflowAndOutflow { boundary } if boundary == "west"),
            "{error}"
        );
    }

    #[test]
    fn rejects_zero_particles() {
        let text = VALID.replace("particle_count = 100", "particle_count = 0");
        let error = input_error(&text);
        assert!(matches!(error, DfntransError::ParticleCount), "{error}");
    }

    fn with_points(points: &str) -> String {
        VALID.replace(
            "method = \"flux_weighted\"\nparticle_count = 100",
            &format!("method = \"points\"\npoints = {points}\nparticles_per_point = 3"),
        )
    }

    #[test]
    fn reads_the_seeding_points() {
        let input = build_input_from_text(&with_points("[[1.0, 2.0, 3.0]]")).unwrap();
        assert_eq!(
            input.seeding,
            Seeding::Points {
                points: vec![[1.0, 2.0, 3.0]],
                particles_per_point: 3
            }
        );
    }

    #[test]
    fn rejects_no_seeding_point() {
        let error = input_error(&with_points("[]"));
        assert!(matches!(error, DfntransError::NoSeedingPoints), "{error}");
    }

    #[test]
    fn rejects_a_seeding_point_of_two_coordinates() {
        let error = input_error(&with_points("[[1.0, 2.0]]"));
        assert!(matches!(error, DfntransError::InvalidToml(_)), "{error}");
    }

    #[test]
    fn rejects_a_porosity_above_one() {
        let error = input_error(&VALID.replace("porosity = 1.0", "porosity = 1.5"));
        assert!(matches!(error, DfntransError::Porosity), "{error}");
    }

    #[test]
    fn rejects_zero_time_steps() {
        let error = input_error(&VALID.replace("max_time_steps = 1000", "max_time_steps = 0"));
        assert!(matches!(error, DfntransError::MaxTimeSteps), "{error}");
    }

    #[test]
    fn rejects_no_outflow_boundary() {
        let error = input_error(&VALID.replace(r#"["east"]"#, "[]"));
        assert!(
            matches!(error, DfntransError::InflowOutflowMissing),
            "{error}"
        );
    }

    #[test]
    fn rejects_a_zero_matrix_porosity() {
        let text = format!("{VALID}\n[tdrw]\nmatrix_porosity = 0.0\nmatrix_diffusivity = 1e-13\n");
        let error = input_error(&text);
        assert!(matches!(error, DfntransError::MatrixPorosity), "{error}");
    }

    #[test]
    fn rejects_a_zero_matrix_diffusivity() {
        let text = format!("{VALID}\n[tdrw]\nmatrix_porosity = 0.02\nmatrix_diffusivity = 0.0\n");
        let error = input_error(&text);
        assert!(matches!(error, DfntransError::MatrixDiffusivity), "{error}");
    }

    #[test]
    fn rejects_a_zero_plane_spacing() {
        let text = format!(
            "{VALID}\n[control_planes]\nshape = \"plane\"\ndirection = \"x\"\nspacing = 0.0\n"
        );
        let error = input_error(&text);
        assert!(
            matches!(error, DfntransError::ControlPlaneSpacing),
            "{error}"
        );
    }

    #[test]
    fn rejects_a_negative_well_thickness() {
        let text = format!(
            "{VALID}\n[control_planes]\nshape = \"cylinder\"\naxis = \"z\"\nspacing = 1.0\n\
             length_to_well = 10.0\nwell_thickness = -1.0\n"
        );
        let error = input_error(&text);
        assert!(
            matches!(error, DfntransError::WellControlSurface),
            "{error}"
        );
    }
}
