//! Particle tracking on the reconstructed velocities (DFNTrans
//! `TrackingPart.c`): predictor-corrector steps in the 2D coordinates of a
//! fracture, the search of the cell of a particle, the move through
//! fracture intersections, and the exit at the out-flow boundary.
//!
//! The port keeps the order of the operations of DFNTrans, thus the same
//! seed gives the same random numbers at the same decisions.

use fehm::Solution;
use tracing::{info, warn};

use crate::DfntransError;
use crate::io::input::{ControlSurfaces, Input, RoutingRule};
use crate::network::{FlowZone, Network, build_network, norm, subtract};
use crate::random::Drand48;
use crate::routing::{IntersectionCell, complete_mixing, streamline_routing};
use crate::seeding::seed_particles;
use crate::tdrw::{MatrixDiffusion, diffusion_time};
use crate::velocity::reconstruct_velocities;

/// The result of a run
#[derive(Clone, Debug)]
pub struct Transport {
    /// The network with the reconstructed velocities
    pub network: Network,
    /// The particles of the output: those which reached the out-flow
    /// boundary, or all with `output.all_particles`
    pub particles: Vec<ParticleResult>,
    /// Number of particles which were placed
    pub particle_count: usize,
    /// Number of particles which reached the out-flow boundary
    pub exit_count: usize,
    /// Volume flux into the network through the in-flow boundary (m³/s)
    pub inflow_volume_flux: f64,
    /// Volume flux out of the network through the out-flow boundary (m³/s)
    pub outflow_volume_flux: f64,
    pub has_matrix_diffusion: bool,
}

/// The path and the times of one particle
#[derive(Clone, Debug)]
pub struct ParticleResult {
    /// Index of the particle in the order of the seeding
    pub index: usize,
    pub has_exited: bool,
    /// The last time step
    pub time_steps: usize,
    pub flux_weight: f64,
    /// Advective travel time (time unit)
    pub travel_time: f64,
    /// With matrix diffusion
    pub diffusion: Option<DiffusionTimes>,
    pub final_position: [f64; 3],
    /// Lagrangian β: Σ 2·length / (speed · aperture) (time unit / m)
    pub beta: f64,
    /// Length of the path (m)
    pub length: f64,
    /// Points of the path; empty without `output.trajectories`
    pub trajectory: Vec<PathPoint>,
    /// The crossings of the control planes or cylinders
    pub control_crossings: Vec<ControlCrossing>,
}

/// Times with matrix diffusion
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiffusionTimes {
    /// Advective time + diffusion time
    pub total_time: f64,
    pub diffusion_time: f64,
}

/// The state of a particle at one point of its path
#[derive(Clone, Copy, Debug)]
pub struct PathPoint {
    pub position: [f64; 3],
    pub velocity: [f64; 3],
    /// Index of the fracture in `Network.fractures`
    pub fracture: usize,
    /// Advective time
    pub time: f64,
    pub aperture: f64,
    /// MPa
    pub pressure: f64,
}

/// The state of a particle where it crosses a control plane or cylinder
#[derive(Clone, Copy, Debug)]
pub struct ControlCrossing {
    pub point: PathPoint,
    pub length: f64,
    pub diffusion: Option<DiffusionTimes>,
}

/// Where a particle is in its cell (DFNTrans `intcell`)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum CellLocation {
    /// 0
    Interior,
    /// 1: a node of the cell is on an intersection
    Intersection,
    /// 2: two nodes of the cell are on the fracture edge
    Boundary,
    /// 3: on an intersection and on the fracture edge
    IntersectionBoundary,
    /// 4: the particle just moved through an intersection
    CrossedIntersection,
}

/// DFNTrans `FLAG_OUT`
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ExitState {
    /// 0
    Inside,
    /// 1: through the out-flow boundary
    Outflow,
    /// 3: out of the fracture through its edge
    OutsideFracture,
}

/// One particle (DFNTrans `struct contam`)
#[derive(Clone, Debug)]
pub(crate) struct Particle {
    pub velocity: [f64; 2],
    pub position: [f64; 2],
    pub previous_position: [f64; 2],
    pub fracture: usize,
    pub cell: Option<usize>,
    pub location: CellLocation,
    /// Interpolation weights of the corners of the cell
    pub weights: [f64; 3],
    /// Advective time
    pub time: f64,
    pub flux_weight: f64,
    pub pressure: f64,
}

impl Particle {
    pub fn new(position: [f64; 2], fracture: usize, flux_weight: f64) -> Self {
        Self {
            velocity: [0.0; 2],
            position,
            previous_position: [0.0; 2],
            fracture,
            cell: None,
            location: CellLocation::Interior,
            weights: [0.0; 3],
            time: 0.0,
            flux_weight,
            pressure: 0.0,
        }
    }
}

/// The state of the run of one particle
struct Run {
    particle: Particle,
    exit: ExitState,
    /// Advective time at the last matrix diffusion step (`t_adv0`)
    diffused_advective_time: f64,
    diffusion: DiffusionTimes,
    trajectory: Vec<PathPoint>,
}

/// The control surfaces in the coordinates of the network
#[derive(Clone, Copy, Debug)]
enum ControlGeometry {
    /// Planes `axis` = start + k · step
    Plane { axis: usize, start: f64, step: f64 },
    /// Cylinders of radius start − k · spacing about the line through the
    /// origin along `axis`, down to `smallest_radius`
    Cylinder {
        axis: usize,
        start: f64,
        spacing: f64,
        smallest_radius: f64,
    },
}

/// Runs the particle tracking on the flow solution of fehm-rs.
pub fn track_particles(
    input: &Input,
    flow: &fehm::Input,
    solution: &Solution,
) -> Result<Transport, DfntransError> {
    let mut network = build_network(flow, solution, &input.inflow, &input.outflow)?;
    info!(
        "Network: {} nodes, {} triangles, {} fractures; {} in-flow and {} out-flow nodes",
        network.nodes.len(),
        network.triangles.len(),
        network.fractures.len(),
        network.inflow_nodes.len(),
        network.outflow_nodes.len()
    );

    reconstruct_velocities(&mut network, input.porosity, input.time_unit.seconds());

    let volume_flux = |nodes: &[usize]| -> f64 {
        nodes
            .iter()
            .map(|&node| network.nodes[node].total_mass_flux() / network.nodes[node].density)
            .sum()
    };
    let inflow_volume_flux = volume_flux(&network.inflow_nodes);
    let outflow_volume_flux = -volume_flux(&network.outflow_nodes);
    info!(
        "In-flow volume flux {inflow_volume_flux:.5e} m³/s, out-flow volume flux {outflow_volume_flux:.5e} m³/s"
    );

    let mut tracker = Tracker::new(&network, input);
    let particles = seed_particles(&mut tracker, &input.seeding, inflow_volume_flux)?;
    if particles.is_empty() {
        return Err(DfntransError::NoParticles(
            "the seeding gives no particle; check the in-flow boundary".to_string(),
        ));
    }
    let particle_count = particles.len();
    info!(
        "Seeding {particle_count} particles, at most {} time steps each",
        input.max_time_steps
    );

    let mut results = Vec::new();
    let mut exit_count = 0;
    let mut next_report = 10;
    for (index, particle) in particles.into_iter().enumerate() {
        let percent = 100 * index / particle_count;
        if percent >= next_report {
            info!("{index} of {particle_count} particles done, {exit_count} exited");
            next_report += 10;
        }
        if let Some(result) = tracker.track(index, particle) {
            if result.has_exited {
                exit_count += 1;
            }
            if result.has_exited || input.output.all_particles {
                results.push(result);
            }
        }
    }
    info!(
        "{exit_count} of {particle_count} particles exited through the out-flow boundary ({:.2} %)",
        100.0 * exit_count as f64 / particle_count as f64
    );

    Ok(Transport {
        network,
        particles: results,
        particle_count,
        exit_count,
        inflow_volume_flux,
        outflow_volume_flux,
        has_matrix_diffusion: input.matrix_diffusion.is_some(),
    })
}

/// The tracking of the particles on one network
pub(crate) struct Tracker<'a> {
    pub network: &'a Network,
    input: &'a Input,
    pub random: Drand48,
    matrix_diffusion: Option<MatrixDiffusion>,
    control: Option<ControlGeometry>,
}

impl<'a> Tracker<'a> {
    fn new(network: &'a Network, input: &'a Input) -> Self {
        let matrix_diffusion = input
            .matrix_diffusion
            .as_ref()
            .map(|diffusion| MatrixDiffusion {
                porosity: diffusion.matrix_porosity,
                diffusivity: diffusion.matrix_diffusivity * input.time_unit.seconds(),
            });
        let control = input.control_surfaces.as_ref().and_then(|surfaces| match *surfaces {
            ControlSurfaces::Plane { direction, spacing } => {
                // From the first in-flow node toward the first out-flow node
                let axis = direction.index();
                let start = network.nodes[network.inflow_nodes[0]].position[axis];
                let end = network.nodes[network.outflow_nodes[0]].position[axis];
                if end == start {
                    warn!("The first in-flow and out-flow nodes have the same coordinate; no control planes");
                    return None;
                }
                Some(ControlGeometry::Plane {
                    axis,
                    start,
                    step: spacing * (end - start).signum(),
                })
            }
            ControlSurfaces::Cylinder {
                axis,
                spacing,
                length_to_well,
                well_thickness,
            } => Some(ControlGeometry::Cylinder {
                axis: axis.index(),
                start: length_to_well,
                spacing,
                smallest_radius: well_thickness / 2.0,
            }),
        });

        Self {
            network,
            input,
            random: Drand48::new(input.seed),
            matrix_diffusion,
            control,
        }
    }

    /// Tracks one particle (the loop on the time steps of `ParticleTrack`).
    /// `None` if its first cell is not found.
    fn track(&mut self, index: usize, particle: Particle) -> Option<ParticleResult> {
        let mut run = Run {
            particle,
            exit: ExitState::Inside,
            diffused_advective_time: 0.0,
            diffusion: DiffusionTimes {
                total_time: 0.0,
                diffusion_time: 0.0,
            },
            trajectory: Vec::new(),
        };
        if run.particle.cell.is_none() && !self.find_initial_cell(&mut run.particle) {
            warn!(
                "Particle {}: no cell found at the initial position {:?} in fracture {}",
                index + 1,
                run.particle.position,
                self.network.fractures[run.particle.fracture].id
            );
            return None;
        }

        let records_trajectory = self.input.output.trajectories;
        let max_steps = self.input.max_time_steps;
        let mut previous_cell = run.particle.cell;
        let mut previous_fracture = run.particle.fracture;
        let start = self.position_3d(&run.particle);
        let mut step_start = start;
        let mut lagrangian_start = start;
        let mut beta = 0.0;
        let mut length = 0.0;
        let mut steps_in_cell = 0;
        let mut fracture_changes = 0;
        let mut last_step = 0;

        let mut crossings = Vec::new();
        let mut next_control = match self.control {
            Some(ControlGeometry::Plane { start, step, .. }) => {
                crossings.push(self.control_crossing(&run, 0.0));
                start + step
            }
            Some(ControlGeometry::Cylinder { start, .. }) => start,
            None => 0.0,
        };

        for step in 0..max_steps {
            let position = self.position_3d(&run.particle);
            length += norm(subtract(position, step_start));
            step_start = position;
            if records_trajectory {
                run.trajectory.push(self.path_point(&run.particle));
            }

            if let Some(geometry) = self.control
                && let Some(crossing) = self.cross_control_surface(
                    &mut run,
                    geometry,
                    &mut next_control,
                    position,
                    length,
                )
            {
                crossings.push(crossing);
            }

            if matches!(
                run.particle.location,
                CellLocation::Intersection | CellLocation::IntersectionBoundary
            ) {
                self.check_distance(&mut run, records_trajectory);
                previous_cell = run.particle.cell;
            }
            if run.particle.cell.is_none() {
                break;
            }

            self.predictor_step(&mut run.particle);
            self.check_new_cell(&mut run);
            if run.exit == ExitState::Outflow {
                break;
            }

            if run.particle.cell.is_some() {
                if !matches!(
                    run.particle.location,
                    CellLocation::Intersection | CellLocation::IntersectionBoundary
                ) {
                    self.corrector_step(&mut run.particle);
                }
                let (tau, step_beta) = self.lagrangian(&run.particle, step_start, lagrangian_start);
                lagrangian_start = step_start;
                beta += step_beta;
                run.particle.time += tau;
            } else if run.exit == ExitState::OutsideFracture {
                self.move_from_fracture_edge(&mut run.particle, previous_cell);
                run.exit = ExitState::Inside;
            }

            if run.particle.cell.is_some() {
                self.check_new_cell(&mut run);
            }
            if run.particle.cell.is_none() && run.exit == ExitState::OutsideFracture {
                run.exit = ExitState::Inside;
                self.move_from_fracture_edge(&mut run.particle, previous_cell);
            }

            if run.particle.cell.is_some() && run.particle.cell != previous_cell {
                steps_in_cell = 0;
                previous_cell = run.particle.cell;
                if previous_fracture != run.particle.fracture {
                    previous_fracture = run.particle.fracture;
                    fracture_changes += 1;
                    // A particle which visits more fractures than the network has goes in a circle
                    if fracture_changes > self.network.fractures.len() {
                        run.exit = ExitState::Inside;
                        break;
                    }
                }
            } else if run.particle.cell == previous_cell {
                // A particle which stays in one cell for 1/10 of the time steps is stuck
                steps_in_cell += 1;
                if steps_in_cell as f64 > max_steps as f64 / 10.0 {
                    run.exit = ExitState::Inside;
                    break;
                }
            }

            last_step = step;
            if run.particle.cell.is_none() {
                run.exit = ExitState::Inside;
                break;
            }
        }

        let has_exited = run.exit == ExitState::Outflow;
        if !(has_exited || self.input.output.all_particles) {
            return Some(ParticleResult {
                index,
                has_exited,
                time_steps: last_step,
                flux_weight: run.particle.flux_weight,
                travel_time: run.particle.time,
                diffusion: None,
                final_position: self.position_3d(&run.particle),
                beta,
                length,
                trajectory: Vec::new(),
                control_crossings: Vec::new(),
            });
        }

        if run.particle.cell.is_some() {
            self.move_to_outflow_boundary(&mut run.particle);
        }
        let position = self.position_3d(&run.particle);
        if records_trajectory {
            run.trajectory.push(self.path_point(&run.particle));
        }
        self.diffuse(&mut run);
        if run.particle.cell.is_some() {
            let (_, step_beta) = self.lagrangian(&run.particle, step_start, lagrangian_start);
            beta += step_beta;
        }
        length += norm(subtract(position, step_start));

        if self.control.is_some() {
            if run.particle.cell.is_none() {
                run.particle.cell = previous_cell;
            }
            self.diffuse(&mut run);
            crossings.push(self.control_crossing(&run, length));
        }

        Some(ParticleResult {
            index,
            has_exited,
            time_steps: last_step,
            flux_weight: run.particle.flux_weight,
            travel_time: run.particle.time,
            diffusion: self.matrix_diffusion.map(|_| run.diffusion),
            final_position: position,
            beta,
            length,
            trajectory: run.trajectory,
            control_crossings: crossings,
        })
    }

    /// The crossing of the next control surface in this time step, if any
    fn cross_control_surface(
        &mut self,
        run: &mut Run,
        geometry: ControlGeometry,
        next: &mut f64,
        position: [f64; 3],
        length: f64,
    ) -> Option<ControlCrossing> {
        match geometry {
            ControlGeometry::Plane { axis, step, .. } => {
                let has_crossed = if step > 0.0 {
                    position[axis] >= *next
                } else {
                    position[axis] <= *next
                };
                if !has_crossed {
                    return None;
                }
                self.diffuse(run);
                *next += step;
            }
            ControlGeometry::Cylinder {
                axis,
                spacing,
                smallest_radius,
                ..
            } => {
                let distance = match axis {
                    0 => position[1].hypot(position[2]),
                    1 => position[0].hypot(position[2]),
                    _ => position[0].hypot(position[1]),
                };
                if !(*next >= smallest_radius && distance <= *next) {
                    return None;
                }
                *next -= spacing;
            }
        }
        Some(self.control_crossing(run, length))
    }

    fn control_crossing(&self, run: &Run, length: f64) -> ControlCrossing {
        ControlCrossing {
            point: self.path_point(&run.particle),
            length,
            diffusion: self.matrix_diffusion.map(|_| run.diffusion),
        }
    }

    fn path_point(&self, particle: &Particle) -> PathPoint {
        PathPoint {
            position: self.position_3d(particle),
            velocity: self.velocity_3d(particle),
            fracture: particle.fracture,
            time: particle.time,
            aperture: particle.cell.map_or(0.0, |cell| {
                self.network.nodes[self.network.triangles[cell].nodes[0]].aperture
            }),
            pressure: particle.pressure,
        }
    }

    /// One matrix diffusion step for the advective time since the last one
    /// (the uses of `TimeDomainRW`)
    fn diffuse(&mut self, run: &mut Run) {
        let Some(parameters) = self.matrix_diffusion else {
            return;
        };
        let advective_time = run.particle.time - run.diffused_advective_time;
        let aperture = self.diffusion_aperture(&run.particle);
        let time = diffusion_time(parameters, aperture, advective_time, &mut self.random);
        run.diffused_advective_time = run.particle.time;
        run.diffusion.diffusion_time += time;
        run.diffusion.total_time += advective_time + time;
    }

    /// The aperture of the first cell node which is not on an intersection,
    /// or of the first node of the fracture without a cell
    fn diffusion_aperture(&self, particle: &Particle) -> f64 {
        let nodes = &self.network.nodes;
        match particle.cell {
            Some(cell) => {
                let [n1, n2, n3] = self.network.triangles[cell].nodes;
                let off_line = |node: usize| !matches!(nodes[node].type_number(), 2 | 12);
                let node = if off_line(n1) {
                    n1
                } else if off_line(n2) {
                    n2
                } else {
                    n3
                };
                nodes[node].aperture
            }
            None => nodes[self.network.fractures[particle.fracture].first_node].aperture,
        }
    }

    pub fn position_3d(&self, particle: &Particle) -> [f64; 3] {
        self.network.fractures[particle.fracture].to_space(particle.position)
    }

    fn velocity_3d(&self, particle: &Particle) -> [f64; 3] {
        self.network.fractures[particle.fracture].vector_to_space(particle.velocity)
    }

    /// The 2D position of a node in the fracture of the particle
    fn node_position(&self, node: usize, particle: &Particle) -> [f64; 2] {
        self.network.nodes[node].position_in(particle.fracture)
    }

    /// Barycentric weights of the particle position in a triangle; weights
    /// near 0 are set to 0 (`CalculateWeights`)
    fn weights(&self, particle: &Particle, [n1, n2, n3]: [usize; 3]) -> [f64; 3] {
        const TOLERANCE: f64 = 10e-5;
        let [x1, y1] = self.node_position(n1, particle);
        let [x2, y2] = self.node_position(n2, particle);
        let [x3, y3] = self.node_position(n3, particle);
        let [x, y] = particle.position;
        let determinant = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);
        let mut first = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / determinant;
        let mut second = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / determinant;
        let snap = |weight: &mut f64| {
            if *weight >= -TOLERANCE && *weight <= TOLERANCE {
                *weight = 0.0;
            }
        };
        snap(&mut first);
        snap(&mut second);
        let mut third = 1.0 - first - second;
        snap(&mut third);
        [first, second, third]
    }

    /// True if the particle is in the cell (within 1e-5 of the weights);
    /// then the cell, the weights and the location of the particle are set
    /// (`InsideCell`)
    pub fn is_inside(&self, particle: &mut Particle, cell: usize) -> bool {
        const TOLERANCE: f64 = 1e-5;
        let nodes = self.network.triangles[cell].nodes;
        let weights = self.weights(particle, nodes);
        if !weights
            .iter()
            .all(|weight| (-TOLERANCE..=1.0 + TOLERANCE).contains(weight))
        {
            return false;
        }

        particle.weights = weights.map(f64::abs);
        particle.cell = Some(cell);
        let has_crossed = particle.location == CellLocation::CrossedIntersection;
        let type_of = |node: usize| self.network.nodes[node].type_number();
        let edge_count = nodes.iter().filter(|&&node| type_of(node) == 10).count();
        particle.location = if edge_count > 1 {
            CellLocation::Boundary
        } else {
            CellLocation::Interior
        };
        if nodes.iter().any(|&node| matches!(type_of(node), 2 | 12)) && !has_crossed {
            particle.location = CellLocation::Intersection;
        }
        if edge_count > 0 && particle.location == CellLocation::Intersection {
            particle.location = CellLocation::IntersectionBoundary;
        }
        true
    }

    /// The first cell around the in-flow nodes which has the particle
    /// (`InitCell`). Unlike DFNTrans, only cells of the fracture of the particle.
    pub fn find_initial_cell(&self, particle: &mut Particle) -> bool {
        for &node in &self.network.inflow_nodes {
            for neighbor in &self.network.nodes[node].neighbors {
                for &cell in &neighbor.triangles {
                    if self.network.triangles[cell].fracture == particle.fracture
                        && self.is_inside(particle, cell)
                    {
                        return true;
                    }
                }
            }
        }
        false
    }

    /// The cells around a node in the fracture of the particle (`NeighborCells`)
    fn search_cells_of_node(&self, particle: &mut Particle, node: usize) {
        for neighbor in &self.network.nodes[node].neighbors {
            // As DFNTrans: the result of the last cell of the edge counts
            let mut is_inside = false;
            for &cell in &neighbor.triangles {
                if self.network.triangles[cell].fracture == particle.fracture {
                    is_inside = self.is_inside(particle, cell);
                }
            }
            if is_inside {
                break;
            }
        }
    }

    /// The cells around the node with the largest weight, then around the
    /// other nodes (`SearchNeighborCells`)
    fn search_neighbor_cells(&self, particle: &mut Particle, [n1, n2, n3]: [usize; 3]) {
        let [w1, w2, w3] = particle.weights;
        let mut first = None;
        if w1 >= w2 && w1 >= w3 {
            first = Some(n1);
        }
        if w2 >= w1 && w2 >= w3 {
            first = Some(n2);
        }
        if w3 >= w2 && w3 >= w1 {
            first = Some(n3);
        }
        if let Some(node) = first {
            self.search_cells_of_node(particle, node);
        }
        for node in [n1, n2, n3] {
            if particle.cell.is_none() && first != Some(node) {
                self.search_cells_of_node(particle, node);
            }
        }
    }

    /// Time step at the particle: interpolation of the corner time steps (`CalculateCurrentDT`)
    fn time_step(&self, particle: &Particle) -> f64 {
        let cell = particle.cell.expect("the particle has a cell");
        (0..3)
            .map(|corner| self.network.corner_time_step(cell, corner) * particle.weights[corner])
            .sum()
    }

    /// Velocity and pressure from the corners, and an Euler step (`PredictorStep`)
    fn predictor_step(&self, particle: &mut Particle) {
        let cell = particle.cell.expect("the particle has a cell");
        let time_step = self.time_step(particle);
        let mut velocity = [0.0; 2];
        let mut pressure = 0.0;
        for corner in 0..3 {
            let corner_velocity = self.network.corner_velocity(cell, corner);
            let weight = particle.weights[corner];
            velocity[0] += weight * corner_velocity[0];
            velocity[1] += weight * corner_velocity[1];
            pressure +=
                weight * self.network.nodes[self.network.triangles[cell].nodes[corner]].pressure;
        }
        particle.velocity = velocity;
        particle.pressure = pressure;
        particle.previous_position = particle.position;
        particle.position = [
            particle.position[0] + time_step * velocity[0],
            particle.position[1] + time_step * velocity[1],
        ];
    }

    /// The step again from the previous position with the time step of the
    /// new position (`CorrectorStep`)
    fn corrector_step(&self, particle: &mut Particle) {
        let time_step = self.time_step(particle);
        particle.position = [
            particle.previous_position[0] + time_step * particle.velocity[0],
            particle.previous_position[1] + time_step * particle.velocity[1],
        ];
    }

    /// After a step: the new cell of the particle, or its exit (`CheckNewCell`)
    fn check_new_cell(&self, run: &mut Run) {
        let particle = &mut run.particle;
        let cell = particle.cell.expect("the particle has a cell");
        let cell_nodes = self.network.triangles[cell].nodes;
        let time_step = self.time_step(particle);
        let weights = self.weights(particle, cell_nodes);
        if weights.iter().all(|&weight| (0.0..=1.0).contains(&weight)) {
            particle.weights = weights;
            return;
        }

        particle.cell = None;
        self.search_neighbor_cells(particle, cell_nodes);
        let nodes = &self.network.nodes;
        let type_of = |node: usize| nodes[node].type_number();
        if cell_nodes
            .iter()
            .any(|&node| nodes[node].zone == Some(FlowZone::Outflow))
        {
            run.exit = ExitState::Outflow;
        }
        if particle.cell.is_some() {
            return;
        }

        let is_on_edge = cell_nodes
            .iter()
            .any(|&node| matches!(type_of(node), 10 | 12));
        let is_on_outflow_edge = cell_nodes
            .iter()
            .any(|&node| matches!(type_of(node), 210 | 212));
        if is_on_edge && !is_on_outflow_edge {
            // Out of the fracture: reflect the step in x, else in y
            for reflection in [[-1.0, 1.0], [1.0, -1.0]] {
                particle.position = [
                    particle.previous_position[0]
                        + time_step * particle.velocity[0] * reflection[0],
                    particle.previous_position[1]
                        + time_step * particle.velocity[1] * reflection[1],
                ];
                self.search_neighbor_cells(particle, cell_nodes);
                if particle.cell.is_some() {
                    break;
                }
            }
            if particle.cell.is_none() {
                run.exit = ExitState::OutsideFracture;
            }
        }

        let is_out_node = |node: usize| matches!(type_of(node), 210 | 212 | 200);
        let is_next_to_out_edge = |node: usize| {
            nodes[node]
                .neighbors
                .iter()
                .any(|neighbor| matches!(type_of(neighbor.node), 210 | 212))
        };
        if cell_nodes.iter().any(|&node| is_out_node(node))
            || cell_nodes.iter().any(|&node| is_next_to_out_edge(node))
        {
            run.exit = ExitState::Outflow;
        }
    }

    /// The fracture of a cell for the particle: its position in the
    /// coordinates of that fracture (`ChangeFracture`)
    fn change_fracture(&self, particle: &mut Particle, cell: usize) {
        let position = self.position_3d(particle);
        let fracture = self.network.triangles[cell].fracture;
        particle.position = self.network.fractures[fracture].to_plane(position);
        particle.fracture = fracture;
        particle.cell = Some(cell);
    }

    /// In a cell with a node on an intersection: if the next step crosses
    /// the intersection, the particle moves to the crossing point and into
    /// a cell of the intersection (`CheckDistance`)
    fn check_distance(&mut self, run: &mut Run, records_trajectory: bool) {
        let Some(cell) = run.particle.cell else {
            return;
        };
        let nodes = &self.network.nodes;
        let cell_nodes = self.network.triangles[cell].nodes;
        let time_step = self.time_step(&run.particle);
        let [px, py] = run.particle.position;
        let step_length =
            (run.particle.velocity[0] * time_step).hypot(run.particle.velocity[1] * time_step);
        let is_line_node = |node: usize| matches!(nodes[node].type_number(), 2 | 12);

        // The intersection edge: two nodes of the cell, or a node of the
        // cell and its neighbor on the intersection
        #[derive(PartialEq)]
        enum Third {
            Missing,
            InCell,
            Node,
        }
        let mut first = None;
        let mut second = None;
        for node in cell_nodes {
            if is_line_node(node) {
                if first.is_some() {
                    second = Some(node);
                } else {
                    first = Some(node);
                }
            }
        }
        let mut third = if first.is_some() && second.is_some() {
            Third::InCell
        } else {
            Third::Missing
        };
        let Some(int1) = first else {
            return;
        };
        if second.is_none() {
            for neighbor in &nodes[int1].neighbors {
                let other = &nodes[neighbor.node];
                if matches!(other.base_type_number(), 2 | 12)
                    && other.is_in_fracture(run.particle.fracture)
                {
                    if second.is_none() {
                        second = Some(neighbor.node);
                    } else {
                        third = Third::Node;
                        break;
                    }
                }
            }
        }
        let Some(int2) = second else {
            return;
        };
        let [cx1, cy1] = self.node_position(int1, &run.particle);
        let [cx2, cy2] = self.node_position(int2, &run.particle);

        if third != Third::Missing {
            // Height of the triangle of the particle on the intersection segment
            let base = (cx2 - cx1).hypot(cy2 - cy1);
            let side2 = (cx2 - px).hypot(cy2 - py);
            let side1 = (cx1 - px).hypot(cy1 - py);
            let half_perimeter = (side1 + side2 + base) / 2.0;
            let area = (half_perimeter
                * (half_perimeter - base)
                * (half_perimeter - side1)
                * (half_perimeter - side2))
                .sqrt();
            let height = area / (base * 0.5);
            if height <= step_length {
                let [px2, py2] = run.particle.previous_position;
                self.predictor_step(&mut run.particle);
                let [px1, py1] = run.particle.position;
                let pr1 = (px1 - cx1) * (py2 - cy1) - (py1 - cy1) * (px2 - cx1);
                let pr2 = (cx1 - px1) * (cy2 - py1) - (cy1 - py1) * (cx2 - px1);
                let pr3 = (px1 - cx2) * (py2 - cy2) - (py1 - cy2) * (px2 - cx2);
                let pr4 = (cx1 - px2) * (cy2 - py2) - (cy1 - py2) * (cx2 - px2);
                if pr1 * pr3 < 0.0 && pr2 * pr4 < 0.0 {
                    let q1 = cx1 * cy2 - cy1 * cx2;
                    let q2 = px1 * py2 - py1 * px2;
                    let q3 = (cx1 - cx2) * (py1 - py2) - (cy1 - cy2) * (px1 - px2);
                    run.particle.position = [
                        ((px1 - px2) * q1 - (cx1 - cx2) * q2) / q3,
                        ((py1 - py2) * q1 - (cy1 - cy2) * q2) / q3,
                    ];
                    self.check_new_cell(run);
                    if run.particle.cell.is_some() {
                        self.cross_intersection(run, cell, [int1, int2], records_trajectory);
                    }
                }
            }
        }

        if third == Third::Missing {
            // The cell is at the end of the intersection
            let middle = [(cx1 + cx2) / 2.0, (cy1 + cy2) / 2.0];
            let squared_distance = |position: [f64; 2]| {
                (position[0] - middle[0]).powi(2) + (position[1] - middle[1]).powi(2)
            };
            let initial_distance = squared_distance(run.particle.position);
            self.predictor_step(&mut run.particle);
            let moves_to_line = initial_distance >= squared_distance(run.particle.position);
            let moves_out = !moves_to_line && {
                // As DFNTrans: the cells are counted only if int2 is the
                // first neighbor of int1
                let mut outflow_cells = 0;
                if self.network.neighbor_index(int1, int2).unwrap_or(0) == 0 {
                    for &other_cell in &nodes[int1].neighbors[0].triangles {
                        if self.network.triangles[other_cell].fracture == run.particle.fracture
                            && self.in_out_flow(&mut run.particle, other_cell, [int1, int2], middle)
                                > 0.0
                        {
                            outflow_cells += 1;
                        }
                    }
                }
                outflow_cells > 1
            };
            if moves_to_line || moves_out {
                run.particle.position = middle;
                self.check_new_cell(run);
                match run.particle.cell {
                    Some(new_cell) => {
                        self.cross_intersection(run, new_cell, [int1, int2], records_trajectory);
                    }
                    None if moves_to_line => {
                        warn!("A particle is lost at the end of an intersection");
                    }
                    None => {}
                }
            }
        }
    }

    /// At the crossing point: the path point and the matrix diffusion of
    /// the fracture, then the move into the next fracture
    fn cross_intersection(
        &mut self,
        run: &mut Run,
        previous_cell: usize,
        [int1, int2]: [usize; 2],
        records_trajectory: bool,
    ) {
        if records_trajectory {
            run.trajectory.push(self.path_point(&run.particle));
        }
        self.diffuse(run);
        self.across_intersection(&mut run.particle, previous_cell, int1, int2);
    }

    /// Positive if the velocity at the middle of the intersection segment
    /// points out of the intersection into `cell` (`InOutFlowCell`)
    fn in_out_flow(
        &self,
        particle: &mut Particle,
        cell: usize,
        edge: [usize; 2],
        middle: [f64; 2],
    ) -> f64 {
        let saved = (particle.position, particle.fracture, particle.cell);
        particle.position = middle;
        let product = self.outflow_product(particle, cell, edge).0;
        (particle.position, particle.fracture, particle.cell) = saved;
        product
    }

    /// Moves the particle into `cell` and gives the product which is
    /// negative if the flow goes out of the intersection into the cell, and
    /// the squared speed there (part of `AcrossIntersection`)
    fn outflow_product(
        &self,
        particle: &mut Particle,
        cell: usize,
        [int1, int2]: [usize; 2],
    ) -> (f64, f64, usize) {
        let nodes = &self.network.nodes;
        let triangle = &self.network.triangles[cell];
        // The last node which is not on an intersection; if all are (not
        // in DFNTrans), the node which is not on the edge
        let third = triangle
            .nodes
            .iter()
            .copied()
            .rev()
            .find(|&node| !matches!(nodes[node].type_number(), 2 | 12))
            .or_else(|| {
                triangle
                    .nodes
                    .iter()
                    .copied()
                    .find(|&node| node != int1 && node != int2)
            })
            .unwrap_or(triangle.nodes[0]);
        let [tnx, tny] = nodes[third].positions_2d[0];
        let [vintx, vinty] = nodes[int1].position_in(triangle.fracture);

        self.change_fracture(particle, cell);
        let weights = self.weights(particle, triangle.nodes);
        let mut velocity = [0.0; 2];
        for (corner, weight) in weights.iter().enumerate() {
            let corner_velocity = self.network.corner_velocity(cell, corner);
            velocity[0] += weight * corner_velocity[0];
            velocity[1] += weight * corner_velocity[1];
        }
        let [px, py] = particle.position;
        let side = (px - tnx) * (py - vinty) - (px - vintx) * (py - tny);
        let flow_side = velocity[0] * (py - vinty) - (px - vintx) * velocity[1];
        (
            flow_side * side,
            velocity[0] * velocity[0] + velocity[1] * velocity[1],
            nodes[third].fractures[0],
        )
    }

    /// The move through the intersection edge (int1, int2): the out-going
    /// cell from the routing rule (`AcrossIntersection`)
    fn across_intersection(
        &mut self,
        particle: &mut Particle,
        previous_cell: usize,
        int1: usize,
        int2: usize,
    ) {
        let previous_fracture = self.network.triangles[previous_cell].fracture;
        let Some(edge) = self.network.neighbor_index(int1, int2) else {
            warn!("Intersection edge {int1}-{int2} not found");
            return;
        };
        let mut cells = Vec::new();
        let mut previous_index = 0;
        for (index, &cell) in self.network.nodes[int1].neighbors[edge]
            .triangles
            .iter()
            .enumerate()
        {
            if cell == previous_cell {
                previous_index = index;
            }
            let (outflow_product, squared_speed, third_node_fracture) =
                self.outflow_product(particle, cell, [int1, int2]);
            cells.push(IntersectionCell {
                cell,
                fracture: self.network.triangles[cell].fracture,
                outflow_product,
                squared_speed,
                third_node_fracture,
            });
        }

        let winner = match self.input.routing {
            RoutingRule::CompleteMixing => {
                complete_mixing(&cells, previous_index, &mut self.random)
            }
            RoutingRule::StreamlineRouting => streamline_routing(
                &cells,
                previous_index,
                previous_cell,
                previous_fracture,
                &mut self.random,
            ),
        };
        let Some(winner) = winner else {
            // No flow out of the intersection
            particle.cell = None;
            return;
        };
        self.change_fracture(particle, cells[winner].cell);
        particle.location = CellLocation::CrossedIntersection;
        particle.previous_position = particle.position;
    }

    /// Moves the particle to the center of a cell (`Moving2Center`)
    pub fn move_to_center(&self, particle: &mut Particle, cell: usize) {
        let nodes = self.network.triangles[cell].nodes;
        let mut center = [0.0; 2];
        for node in nodes {
            let position = self.node_position(node, particle);
            center[0] += position[0];
            center[1] += position[1];
        }
        particle.position = [center[0] / 3.0, center[1] / 3.0];
        particle.cell = Some(cell);
        self.is_inside(particle, cell);
    }

    /// A particle out of the fracture moves to the center of a cell next to
    /// its previous cell, preferably an interior one (`Moving2NextCellBound`)
    fn move_from_fracture_edge(&self, particle: &mut Particle, previous_cell: Option<usize>) {
        let Some(previous_cell) = previous_cell else {
            return;
        };
        let nodes = &self.network.nodes;
        let [n1, n2, n3] = self.network.triangles[previous_cell].nodes;
        let center = if nodes[n1].type_number() == 0 {
            n1
        } else if nodes[n2].type_number() == 0 {
            n2
        } else {
            n3
        };
        let type_sum = |cell: usize| -> u32 {
            self.network.triangles[cell]
                .nodes
                .iter()
                .map(|&node| nodes[node].type_number())
                .sum()
        };

        let mut interior = None;
        let mut one_edge_node = None;
        let mut sum_twelve = None;
        let mut on_boundary = None;
        'search: for neighbor in &nodes[center].neighbors {
            for &cell in &neighbor.triangles {
                if self.network.triangles[cell].fracture != particle.fracture
                    || cell == previous_cell
                {
                    continue;
                }
                match type_sum(cell) {
                    0 => {
                        interior = Some(cell);
                        break 'search;
                    }
                    10 => one_edge_node = Some(cell),
                    12 => sum_twelve = Some(cell),
                    sum if sum > 15 => on_boundary = Some(cell),
                    _ => {}
                }
            }
        }
        if let Some(cell) = interior.or(one_edge_node).or(sum_twelve).or(on_boundary) {
            particle.cell = Some(cell);
            self.move_to_center(particle, cell);
        }
    }

    /// The exit point on the out-flow boundary: the crossing of the last
    /// step with the out-flow edge, and the time to it (`FinalPosition`)
    fn move_to_outflow_boundary(&self, particle: &mut Particle) {
        let nodes = &self.network.nodes;
        let Some(cell) = particle.cell else {
            return;
        };
        let [n1, n2, n3] = self.network.triangles[cell].nodes;
        let is_out = |node: usize| nodes[node].zone == Some(FlowZone::Outflow);
        // The selection of the two nodes of DFNTrans
        let mut first = None;
        let mut second = None;
        if is_out(n1) {
            first = Some(n1);
        }
        if is_out(n2) {
            if first.is_none() {
                second = Some(n2);
            } else {
                first = Some(n2);
            }
        }
        if is_out(n3) && second.is_none() {
            second = Some(n3);
        }

        let edge = match (first, second) {
            (Some(first), Some(second)) => Some((first, second, None)),
            (Some(center), None) | (None, Some(center)) => nodes[center]
                .neighbors
                .iter()
                .find(|neighbor| matches!(nodes[neighbor.node].type_number(), 212 | 210 | 202))
                .map(|neighbor| (center, neighbor.node, Some(neighbor.triangles[0]))),
            (None, None) => None,
        };
        let Some((first, second, new_cell)) = edge else {
            return;
        };

        let [cx1, cy1] = self.node_position(first, particle);
        let [cx2, cy2] = self.node_position(second, particle);
        let a_edge = cy2 - cy1;
        let b_edge = cx1 - cx2;
        let c_edge = a_edge * cx1 + b_edge * cy1;
        let [px1, py1] = particle.position;
        let [px2, py2] = particle.previous_position;
        let a_step = py2 - py1;
        let b_step = px1 - px2;
        let c_step = a_step * px1 + b_step * py1;
        let determinant = a_step * b_edge - a_edge * b_step;
        let crossing = [
            (b_edge * c_step - b_step * c_edge) / determinant,
            (a_step * c_edge - a_edge * c_step) / determinant,
        ];
        let distance = (crossing[0] - px1).hypot(crossing[1] - py1);
        let speed = particle.velocity[0].hypot(particle.velocity[1]);
        particle.previous_position = particle.position;
        particle.position = crossing;
        particle.time += distance / speed;
        if let Some(new_cell) = new_cell {
            particle.cell = Some(new_cell);
        }
    }

    /// Lagrangian time τ and β of the step from `previous` to `current`
    /// (`CalculateLagrangian`)
    fn lagrangian(&self, particle: &Particle, current: [f64; 3], previous: [f64; 3]) -> (f64, f64) {
        let velocity = self.velocity_3d(particle);
        let squared_distance = subtract(current, previous)
            .iter()
            .map(|d| d * d)
            .sum::<f64>();
        let squared_speed = velocity.iter().map(|v| v * v).sum::<f64>();
        if squared_distance <= 0.0 {
            return (0.0, 0.0);
        }
        let cell = particle.cell.expect("the particle has a cell");
        let aperture = self.network.nodes[self.network.triangles[cell].nodes[0]].aperture;
        let tau = (squared_distance / squared_speed).sqrt();
        let beta = squared_distance.sqrt() / (squared_speed.sqrt() * (aperture * 0.5));
        (tau, beta)
    }
}
