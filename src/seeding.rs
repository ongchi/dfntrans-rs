//! Initial positions and flux weights of the particles (DFNTrans
//! `InitialPartPositions.c`).
//!
//! DFNTrans places the particles on the straight segment between the first
//! and the last in-flow node of each fracture. The port uses the in-flow
//! lines: the chains of mesh edges between two in-flow nodes in one
//! fracture, of any shape.

use std::collections::BTreeMap;

use crate::DfntransError;
use crate::io::input::Seeding;
use crate::network::{FlowZone, Network};
use crate::tracking::{Particle, Tracker};

/// A chain of edges between in-flow nodes in one fracture
#[derive(Clone, Debug)]
pub struct InflowLine {
    pub fracture: usize,
    pub nodes: Vec<usize>,
}

impl InflowLine {
    fn points(&self, network: &Network) -> Vec<[f64; 2]> {
        self.nodes
            .iter()
            .map(|&node| network.nodes[node].position_in(self.fracture))
            .collect()
    }

    /// Distance along the line of each node
    fn distances(&self, network: &Network) -> Vec<f64> {
        let points = self.points(network);
        let mut distances = vec![0.0];
        for pair in points.windows(2) {
            let last = distances[distances.len() - 1];
            distances.push(last + (pair[1][0] - pair[0][0]).hypot(pair[1][1] - pair[0][1]));
        }
        distances
    }

    fn length(&self, network: &Network) -> f64 {
        self.distances(network).last().copied().unwrap_or(0.0)
    }

    /// The point at a distance along the line
    fn point_at(&self, network: &Network, distance: f64) -> [f64; 2] {
        let points = self.points(network);
        let distances = self.distances(network);
        for segment in 0..points.len().saturating_sub(1) {
            let (start, end) = (distances[segment], distances[segment + 1]);
            if distance <= end || segment + 2 == points.len() {
                let fraction = if end > start {
                    ((distance - start) / (end - start)).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                let (a, b) = (points[segment], points[segment + 1]);
                return [
                    a[0] + fraction * (b[0] - a[0]),
                    a[1] + fraction * (b[1] - a[1]),
                ];
            }
        }
        points[0]
    }
}

/// The in-flow lines of all fractures, in increasing order of the
/// fracture and of the first node
pub fn inflow_lines(network: &Network) -> Vec<InflowLine> {
    let is_inflow = |node: usize| network.nodes[node].zone == Some(FlowZone::Inflow);
    let mut lines = Vec::new();
    for fracture in 0..network.fractures.len() {
        // The in-flow edges of this fracture
        let mut adjacent: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for &node in &network.inflow_nodes {
            for neighbor in &network.nodes[node].neighbors {
                if neighbor.node > node
                    && is_inflow(neighbor.node)
                    && neighbor
                        .triangles
                        .iter()
                        .any(|&triangle| network.triangles[triangle].fracture == fracture)
                {
                    adjacent.entry(node).or_default().push(neighbor.node);
                    adjacent.entry(neighbor.node).or_default().push(node);
                }
            }
        }

        let mut used: Vec<[usize; 2]> = Vec::new();
        let is_used =
            |used: &[[usize; 2]], a: usize, b: usize| used.contains(&[a.min(b), a.max(b)]);
        // Lines from the ends and branch points first, then closed lines
        let starts: Vec<usize> = adjacent
            .iter()
            .filter(|(_, others)| others.len() != 2)
            .map(|(&node, _)| node)
            .chain(adjacent.keys().copied())
            .collect();
        for start in starts {
            while let Some(&next) = adjacent[&start]
                .iter()
                .find(|&&other| !is_used(&used, start, other))
            {
                let mut nodes = vec![start];
                let (mut previous, mut current) = (start, next);
                used.push([previous.min(current), previous.max(current)]);
                nodes.push(current);
                while adjacent[&current].len() == 2 {
                    let Some(&following) = adjacent[&current]
                        .iter()
                        .find(|&&other| other != previous && !is_used(&used, current, other))
                    else {
                        break;
                    };
                    used.push([current.min(following), current.max(following)]);
                    nodes.push(following);
                    (previous, current) = (current, following);
                }
                lines.push(InflowLine { fracture, nodes });
            }
        }
    }
    lines
}

/// The particles of a seeding method
pub(crate) fn seed_particles(
    tracker: &mut Tracker,
    seeding: &Seeding,
    inflow_volume_flux: f64,
) -> Result<Vec<Particle>, DfntransError> {
    let network = tracker.network;
    let particles = match *seeding {
        Seeding::UniformPerLine {
            particles_per_line,
            weight_by_flux,
        } => {
            let mut particles = Vec::new();
            for line in inflow_lines(network) {
                let length = line.length(network);
                for index in 0..particles_per_line {
                    let distance = length * (index as f64 + 0.5) / particles_per_line as f64;
                    particles.push(Particle::new(
                        line.point_at(network, distance),
                        line.fracture,
                        0.0,
                    ));
                }
            }
            if weight_by_flux {
                set_flux_weights(tracker, &mut particles);
            }
            particles
        }
        Seeding::Equidistant {
            particles_per_line,
            weight_by_flux,
        } => {
            let lines = inflow_lines(network);
            let total_length: f64 = lines.iter().map(|line| line.length(network)).sum();
            let spacing = total_length / (particles_per_line * lines.len()) as f64;
            let mut particles = Vec::new();
            for line in &lines {
                let length = line.length(network);
                let count = (length / spacing) as usize;
                let distances: Vec<f64> = if count < 2 {
                    vec![length / 2.0]
                } else {
                    (0..count)
                        .map(|index| spacing * (index as f64 + 0.5))
                        .collect()
                };
                for distance in distances {
                    particles.push(Particle::new(
                        line.point_at(network, distance),
                        line.fracture,
                        0.0,
                    ));
                }
            }
            if weight_by_flux {
                set_flux_weights(tracker, &mut particles);
            }
            particles
        }
        Seeding::Random { particle_count } => random_particles(tracker, particle_count)?,
        Seeding::FluxWeighted { particle_count } => {
            flux_weighted_particles(network, particle_count, inflow_volume_flux)?
        }
        Seeding::AtNodes { particles_per_node } => {
            let mut particles = Vec::new();
            for &node in &network.inflow_nodes {
                let node = &network.nodes[node];
                for _ in 0..particles_per_node {
                    particles.push(Particle::new(
                        node.positions_2d[0],
                        node.fractures[0],
                        node.aperture,
                    ));
                }
            }
            normalize_weights(&mut particles);
            particles
        }
    };
    Ok(particles)
}

fn normalize_weights(particles: &mut [Particle]) {
    let total: f64 = particles.iter().map(|particle| particle.flux_weight).sum();
    if total > 0.0 {
        for particle in particles {
            particle.flux_weight /= total;
        }
    }
}

/// Weight of each particle: the in-flow flux of the in-flow nodes of its
/// cell, interpolated (`FlowInWeight`)
fn set_flux_weights(tracker: &Tracker, particles: &mut [Particle]) {
    let network = tracker.network;
    let absolute_flux = |node: usize| -> f64 {
        network.nodes[node]
            .neighbors
            .iter()
            .map(|neighbor| neighbor.mass_flux.abs())
            .sum()
    };
    for particle in particles.iter_mut() {
        particle.flux_weight = 0.0;
        if !tracker.find_initial_cell(particle) {
            continue;
        }
        let Some(cell) = particle.cell else {
            continue;
        };
        let inflow_corners: Vec<usize> = (0..3)
            .filter(|&corner| {
                network.nodes[network.triangles[cell].nodes[corner]].zone == Some(FlowZone::Inflow)
            })
            .collect();
        let node_of = |corner: usize| network.triangles[cell].nodes[corner];
        particle.flux_weight = match inflow_corners[..] {
            [only] => absolute_flux(node_of(only)),
            [first, .., last] => {
                particle.weights[first] * absolute_flux(node_of(first))
                    + particle.weights[last] * absolute_flux(node_of(last))
            }
            [] => 0.0,
        };
    }
    normalize_weights(particles);
}

/// Particles at the centers of random cells which have no out-flow node,
/// weighted by aperture (`init_random`)
fn random_particles(tracker: &mut Tracker, count: usize) -> Result<Vec<Particle>, DfntransError> {
    let network = tracker.network;
    let is_allowed = |cell: usize| {
        network.triangles[cell]
            .nodes
            .iter()
            .all(|&node| network.nodes[node].zone != Some(FlowZone::Outflow))
    };
    if !(0..network.triangles.len().saturating_sub(1)).any(is_allowed) {
        return Err(DfntransError::NoParticles(
            "all cells (but the last one) have a node on the out-flow boundary".to_string(),
        ));
    }

    let cell_count = network.triangles.len();
    let mut particles = Vec::with_capacity(count);
    while particles.len() < count {
        // As DFNTrans: the cell number from 1, and 0 is no cell; thus the
        // last cell is never selected
        let number = (tracker.random.next_f64() * cell_count as f64) as usize;
        let Some(cell) = number.checked_sub(1) else {
            continue;
        };
        if !is_allowed(cell) {
            continue;
        }
        let triangle = &network.triangles[cell];
        let mut particle = Particle::new(
            [0.0; 2],
            triangle.fracture,
            network.nodes[triangle.nodes[0]].aperture,
        );
        tracker.move_to_center(&mut particle, cell);
        particles.push(particle);
    }
    normalize_weights(&mut particles);
    Ok(particles)
}

/// Particles on the in-flow lines in proportion to the in-flow volume flux
/// of each node: at least one for each node, and one for each part
/// `total flux / count` of its flux; each particle has the weight
/// `total flux / count` (`init_fluxw`, `InitParticles_flux`)
fn flux_weighted_particles(
    network: &Network,
    count: usize,
    inflow_volume_flux: f64,
) -> Result<Vec<Particle>, DfntransError> {
    if inflow_volume_flux.is_nan() || inflow_volume_flux <= 0.0 {
        return Err(DfntransError::NoParticles(format!(
            "the in-flow volume flux is {inflow_volume_flux:e} m³/s, not above 0"
        )));
    }
    let weight = inflow_volume_flux / count as f64;
    let lines = inflow_lines(network);

    // The number of lines of each (node, fracture): the flux of a node is
    // divided among them
    let mut line_counts: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for line in &lines {
        for &node in &line.nodes {
            *line_counts.entry((node, line.fracture)).or_default() += 1;
        }
    }

    let mut particles = Vec::new();
    for line in &lines {
        let distances = line.distances(network);
        let last = line.nodes.len() - 1;
        for (index, &node) in line.nodes.iter().enumerate() {
            let flux = fracture_mass_flux(network, node, line.fracture)
                / line_counts[&(node, line.fracture)] as f64
                / network.nodes[node].density;
            let particle_count = ((flux.abs() / weight).ceil() as usize).max(1);
            let start = if index == 0 {
                0.0
            } else {
                (distances[index - 1] + distances[index]) / 2.0
            };
            let end = if index == last {
                distances[last]
            } else {
                (distances[index] + distances[index + 1]) / 2.0
            };
            for particle in 0..particle_count {
                let distance =
                    start + (end - start) * (particle as f64 + 0.5) / particle_count as f64;
                particles.push(Particle::new(
                    line.point_at(network, distance),
                    line.fracture,
                    weight,
                ));
            }
        }
    }
    Ok(particles)
}

/// The mass flux out of a node into one fracture: the flux of each edge,
/// times the part of the triangles of the edge which are in the fracture
fn fracture_mass_flux(network: &Network, node: usize, fracture: usize) -> f64 {
    network.nodes[node]
        .neighbors
        .iter()
        .map(|neighbor| {
            let in_fracture = neighbor
                .triangles
                .iter()
                .filter(|&&triangle| network.triangles[triangle].fracture == fracture)
                .count();
            neighbor.mass_flux * in_fracture as f64 / neighbor.triangles.len().max(1) as f64
        })
        .sum()
}
