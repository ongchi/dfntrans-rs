//! Comparison with the C DFNTrans of dfnWorks. The files
//! `tests/data/reference/<case>.partime.dat` are its output for the cases
//! (see `scripts/compare_with_dfntrans.sh`). Each particle must take the
//! same path: the same random decisions at the intersection, a travel time
//! within 1e-4 (relative; DFNTrans has one density, dfntrans-rs the density
//! of each node) and the same exit point within 1e-3 m.

use std::path::Path;

use dfntrans::tracking::ParticleResult;
use dfntrans::{read_input, track_particles};

const TOLERANCE: f64 = 1e-4;

fn particles(case: &str) -> Vec<ParticleResult> {
    let input =
        read_input(&Path::new("tests/data/reference").join(format!("{case}.toml"))).unwrap();
    let flow = fehm::read_input(&input.flow_file).unwrap();
    let solution = fehm::solve(&flow).unwrap();
    track_particles(&input, &flow, &solution).unwrap().particles
}

/// The rows of `partime.dat`
fn reference(case: &str) -> Vec<Vec<f64>> {
    let path = Path::new("tests/data/reference").join(format!("{case}.partime.dat"));
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .skip(1)
        .map(|line| {
            line.split_whitespace()
                .map(|value| value.parse().unwrap())
                .collect()
        })
        .collect()
}

fn relative_difference(a: f64, b: f64) -> f64 {
    (a - b).abs() / a.abs().max(b.abs())
}

fn assert_same_paths(case: &str) {
    let rows = reference(case);
    let particles = particles(case);
    assert_eq!(rows.len(), particles.len());
    for (row, particle) in rows.iter().zip(&particles) {
        // Steps, flux weight, travel time, x, y, z, beta, length
        assert!(
            relative_difference(row[2], particle.travel_time) < TOLERANCE,
            "particle {}: {} {}",
            particle.index + 1,
            row[2],
            particle.travel_time
        );
        for axis in 0..3 {
            assert!(
                (row[3 + axis] - particle.final_position[axis]).abs() < 1e-3,
                "particle {}: {:?} {:?}",
                particle.index + 1,
                &row[3..6],
                particle.final_position
            );
        }
    }
}

#[test]
fn complete_mixing_gives_the_paths_of_dfntrans() {
    assert_same_paths("complete_mixing");
}

#[test]
fn streamline_routing_gives_the_paths_of_dfntrans() {
    assert_same_paths("streamline_routing");
}

#[test]
fn matrix_diffusion_gives_the_times_of_dfntrans() {
    let rows = reference("matrix_diffusion");
    let particles = particles("matrix_diffusion");
    assert_eq!(rows.len(), particles.len());
    for (row, particle) in rows.iter().zip(&particles) {
        // Steps, flux weight, advective time, total time, diffusion time, beta, length
        let diffusion = particle.diffusion.unwrap();
        for (expected, actual) in [
            (row[2], particle.travel_time),
            (row[3], diffusion.total_time),
        ] {
            assert!(
                relative_difference(expected, actual) < TOLERANCE,
                "particle {}: {expected} {actual}",
                particle.index + 1
            );
        }
    }
}
