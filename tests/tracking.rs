//! Checks of the particle tracking on the small meshes of fehm-rs.

use std::cell::RefCell;
use std::io::Write;
use std::path::Path;
use std::rc::Rc;

use dfntrans::network::PointType;
use dfntrans::{
    DfntransError, Input, Transport, build_input_from_text, read_input, track_particles,
    write_output,
};

fn run_input(input: &Input) -> Result<Transport, DfntransError> {
    let flow = fehm::read_input(&input.flow_file)?;
    let solution = fehm::solve(&flow)?;
    track_particles(input, &flow, &solution)
}

fn run(case: &str) -> (Input, Transport) {
    let input = read_input(Path::new("tests/data").join(case).as_path()).unwrap();
    let transport = run_input(&input).unwrap();
    (input, transport)
}

/// The input file of a case with other text
fn changed_input(case: &str, from: &str, to: &str) -> Input {
    let path = Path::new("tests/data").join(case);
    let text = std::fs::read_to_string(&path).unwrap();
    assert!(text.contains(from), "{from}");
    let mut input = build_input_from_text(&text.replace(from, to)).unwrap();
    input.flow_file = Path::new("tests/data").join(input.flow_file);
    input
}

#[test]
fn uniform_flow_gives_one_velocity_at_all_nodes() {
    let (_, transport) = run("single_fracture.toml");
    let speeds: Vec<f64> = transport
        .network
        .nodes
        .iter()
        .flat_map(|node| &node.sectors)
        .map(|sector| sector.velocity[0].hypot(sector.velocity[1]))
        .collect();
    let mean = speeds.iter().sum::<f64>() / speeds.len() as f64;
    for speed in speeds {
        // The density changes with the pressure by about 1e-5
        assert!((speed / mean - 1.0).abs() < 1e-4, "{speed} {mean}");
    }
}

#[test]
fn uniform_flow_gives_the_travel_time_of_the_mean_velocity() {
    let (input, transport) = run("single_fracture.toml");
    // Fracture 10 m long, 5 m high, aperture 1e-4 m
    let speed = transport.inflow_volume_flux / (5.0 * 1e-4) * input.time_unit.seconds();
    let expected = 10.0 / speed;
    assert_eq!(transport.exit_count, 20);
    for particle in &transport.particles {
        assert!(
            (particle.travel_time / expected - 1.0).abs() < 1e-3,
            "{} {expected}",
            particle.travel_time
        );
    }
}

#[test]
fn uniform_flow_keeps_the_height_of_a_particle() {
    let (_, transport) = run("single_fracture.toml");
    for particle in &transport.particles {
        // Particles at 20 equal distances on the west edge, z = 0..5
        let start = 5.0 * (particle.index as f64 + 0.5) / 20.0;
        assert!(
            (particle.final_position[2] - start).abs() < 1e-6,
            "{:?}",
            particle.final_position
        );
        assert!((particle.final_position[0] - 10.0).abs() < 1e-9);
    }
}

/// The single fracture case with particles at given points (y = 0 is the
/// plane of the fracture)
fn single_fracture_with_points(points: &str) -> Input {
    changed_input(
        "single_fracture.toml",
        "method = \"uniform_per_line\"\nparticles_per_line = 20",
        &format!("method = \"points\"\npoints = {points}\nparticles_per_point = 2"),
    )
}

#[test]
fn particles_at_points_take_the_travel_time_from_their_points() {
    let input = single_fracture_with_points("[[4.0, 0.0, 2.5], [7.5, 1e-3, 1.2]]");
    let transport = run_input(&input).unwrap();
    let speed = transport.inflow_volume_flux / (5.0 * 1e-4) * input.time_unit.seconds();
    let starts = [4.0, 4.0, 7.5, 7.5];
    assert_eq!(transport.exit_count, 4);
    // The error of the time steps is absolute: the tolerance of the uniform
    // flow test, on the time of the full length
    let tolerance = 1e-3 * 10.0 / speed;
    for particle in &transport.particles {
        let expected = (10.0 - starts[particle.index]) / speed;
        assert!(
            (particle.travel_time - expected).abs() < tolerance,
            "{} {expected}",
            particle.travel_time
        );
    }
}

#[test]
fn particles_at_points_keep_the_height_of_their_points() {
    let input = single_fracture_with_points("[[4.0, 0.0, 2.5], [7.5, 1e-3, 1.2]]");
    let transport = run_input(&input).unwrap();
    let heights = [2.5, 2.5, 1.2, 1.2];
    for particle in &transport.particles {
        assert!(
            (particle.final_position[2] - heights[particle.index]).abs() < 1e-6,
            "{:?}",
            particle.final_position
        );
    }
}

#[test]
fn seeding_point_far_from_the_mesh_is_an_error() {
    let input = single_fracture_with_points("[[4.0, 0.0, 2.5], [4.0, 3.0, 2.5]]");
    let result = run_input(&input);
    assert!(matches!(
        result,
        Err(DfntransError::SeedingPointOutsideMesh { index: 2, .. })
    ));
}

#[test]
fn particles_at_an_intersection_go_into_both_outgoing_fractures() {
    let (_, transport) = run("crossing_fractures.toml");
    assert_eq!(transport.exit_count, transport.particle_count);
    let east = transport
        .particles
        .iter()
        .filter(|particle| (particle.final_position[0] - 5.0).abs() < 1e-6)
        .count();
    let north = transport
        .particles
        .iter()
        .filter(|particle| (particle.final_position[1] - 5.0).abs() < 1e-6)
        .count();
    assert!(east > 0 && north > 0, "{east} {north}");
    assert_eq!(east + north, transport.particle_count);
}

#[test]
fn intersection_nodes_are_found_from_the_mesh() {
    let (_, transport) = run("crossing_fractures.toml");
    let count = |point_type: PointType| {
        transport
            .network
            .nodes
            .iter()
            .filter(|node| node.point_type == point_type)
            .count()
    };
    // The line x = y = 0, z = -5..5: 11 nodes, 2 of them at the ends
    assert_eq!(count(PointType::Interface), 9);
    assert_eq!(count(PointType::ExteriorInterface), 2);
    // A node on the line has two sectors in each fracture
    for node in &transport.network.nodes {
        let expected = if node.fractures.len() == 2 { 4 } else { 1 };
        assert_eq!(node.sectors.len(), expected);
    }
}

#[test]
fn plane_coordinates_give_the_point_back_in_3d() {
    let (_, transport) = run("crossing_fractures.toml");
    let network = &transport.network;
    for node in &network.nodes {
        for (&fracture, &position) in node.fractures.iter().zip(&node.positions_2d) {
            let back = network.fractures[fracture].to_space(position);
            for (coordinate, expected) in back.iter().zip(node.position) {
                assert!((coordinate - expected).abs() < 1e-12);
            }
        }
    }
}

#[test]
fn particles_cross_each_control_plane_once_in_order() {
    let input = changed_input(
        "single_fracture.toml",
        "[output]",
        "[control_planes]\nshape = \"plane\"\ndirection = \"x\"\nspacing = 1.0\n\n[output]",
    );
    let transport = run_input(&input).unwrap();
    for particle in &transport.particles {
        let crossings = &particle.control_crossings;
        // The start (x = 0), the planes x = 1..9 and the exit (x = 10)
        assert_eq!(crossings.len(), 11);
        for pair in crossings.windows(2) {
            assert!(pair[1].point.time > pair[0].point.time);
            assert!(pair[1].point.position[0] > pair[0].point.position[0]);
        }
        for (plane, crossing) in crossings.iter().enumerate().skip(1).take(9) {
            // The first position at or behind the plane, one step of about 5 mm
            let x = crossing.point.position[0];
            assert!(x >= plane as f64 && x < plane as f64 + 0.05, "{plane} {x}");
        }
    }
}

#[test]
fn matrix_diffusion_adds_time() {
    let input = changed_input(
        "single_fracture.toml",
        "[output]",
        "[tdrw]\nmatrix_porosity = 0.02\nmatrix_diffusivity = 1e-11\n\n[output]",
    );
    let transport = run_input(&input).unwrap();
    for particle in &transport.particles {
        let diffusion = particle.diffusion.unwrap();
        assert!(diffusion.diffusion_time > 0.0);
        assert!(
            (diffusion.total_time - particle.travel_time - diffusion.diffusion_time).abs()
                < 1e-12 * diffusion.total_time
        );
    }
}

#[test]
fn every_seeding_method_places_particles_which_mostly_exit() {
    for seeding in [
        "method = \"uniform_per_line\"\nparticles_per_line = 5",
        "method = \"equidistant\"\nparticles_per_line = 5",
        "method = \"random\"\nparticle_count = 20",
        "method = \"flux_weighted\"\nparticle_count = 20",
        "method = \"at_nodes\"\nparticles_per_node = 2",
    ] {
        let input = changed_input(
            "crossing_fractures.toml",
            "method = \"flux_weighted\"\nparticle_count = 200",
            seeding,
        );
        let transport = run_input(&input).unwrap();
        assert!(transport.particle_count > 0, "{seeding}");
        // Random cells can be in a corner without flow
        assert!(
            10 * transport.exit_count >= 9 * transport.particle_count,
            "{seeding}"
        );
        let total_weight: f64 = transport
            .particles
            .iter()
            .map(|particle| particle.flux_weight)
            .sum();
        assert!(total_weight > 0.0, "{seeding}");
    }
}

#[test]
fn unknown_boundary_name_is_an_error() {
    let input = changed_input(
        "single_fracture.toml",
        "inflow = [\"west\"]",
        "inflow = [\"south\"]",
    );
    assert!(
        matches!(run_input(&input), Err(DfntransError::UnknownBoundary { boundary }) if boundary == "south")
    );
}

/// The output files of a run, in memory: (name, text)
fn output_files(input: &Input, transport: &Transport) -> Vec<(String, String)> {
    let buffers: Rc<RefCell<Vec<(String, Buffer)>>> = Rc::default();
    write_output(input, transport, &mut |name| {
        let buffer = Rc::new(RefCell::new(Vec::new()));
        buffers
            .borrow_mut()
            .push((name.to_string(), Rc::clone(&buffer)));
        Ok(Box::new(MemoryFile(buffer)))
    })
    .unwrap();
    buffers
        .borrow()
        .iter()
        .map(|(name, buffer)| {
            (
                name.clone(),
                String::from_utf8(buffer.borrow().clone()).unwrap(),
            )
        })
        .collect()
}

/// The content of a file in memory
type Buffer = Rc<RefCell<Vec<u8>>>;

/// A writer into a shared buffer
struct MemoryFile(Buffer);

impl Write for MemoryFile {
    fn write(&mut self, buffer: &[u8]) -> std::io::Result<usize> {
        self.0.borrow_mut().extend_from_slice(buffer);
        Ok(buffer.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

#[test]
fn output_files_have_one_entry_for_each_particle() {
    let (input, transport) = run("single_fracture.toml");
    let files = output_files(&input, &transport);
    let names: Vec<&str> = files.iter().map(|(name, _)| name.as_str()).collect();
    assert_eq!(
        names,
        ["travel_times.dat", "inflow_flux.dat", "trajectories.vtp"]
    );

    // A header line and one line for each of the 20 particles
    assert_eq!(files[0].1.lines().count(), 21);
    let trajectories = &files[2].1;
    assert!(trajectories.contains(r#"NumberOfLines="20""#));
    let point_count: usize = transport
        .particles
        .iter()
        .map(|particle| particle.trajectory.len())
        .sum();
    assert!(trajectories.contains(&format!(r#"NumberOfPoints="{point_count}""#)));
    assert!(trajectories.trim_end().ends_with("</VTKFile>"));
}

/// Intersections which end next to each other at the outline of a fracture:
/// the outline neighbors of a node on an intersection are on other
/// intersections. The halves of the node are divided at its intersection
/// edge, not at an outline edge, and a sector of faces on one line gets no
/// singular fit (it gave about 1e40).
#[test]
fn sector_velocities_stay_near_the_darcy_velocities_where_intersections_end_at_an_outline() {
    let input = read_input(Path::new("tests/data/outline_intersections.toml")).unwrap();
    let flow = fehm::read_input(&input.flow_file).unwrap();
    let solution = fehm::solve(&flow).unwrap();
    let largest_darcy_speed =
        fehm::velocity::calculate_darcy_velocities(&flow, &solution.pressures)
            .iter()
            .map(|velocity| {
                velocity
                    .iter()
                    .map(|component| component * component)
                    .sum::<f64>()
                    .sqrt()
            })
            .fold(0.0, f64::max);
    let transport = track_particles(&input, &flow, &solution).unwrap();
    // Seconds and porosity 1: a sector velocity has the units of the Darcy velocity
    let largest_sector_speed = transport
        .network
        .nodes
        .iter()
        .flat_map(|node| &node.sectors)
        .map(|sector| sector.velocity[0].hypot(sector.velocity[1]))
        .fold(0.0, f64::max);
    assert!(
        largest_sector_speed < 100.0 * largest_darcy_speed,
        "largest sector speed {largest_sector_speed:e}, largest Darcy speed {largest_darcy_speed:e}"
    );
}
