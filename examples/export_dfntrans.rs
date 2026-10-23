//! Writes the input files of the C DFNTrans of dfnWorks for a case, and
//! the output of dfntrans-rs for the same case, for a comparison of the two.
//!
//! Usage: cargo run --release --example export_dfntrans -- <input.toml> <folder>
//!
//! DFNTrans needs the nodes sorted by their first fracture and the cells by
//! their fracture, thus the mesh is sorted first; dfntrans-rs runs on the
//! sorted mesh too. `<folder>/rust` gets the output of dfntrans-rs, and
//! `<folder>/c` the DFNTrans files: `control.dat` (run `DFNTrans
//! control.dat` in that folder), `params.txt`, `full_mesh.inp`,
//! `full_mesh.stor`, `run.fin`, `boundary.zone` and `aperture.dat`.
//!
//! DFNTrans has one fluid density; it gets the mean of the node densities.
//! Only the seeding methods `at_nodes` and `random` give the same initial
//! positions in both codes.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::process::ExitCode;

use dfntrans::io::input::{Input, RoutingRule, Seeding, TimeUnit};
use dfntrans::network::{Network, PI};
use dfntrans::{read_input, track_particles, write_output_to_folder};

fn main() -> ExitCode {
    let arguments: Vec<String> = std::env::args().collect();
    let [_, input_file, folder] = arguments.as_slice() else {
        eprintln!("Usage: export_dfntrans <input.toml> <folder>");
        return ExitCode::FAILURE;
    };
    match run(Path::new(input_file), Path::new(folder)) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("{error}");
            ExitCode::FAILURE
        }
    }
}

fn run(input_file: &Path, folder: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let input = read_input(input_file)?;
    if !matches!(
        input.seeding,
        Seeding::AtNodes { .. } | Seeding::Random { .. }
    ) {
        return Err("only the seeding methods at_nodes and random have the same initial positions in DFNTrans".into());
    }
    let mut flow = fehm::read_input(&input.flow_file)?;
    flow.mesh = sorted_mesh(&flow.mesh);
    let solution = fehm::solve(&flow)?;
    let transport = track_particles(&input, &flow, &solution)?;
    write_output_to_folder(&input, &transport, &folder.join("rust"))?;

    let network = &transport.network;
    if let Some(node) = network
        .nodes
        .iter()
        .position(|node| node.fractures.len() > 2)
    {
        return Err(format!(
            "node {node} is on more than 2 fractures; DFNTrans does not support it"
        )
        .into());
    }

    let c_folder = folder.join("c");
    std::fs::create_dir_all(&c_folder)?;
    let create = |name: &str| -> std::io::Result<BufWriter<File>> {
        Ok(BufWriter::new(File::create(c_folder.join(name))?))
    };
    write_params(&mut create("params.txt")?, network)?;
    write_inp(&mut create("full_mesh.inp")?, network)?;
    let edge_count = write_stor(&mut create("full_mesh.stor")?, network)?;
    write_fin(&mut create("run.fin")?, network, edge_count)?;
    write_zones(&mut create("boundary.zone")?, network)?;
    write_apertures(&mut create("aperture.dat")?, network, &flow.mesh)?;
    let density =
        network.nodes.iter().map(|node| node.density).sum::<f64>() / network.nodes.len() as f64;
    write_control(&mut create("control.dat")?, &input, density)?;
    println!(
        "Wrote {} and {}",
        folder.join("rust").display(),
        c_folder.display()
    );
    Ok(())
}

/// The mesh with the nodes sorted by their first fracture and the
/// triangles by their fracture (the fractures in increasing order of `fracture_id`)
fn sorted_mesh(mesh: &fehm::mesh::Mesh) -> fehm::mesh::Mesh {
    let mut first_fracture = vec![i64::MAX; mesh.points.len()];
    for (triangle, &id) in mesh.triangles.iter().zip(&mesh.fracture_ids) {
        for &node in triangle {
            first_fracture[node] = first_fracture[node].min(id);
        }
    }
    let mut node_order: Vec<usize> = (0..mesh.points.len()).collect();
    node_order.sort_by_key(|&node| (first_fracture[node], node));
    let mut new_index = vec![0; mesh.points.len()];
    for (new, &old) in node_order.iter().enumerate() {
        new_index[old] = new;
    }
    let mut cell_order: Vec<usize> = (0..mesh.triangles.len()).collect();
    cell_order.sort_by_key(|&cell| (mesh.fracture_ids[cell], cell));

    fehm::mesh::Mesh {
        points: node_order.iter().map(|&old| mesh.points[old]).collect(),
        triangles: cell_order
            .iter()
            .map(|&cell| mesh.triangles[cell].map(|node| new_index[node]))
            .collect(),
        fracture_ids: cell_order
            .iter()
            .map(|&cell| mesh.fracture_ids[cell])
            .collect(),
        apertures: cell_order
            .iter()
            .map(|&cell| mesh.apertures[cell])
            .collect(),
        permeabilities: cell_order
            .iter()
            .map(|&cell| mesh.permeabilities[cell])
            .collect(),
        boundaries: mesh
            .boundaries
            .iter()
            .map(|boundary| {
                let mut nodes: Vec<usize> =
                    boundary.nodes.iter().map(|&node| new_index[node]).collect();
                nodes.sort_unstable();
                fehm::mesh::BoundaryNodes {
                    name: boundary.name.clone(),
                    nodes,
                }
            })
            .collect(),
    }
}

/// `params.txt`: the rotation angle (degrees) and axes of each fracture
fn write_params(writer: &mut impl Write, network: &Network) -> std::io::Result<()> {
    writeln!(writer, "{}", network.fractures.len())?;
    for _ in 0..7 {
        writeln!(writer, "0")?;
    }
    for (index, fracture) in network.fractures.iter().enumerate() {
        let degrees = fracture.angle * 180.0 / PI;
        let [x, y] = fracture.rotation_axis();
        writeln!(
            writer,
            "{} {degrees:.9} {x:.9} {y:.9} 0 {:.9} {:.9} 0 0",
            index + 1,
            -x,
            -y
        )?;
    }
    writer.flush()
}

/// `full_mesh.inp` (AVS UCD): points, triangles and the node arrays imt1 and itp1
fn write_inp(writer: &mut impl Write, network: &Network) -> std::io::Result<()> {
    writeln!(
        writer,
        "{} {} 3 0 0",
        network.nodes.len(),
        network.triangles.len()
    )?;
    for (index, node) in network.nodes.iter().enumerate() {
        let [x, y, z] = node.position;
        writeln!(writer, "{} {x:.17e} {y:.17e} {z:.17e}", index + 1)?;
    }
    for (index, triangle) in network.triangles.iter().enumerate() {
        let [a, b, c] = triangle.nodes.map(|node| node + 1);
        writeln!(
            writer,
            "{} {} tri {a} {b} {c}",
            index + 1,
            triangle.fracture + 1
        )?;
    }
    // DFNTrans skips the rest of the line after itp1, thus a third column
    writeln!(writer, "3 1 1 1")?;
    writeln!(writer, "imt1, integer")?;
    writeln!(writer, "itp1, integer")?;
    writeln!(writer, "isn1, integer")?;
    for (index, node) in network.nodes.iter().enumerate() {
        writeln!(
            writer,
            "{} {} {} 0",
            index + 1,
            node.fractures[0] + 1,
            node.base_type_number()
        )?;
    }
    writer.flush()
}

/// The neighbors of a node with the node itself, in increasing order:
/// (node, area coefficient, mass flux)
fn stor_row(network: &Network, node: usize) -> Vec<(usize, f64, f64)> {
    let mut row: Vec<(usize, f64, f64)> = network.nodes[node]
        .neighbors
        .iter()
        .map(|neighbor| {
            (
                neighbor.node,
                -neighbor.area_coefficient,
                neighbor.mass_flux,
            )
        })
        .collect();
    row.push((node, 0.0, 0.0));
    row.sort_by_key(|&(other, _, _)| other);
    row
}

/// `full_mesh.stor` (FEHM sparse matrix format); returns the number of entries
fn write_stor(writer: &mut impl Write, network: &Network) -> std::io::Result<usize> {
    let node_count = network.nodes.len();
    let rows: Vec<Vec<(usize, f64, f64)>> = (0..node_count)
        .map(|node| stor_row(network, node))
        .collect();
    let entry_count: usize = rows.iter().map(Vec::len).sum();
    let largest_row = rows.iter().map(Vec::len).max().unwrap_or(0);

    writeln!(writer, "fehmstor ascir8i4 dfntrans-rs export")?;
    writeln!(writer, "dfntrans-rs")?;
    writeln!(
        writer,
        "{entry_count} {node_count} {} 1 {largest_row}",
        entry_count + node_count + 1
    )?;
    write_values(
        writer,
        network
            .nodes
            .iter()
            .map(|node| format!("{:.17e}", node.volume)),
    )?;
    let mut pointer = node_count + 1;
    let mut pointers = vec![pointer];
    for row in &rows {
        pointer += row.len();
        pointers.push(pointer);
    }
    write_values(writer, pointers.iter().map(ToString::to_string))?;
    write_values(
        writer,
        rows.iter()
            .flatten()
            .map(|&(other, _, _)| (other + 1).to_string()),
    )?;
    write_values(writer, (1..=entry_count).map(|entry| entry.to_string()))?;
    write_values(writer, (0..=node_count).map(|_| "0".to_string()))?;
    write_values(writer, (0..node_count).map(|node| (node + 1).to_string()))?;
    write_values(
        writer,
        rows.iter()
            .flatten()
            .map(|&(_, area, _)| format!("{area:.17e}")),
    )?;
    writer.flush()?;
    Ok(entry_count)
}

/// `run.fin`: the pressures and the fluxes of each matrix entry
fn write_fin(
    writer: &mut impl Write,
    network: &Network,
    entry_count: usize,
) -> std::io::Result<()> {
    for line in ["dfntrans-rs export", "fehm-rs", "steady", "0"] {
        writeln!(writer, "{line}")?;
    }
    writeln!(writer, "pressure")?;
    write_values(
        writer,
        network
            .nodes
            .iter()
            .map(|node| format!("{:.8e}", node.pressure)),
    )?;
    writeln!(writer, "flux")?;
    writeln!(writer, "{entry_count}")?;
    let fluxes = (0..network.nodes.len())
        .flat_map(|node| stor_row(network, node).into_iter().map(|(_, _, flux)| flux));
    write_values(writer, fluxes.map(|flux| format!("{flux:.8e}")))?;
    writer.flush()
}

/// `boundary.zone`: zone 1 is the in-flow boundary, zone 2 the out-flow boundary
fn write_zones(writer: &mut impl Write, network: &Network) -> std::io::Result<()> {
    writeln!(writer, "zone")?;
    for (number, name, nodes) in [
        (1, "inflow", &network.inflow_nodes),
        (2, "outflow", &network.outflow_nodes),
    ] {
        writeln!(writer, "{number:05}  {name}")?;
        writeln!(writer, "nnum")?;
        writeln!(writer, "{}", nodes.len())?;
        write_values(writer, nodes.iter().map(|node| (node + 1).to_string()))?;
    }
    writeln!(writer, "stop")?;
    writer.flush()
}

/// `aperture.dat`: the aperture of each fracture (DFNTrans `aperture_type: frac`)
fn write_apertures(
    writer: &mut impl Write,
    network: &Network,
    mesh: &fehm::mesh::Mesh,
) -> std::io::Result<()> {
    writeln!(writer, "aperture")?;
    for fracture in 0..network.fractures.len() {
        let aperture = network
            .triangles
            .iter()
            .zip(&mesh.apertures)
            .find(|(triangle, _)| triangle.fracture == fracture)
            .map_or(0.0, |(_, &aperture)| aperture);
        // DFNTrans: fracture number = -(first value) - 6
        writeln!(writer, "-{} 0 0 {aperture:.17e}", fracture + 7)?;
    }
    writer.flush()
}

fn write_control(writer: &mut impl Write, input: &Input, density: f64) -> std::io::Result<()> {
    let yes_no = |value: bool| if value { "yes" } else { "no" };
    let time_unit = match input.time_unit {
        TimeUnit::Seconds => "seconds",
        TimeUnit::Minutes => "minutes",
        TimeUnit::Hours => "hours",
        TimeUnit::Days => "days",
        TimeUnit::Years => "years",
    };
    let (is_random, random_count, is_at_nodes, node_count) = match input.seeding {
        Seeding::Random { particle_count } => (true, particle_count, false, 1),
        Seeding::AtNodes { particles_per_node } => (false, 1, true, particles_per_node),
        _ => unreachable!("checked in run"),
    };
    let mut lines = vec![
        "param: params.txt".to_string(),
        "inp: full_mesh.inp".to_string(),
        "stor: full_mesh.stor".to_string(),
        "boundary: boundary.zone".to_string(),
        "in-flow-boundary: 1".to_string(),
        "out-flow-boundary: 2".to_string(),
        "PFLOTRAN: no".to_string(),
        "FEHM: yes".to_string(),
        "FEHM_fin: run.fin".to_string(),
        "out_dir: output".to_string(),
        "out_path: trajectories".to_string(),
        "out_time: partime.dat".to_string(),
        "out_grid: no".to_string(),
        "out_3dflow: yes".to_string(),
        "out_curv: no".to_string(),
        "out_avs: no".to_string(),
        "out_traj: no".to_string(),
        format!(
            "allparticles_output: {}",
            yes_no(input.output.all_particles)
        ),
        "ControlPlane: no".to_string(),
        "init_nf: no".to_string(),
        "init_eqd: no".to_string(),
        "init_oneregion: no".to_string(),
        format!("init_random: {}", yes_no(is_random)),
        format!("in_randpart: {random_count}"),
        "init_matrix: no".to_string(),
        "init_fluxw: no".to_string(),
        format!("init_well: {}", yes_no(is_at_nodes)),
        format!("init_nodepart: {node_count}"),
        "aperture: yes".to_string(),
        "aperture_type: frac".to_string(),
        "aperture_file: aperture.dat".to_string(),
        "thickness: 1.0".to_string(),
        format!("porosity: {}", input.porosity),
        format!("density: {density:.17e}"),
        "satur: 1.0".to_string(),
        format!("timesteps: {}", input.max_time_steps),
        format!("time_units: {time_unit}"),
        "flux_weight: no".to_string(),
        format!("seed: {}", input.seed),
        format!(
            "streamline_routing: {}",
            yes_no(input.routing == RoutingRule::StreamlineRouting)
        ),
    ];
    if let Some(diffusion) = &input.matrix_diffusion {
        lines.push("tdrw: yes".to_string());
        lines.push(format!("tdrw_porosity: {}", diffusion.matrix_porosity));
        lines.push(format!(
            "tdrw_diffcoeff: {:e}",
            diffusion.matrix_diffusivity
        ));
    } else {
        lines.push("tdrw: no".to_string());
    }
    lines.push("END".to_string());
    for line in lines {
        writeln!(writer, "{line}")?;
    }
    writer.flush()
}

/// Values, 10 on a line
fn write_values(
    writer: &mut impl Write,
    values: impl Iterator<Item = String>,
) -> std::io::Result<()> {
    let values: Vec<String> = values.collect();
    for line in values.chunks(10) {
        writeln!(writer, "{}", line.join(" "))?;
    }
    Ok(())
}
