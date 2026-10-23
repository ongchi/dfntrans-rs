//! The output files of a run.

use std::fmt::Display;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::DfntransError;
use crate::io::input::Input;
use crate::network::Sector;
use crate::tracking::{DiffusionTimes, ParticleResult, PathPoint, Transport};

/// Travel time, final position, β and length of each particle (DFNTrans `partime.dat`)
pub const TRAVEL_TIMES_FILE: &str = "travel_times.dat";
/// The total in-flow volume flux (DFNTrans `inputflux_m3s`)
pub const INFLOW_FLUX_FILE: &str = "inflow_flux.dat";
/// The path of each particle (with `output.trajectories`)
pub const TRAJECTORIES_FILE: &str = "trajectories.vtp";
/// The crossings of the control planes or cylinders (with `[control_planes]`)
pub const CONTROL_PLANES_FILE: &str = "control_planes.vtp";
/// The reconstructed velocities (with `output.velocities`)
pub const VELOCITIES_FILE: &str = "velocities.vtp";

/// Writes the output files into a folder. The folder is made if it does not exist.
pub fn write_output_to_folder(
    input: &Input,
    transport: &Transport,
    folder: &Path,
) -> Result<(), DfntransError> {
    std::fs::create_dir_all(folder)?;
    write_output(input, transport, &mut |name| {
        Ok(Box::new(BufWriter::new(File::create(folder.join(name))?)))
    })
}

/// Writes the output files. `create_file` gives the writer of each file from its name.
pub fn write_output(
    input: &Input,
    transport: &Transport,
    create_file: &mut dyn FnMut(&str) -> std::io::Result<Box<dyn Write>>,
) -> Result<(), DfntransError> {
    let mut writer = create_file(TRAVEL_TIMES_FILE)?;
    write_travel_times(&mut writer, input, transport)?;
    writer.flush()?;

    let mut writer = create_file(INFLOW_FLUX_FILE)?;
    writeln!(writer, "{:12.5e}", transport.inflow_volume_flux)?;
    writer.flush()?;

    if input.output.trajectories {
        let mut writer = create_file(TRAJECTORIES_FILE)?;
        write_trajectories(&mut writer, transport)?;
        writer.flush()?;
    }
    if input.control_surfaces.is_some() {
        let mut writer = create_file(CONTROL_PLANES_FILE)?;
        write_control_crossings(&mut writer, transport)?;
        writer.flush()?;
    }
    if input.output.velocities {
        let mut writer = create_file(VELOCITIES_FILE)?;
        write_velocities(&mut writer, transport)?;
        writer.flush()?;
    }
    Ok(())
}

fn write_travel_times(
    writer: &mut dyn Write,
    input: &Input,
    transport: &Transport,
) -> std::io::Result<()> {
    let unit = input.time_unit.name();
    if transport.has_matrix_diffusion {
        writeln!(
            writer,
            "# particle, time steps, flux weight, advective time [{unit}], advective + diffusion time [{unit}], diffusion time [{unit}], x, y, z final position [m], beta [{unit}/m], length [m]"
        )?;
    } else {
        writeln!(
            writer,
            "# particle, time steps, flux weight, travel time [{unit}], x, y, z final position [m], beta [{unit}/m], length [m]"
        )?;
    }
    for particle in &transport.particles {
        write!(
            writer,
            "{} {} {:.12e} {:.12e}",
            particle.index + 1,
            particle.time_steps,
            particle.flux_weight,
            particle.travel_time
        )?;
        if let Some(diffusion) = particle.diffusion {
            write!(
                writer,
                " {:.12e} {:.12e}",
                diffusion.total_time, diffusion.diffusion_time
            )?;
        }
        let [x, y, z] = particle.final_position;
        writeln!(
            writer,
            " {x:.12e} {y:.12e} {z:.12e} {:.12e} {:.12e}",
            particle.beta, particle.length
        )?;
    }
    Ok(())
}

/// ASCII VTK PolyData: one polyline for each particle
fn write_trajectories(writer: &mut dyn Write, transport: &Transport) -> std::io::Result<()> {
    let points: Vec<&PathPoint> = transport
        .particles
        .iter()
        .flat_map(|particle| &particle.trajectory)
        .collect();
    let lines: Vec<&ParticleResult> = transport.particles.iter().collect();

    write_header(writer, points.len(), "NumberOfLines", lines.len())?;
    write_points(writer, points.iter().map(|point| point.position))?;
    write_path_point_data(writer, transport, &points)?;

    writeln!(writer, "<CellData>")?;
    write_data_array(
        writer,
        "particle",
        "Int64",
        lines.iter().map(|particle| particle.index + 1),
    )?;
    write_data_array(
        writer,
        "exited",
        "UInt8",
        lines.iter().map(|particle| u8::from(particle.has_exited)),
    )?;
    writeln!(writer, "</CellData>")?;

    writeln!(writer, "<Lines>")?;
    write_data_array(writer, "connectivity", "Int64", 0..points.len())?;
    write_data_array(
        writer,
        "offsets",
        "Int64",
        lines.iter().scan(0, |end, particle| {
            *end += particle.trajectory.len();
            Some(*end)
        }),
    )?;
    writeln!(writer, "</Lines>")?;
    write_footer(writer)
}

/// ASCII VTK PolyData: one vertex for each crossing
fn write_control_crossings(writer: &mut dyn Write, transport: &Transport) -> std::io::Result<()> {
    let crossings: Vec<(usize, _)> = transport
        .particles
        .iter()
        .flat_map(|particle| {
            particle
                .control_crossings
                .iter()
                .map(move |crossing| (particle.index, crossing))
        })
        .collect();
    let points: Vec<&PathPoint> = crossings
        .iter()
        .map(|(_, crossing)| &crossing.point)
        .collect();

    write_header(writer, points.len(), "NumberOfVerts", points.len())?;
    write_points(writer, points.iter().map(|point| point.position))?;
    write_path_point_data(writer, transport, &points)?;
    writeln!(writer, "<CellData>")?;
    write_data_array(
        writer,
        "particle",
        "Int64",
        crossings.iter().map(|(index, _)| index + 1),
    )?;
    write_data_array(
        writer,
        "length",
        "Float64",
        crossings.iter().map(|(_, crossing)| crossing.length),
    )?;
    let diffusion: Vec<DiffusionTimes> = crossings
        .iter()
        .filter_map(|(_, crossing)| crossing.diffusion)
        .collect();
    if !diffusion.is_empty() {
        write_data_array(
            writer,
            "total_time",
            "Float64",
            diffusion.iter().map(|times| times.total_time),
        )?;
        write_data_array(
            writer,
            "diffusion_time",
            "Float64",
            diffusion.iter().map(|times| times.diffusion_time),
        )?;
    }
    writeln!(writer, "</CellData>")?;
    write_vertices(writer, points.len())?;
    write_footer(writer)
}

/// ASCII VTK PolyData: one vertex for each velocity of a node
fn write_velocities(writer: &mut dyn Write, transport: &Transport) -> std::io::Result<()> {
    let network = &transport.network;
    let sectors: Vec<(usize, &Sector)> = network
        .nodes
        .iter()
        .enumerate()
        .flat_map(|(index, node)| node.sectors.iter().map(move |sector| (index, sector)))
        .collect();

    write_header(writer, sectors.len(), "NumberOfVerts", sectors.len())?;
    write_points(
        writer,
        sectors
            .iter()
            .map(|&(node, _)| network.nodes[node].position),
    )?;
    writeln!(writer, "<PointData>")?;
    write_vectors(
        writer,
        "velocity",
        sectors
            .iter()
            .map(|(_, sector)| network.fractures[sector.fracture].vector_to_space(sector.velocity)),
    )?;
    write_data_array(
        writer,
        "node",
        "Int64",
        sectors.iter().map(|&(node, _)| node),
    )?;
    write_data_array(
        writer,
        "fracture_id",
        "Int64",
        sectors
            .iter()
            .map(|(_, sector)| network.fractures[sector.fracture].id),
    )?;
    write_data_array(
        writer,
        "time_step",
        "Float64",
        sectors.iter().map(|(_, sector)| sector.time_step),
    )?;
    writeln!(writer, "</PointData>")?;
    write_vertices(writer, sectors.len())?;
    write_footer(writer)
}

fn write_path_point_data(
    writer: &mut dyn Write,
    transport: &Transport,
    points: &[&PathPoint],
) -> std::io::Result<()> {
    writeln!(writer, "<PointData>")?;
    write_data_array(
        writer,
        "time",
        "Float64",
        points.iter().map(|point| point.time),
    )?;
    write_vectors(
        writer,
        "velocity",
        points.iter().map(|point| point.velocity),
    )?;
    write_data_array(
        writer,
        "fracture_id",
        "Int64",
        points
            .iter()
            .map(|point| transport.network.fractures[point.fracture].id),
    )?;
    write_data_array(
        writer,
        "aperture",
        "Float64",
        points.iter().map(|point| point.aperture),
    )?;
    write_data_array(
        writer,
        "pressure",
        "Float64",
        points.iter().map(|point| point.pressure),
    )?;
    writeln!(writer, "</PointData>")
}

fn write_header(
    writer: &mut dyn Write,
    point_count: usize,
    cell_kind: &str,
    cell_count: usize,
) -> std::io::Result<()> {
    writeln!(writer, r#"<?xml version="1.0"?>"#)?;
    writeln!(writer, r#"<VTKFile type="PolyData" version="1.0">"#)?;
    writeln!(writer, "<PolyData>")?;
    writeln!(
        writer,
        r#"<Piece NumberOfPoints="{point_count}" {cell_kind}="{cell_count}">"#
    )
}

fn write_footer(writer: &mut dyn Write) -> std::io::Result<()> {
    writeln!(writer, "</Piece>")?;
    writeln!(writer, "</PolyData>")?;
    writeln!(writer, "</VTKFile>")
}

fn write_points(
    writer: &mut dyn Write,
    points: impl Iterator<Item = [f64; 3]>,
) -> std::io::Result<()> {
    writeln!(writer, "<Points>")?;
    write_vectors(writer, "Points", points)?;
    writeln!(writer, "</Points>")
}

fn write_vertices(writer: &mut dyn Write, count: usize) -> std::io::Result<()> {
    writeln!(writer, "<Verts>")?;
    write_data_array(writer, "connectivity", "Int64", 0..count)?;
    write_data_array(writer, "offsets", "Int64", 1..=count)?;
    writeln!(writer, "</Verts>")
}

fn write_vectors(
    writer: &mut dyn Write,
    name: &str,
    vectors: impl Iterator<Item = [f64; 3]>,
) -> std::io::Result<()> {
    writeln!(
        writer,
        r#"<DataArray type="Float64" Name="{name}" NumberOfComponents="3" format="ascii">"#
    )?;
    for [x, y, z] in vectors {
        writeln!(writer, "{x} {y} {z}")?;
    }
    writeln!(writer, "</DataArray>")
}

fn write_data_array<T: Display>(
    writer: &mut dyn Write,
    name: &str,
    vtk_type: &str,
    values: impl Iterator<Item = T>,
) -> std::io::Result<()> {
    writeln!(
        writer,
        r#"<DataArray type="{vtk_type}" Name="{name}" format="ascii">"#
    )?;
    for value in values {
        writeln!(writer, "{value}")?;
    }
    writeln!(writer, "</DataArray>")
}
