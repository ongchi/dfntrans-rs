//! The work flow dfngen-rs → fehm-rs → dfntrans-rs on an irregular mesh.
//! The network comes from the case `mesh` of dfngen-rs with 300
//! fractures, only the largest cluster and no triple intersections (88
//! fractures, about 31 000 nodes); the flow goes from the face x = -5 to the
//! face x = 5 of the box domain. Slow without `--release`.

use std::cell::RefCell;
use std::io::Write;
use std::rc::Rc;

use dfntrans::{build_input_from_text, track_particles};

const DFNGEN_DATA: &str = "../dfngen-rs/tests/data";

/// The text of `dfn_mesh.vtu` of the network
fn dfngen_mesh() -> String {
    let case = std::fs::read_to_string(format!("{DFNGEN_DATA}/mesh/mesh.toml")).unwrap();
    let case = [
        ("seed = 12345", "seed = 3"),
        ("n_poly = 100", "n_poly = 300"),
        (
            "triple_intersections = true",
            "triple_intersections = false",
        ),
        ("keep_only_largest = false", "keep_only_largest = true"),
    ]
    .iter()
    .fold(case, |text, (from, to)| {
        assert!(text.contains(from), "{from}");
        text.replace(from, to)
    });
    let domain = std::fs::read_to_string(format!("{DFNGEN_DATA}/box.vtu")).unwrap();
    let (input, families) = dfngen::build_input_from_text(&case, &domain, None).unwrap();
    let network = dfngen::generate(input, families, &mut |_| {}).unwrap();

    let mesh: Rc<RefCell<Vec<u8>>> = Rc::default();
    dfngen::write_output(&network, &mut |name| {
        if name == "dfn_mesh.vtu" {
            Ok(Box::new(MemoryFile(Rc::clone(&mesh))))
        } else {
            Ok(Box::new(std::io::sink()))
        }
    })
    .unwrap();
    String::from_utf8(mesh.borrow().clone()).unwrap()
}

/// A writer into a shared buffer
struct MemoryFile(Rc<RefCell<Vec<u8>>>);

impl Write for MemoryFile {
    fn write(&mut self, buffer: &[u8]) -> std::io::Result<usize> {
        self.0.borrow_mut().extend_from_slice(buffer);
        Ok(buffer.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// The flow input: a fixed pressure on the faces x = -5 (`west`) and x = 5
/// (`east`), the points of `boundary__mesh_surface` on them
fn flow_input(mesh_text: &str) -> fehm::Input {
    let flow_text = "mesh = \"dfn_mesh.vtu\"\ntemperature = 20.0\ngravity = false\n\n\
                     [[boundary]]\nname = \"mesh_surface\"\npressure = 1.0\n";
    let mut flow = fehm::build_input_from_text(flow_text, mesh_text).unwrap();
    let surface = flow
        .mesh
        .boundaries
        .iter()
        .find(|boundary| boundary.name == "mesh_surface")
        .unwrap()
        .nodes
        .clone();
    let face = |x: f64| -> Vec<usize> {
        surface
            .iter()
            .copied()
            .filter(|&node| (flow.mesh.points[node][0] - x).abs() < 1e-6)
            .collect()
    };
    let (west, east) = (face(-5.0), face(5.0));
    assert!(!west.is_empty() && !east.is_empty());
    flow.boundaries = vec![
        fehm::io::input::PressureBoundary {
            name: "west".to_string(),
            pressure: 1.1,
            nodes: west.clone(),
        },
        fehm::io::input::PressureBoundary {
            name: "east".to_string(),
            pressure: 1.0,
            nodes: east.clone(),
        },
    ];
    for (name, nodes) in [("west", west), ("east", east)] {
        flow.mesh.boundaries.push(fehm::mesh::BoundaryNodes {
            name: name.to_string(),
            nodes,
        });
    }
    flow
}

#[test]
fn particles_cross_a_dfngen_network() {
    let flow = flow_input(&dfngen_mesh());
    let solution = fehm::solve(&flow).unwrap();
    let input = build_input_from_text(
        "flow = \"flow.toml\"\nseed = 12345\ntime_unit = \"days\"\nporosity = 1.0\n\
         max_time_steps = 1000000\ninflow = [\"west\"]\noutflow = [\"east\"]\n\
         routing = \"complete_mixing\"\n\n[seeding]\nmethod = \"flux_weighted\"\nparticle_count = 500\n",
    )
    .unwrap();
    let transport = track_particles(&input, &flow, &solution).unwrap();

    assert!(transport.network.fractures.len() > 10);
    // Many nodes on intersections, thus many crossings of fractures
    let interface_nodes = transport
        .network
        .nodes
        .iter()
        .filter(|node| node.fractures.len() > 1)
        .count();
    assert!(interface_nodes > 1000, "{interface_nodes}");
    assert!(transport.particle_count >= 500);
    assert_eq!(transport.exit_count, transport.particle_count);
    // The in-flow and the out-flow volume flux agree (mass balance of fehm-rs)
    let balance = (transport.inflow_volume_flux - transport.outflow_volume_flux).abs()
        / transport.inflow_volume_flux;
    assert!(balance < 1e-3, "{balance}");
    for particle in &transport.particles {
        // On the face, or in the last cell if no out-flow edge is found (as
        // DFNTrans `FinalPosition`); the cells there have a size of 0.1 m
        assert!((particle.final_position[0] - 5.0).abs() < 0.1);
        assert!(particle.travel_time > 0.0 && particle.length > 10.0 * 0.99);
    }
}
