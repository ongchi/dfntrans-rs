//! The fracture network of the particle tracking: the nodes, triangles and
//! fractures of the mesh with the flow solution. It replaces the data which
//! DFNTrans reads from the LaGriT, FEHM and dfnGen files (`ReadGridInit.c`)
//! and the rotation of the fractures into the xy plane (`RotateFracture.c`).

use std::collections::{BTreeMap, HashMap};

use fehm::Solution;
use fehm::eos;

use crate::DfntransError;

/// The value of π of DFNTrans (`FuncDef.h`); the thresholds of its angle
/// checks use it
#[allow(clippy::approx_constant)]
pub const PI: f64 = 3.141_592_653_59;

/// The mesh with the flow solution
#[derive(Clone, Debug)]
pub struct Network {
    pub nodes: Vec<Node>,
    pub triangles: Vec<Triangle>,
    pub fractures: Vec<Fracture>,
    /// The nodes in the in-flow boundary, in increasing order
    pub inflow_nodes: Vec<usize>,
    /// The nodes in the out-flow boundary, in increasing order
    pub outflow_nodes: Vec<usize>,
}

/// One fracture: the triangles with one value of the `fracture_id` array
#[derive(Clone, Debug)]
pub struct Fracture {
    /// The value of the `fracture_id` cell array
    pub id: i64,
    /// Unit normal of the plane, with z ≥ 0
    pub normal: [f64; 3],
    /// Rotation angle (rad) of the normal onto the z axis (DFNTrans `theta`, in degrees there)
    pub angle: f64,
    /// Rotation of the 3D coordinates into the xy plane (`rot2mat`); `None`
    /// for a fracture parallel to the xy plane
    rotation_to_plane: Option<[[f64; 3]; 3]>,
    /// Rotation back into 3D (`rot3mat`)
    rotation_to_space: Option<[[f64; 3]; 3]>,
    /// The z coordinate after the rotation of the first node of the fracture:
    /// the 3D position of a 2D point is the rotation of (x, y, this value)
    third_coordinate: f64,
    /// The first node of the fracture (`firstnode`): the node with the
    /// smallest index which has this fracture as its first fracture
    pub first_node: usize,
}

impl Fracture {
    /// The axis of the rotation into the xy plane (`nvect_xy` of `poly_info.dat`)
    pub fn rotation_axis(&self) -> [f64; 2] {
        normalized_2d([self.normal[1], -self.normal[0]])
    }

    /// 2D coordinates of a point in the plane of the fracture (`Convertto2d`)
    pub fn to_plane(&self, point: [f64; 3]) -> [f64; 2] {
        match &self.rotation_to_plane {
            Some(rotation) => {
                let rotated = multiply(rotation, point);
                [rotated[0], rotated[1]]
            }
            None => [point[0], point[1]],
        }
    }

    /// The z coordinate of a point after the rotation into the xy plane
    fn rotated_z(&self, point: [f64; 3]) -> f64 {
        match &self.rotation_to_plane {
            Some(rotation) => multiply(rotation, point)[2],
            None => point[2],
        }
    }

    /// 3D position of a 2D point of the fracture (`CalculatePosition3D`)
    pub fn to_space(&self, point: [f64; 2]) -> [f64; 3] {
        match &self.rotation_to_space {
            Some(rotation) => multiply(rotation, [point[0], point[1], self.third_coordinate]),
            None => [point[0], point[1], self.third_coordinate],
        }
    }

    /// 3D vector of a 2D vector of the fracture (`CalculateVelocity3D`)
    pub fn vector_to_space(&self, vector: [f64; 2]) -> [f64; 3] {
        match &self.rotation_to_space {
            Some(rotation) => multiply(rotation, [vector[0], vector[1], 0.0]),
            None => [vector[0], vector[1], 0.0],
        }
    }
}

/// The node types of LaGriT (`itp1`), which DFNTrans uses
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PointType {
    /// 0
    Interior,
    /// 10: on an edge of the mesh (an edge of one triangle)
    Exterior,
    /// 2: on an intersection of fractures
    Interface,
    /// 12: on an intersection and on an edge of the mesh
    ExteriorInterface,
}

/// The in-flow or out-flow boundary of a node
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FlowZone {
    Inflow,
    Outflow,
}

#[derive(Clone, Debug)]
pub struct Node {
    pub position: [f64; 3],
    /// The fractures of the triangles of the node, in increasing order
    pub fractures: Vec<usize>,
    /// The 2D coordinates of the node in each fracture of `fractures`
    pub positions_2d: Vec<[f64; 2]>,
    pub point_type: PointType,
    pub zone: Option<FlowZone>,
    /// The connected nodes, in increasing order
    pub neighbors: Vec<Neighbor>,
    /// MPa
    pub pressure: f64,
    /// kg/m³
    pub density: f64,
    /// Control volume (m³): Voronoi area · aperture
    pub volume: f64,
    /// Largest aperture of the triangles of the first fracture of the node (m)
    pub aperture: f64,
    /// The velocities of the node, which `velocity` calculates: one for each
    /// sector of the triangles around the node
    pub sectors: Vec<Sector>,
}

impl Node {
    /// The 2D coordinates of the node in a fracture; in the first fracture
    /// of the node if it is not in that fracture (as `Xindex`, `Yindex`)
    pub fn position_in(&self, fracture: usize) -> [f64; 2] {
        let slot = self
            .fractures
            .iter()
            .position(|&own| own == fracture)
            .unwrap_or(0);
        self.positions_2d[slot]
    }

    pub fn is_in_fracture(&self, fracture: usize) -> bool {
        self.fractures.contains(&fracture)
    }

    /// The node type number of DFNTrans without the boundary zones: 0, 10, 2 or 12
    pub fn base_type_number(&self) -> u32 {
        match self.point_type {
            PointType::Interior => 0,
            PointType::Exterior => 10,
            PointType::Interface => 2,
            PointType::ExteriorInterface => 12,
        }
    }

    /// The node type number of DFNTrans (`typeN`): the base type plus 300 in
    /// the in-flow boundary or 200 in the out-flow boundary. The port keeps
    /// the conditions of DFNTrans on these numbers.
    pub fn type_number(&self) -> u32 {
        self.base_type_number()
            + match self.zone {
                None => 0,
                Some(FlowZone::Inflow) => 300,
                Some(FlowZone::Outflow) => 200,
            }
    }

    pub fn is_on_interface(&self) -> bool {
        matches!(
            self.point_type,
            PointType::Interface | PointType::ExteriorInterface
        )
    }

    /// Sum of the mass fluxes out of the node into its neighbors (kg/s)
    pub fn total_mass_flux(&self) -> f64 {
        self.neighbors
            .iter()
            .map(|neighbor| neighbor.mass_flux)
            .sum()
    }
}

/// One connection of a node
#[derive(Clone, Debug)]
pub struct Neighbor {
    pub node: usize,
    /// Mass flux from the node to the neighbor (kg/s)
    pub mass_flux: f64,
    /// |Σ ½·cot · aperture| of the connection (m): the Voronoi face area /
    /// the edge length (the `.stor` coefficient of dfnWorks)
    pub area_coefficient: f64,
    /// The triangles with this edge, in increasing order
    pub triangles: Vec<usize>,
}

/// One velocity of a node: the velocity in the triangles of one sector
/// around the node (DFNTrans `velocity[4]`)
#[derive(Clone, Debug)]
pub struct Sector {
    pub fracture: usize,
    /// Darcy velocity / porosity in the 2D coordinates of the fracture (m per time unit)
    pub velocity: [f64; 2],
    /// Time step of a particle at the node (time unit)
    pub time_step: f64,
}

#[derive(Clone, Debug)]
pub struct Triangle {
    pub nodes: [usize; 3],
    pub fracture: usize,
    /// For each corner: the sector of the node (`veloc_ind`)
    pub sectors: [usize; 3],
}

impl Network {
    /// The sector velocity of a corner of a triangle
    pub fn corner_velocity(&self, triangle: usize, corner: usize) -> [f64; 2] {
        let triangle = &self.triangles[triangle];
        self.nodes[triangle.nodes[corner]].sectors[triangle.sectors[corner]].velocity
    }

    /// The sector time step of a corner of a triangle
    pub fn corner_time_step(&self, triangle: usize, corner: usize) -> f64 {
        let triangle = &self.triangles[triangle];
        self.nodes[triangle.nodes[corner]].sectors[triangle.sectors[corner]].time_step
    }

    /// The neighbor entry of `other` in the list of `node`
    pub fn neighbor_index(&self, node: usize, other: usize) -> Option<usize> {
        self.nodes[node]
            .neighbors
            .iter()
            .position(|neighbor| neighbor.node == other)
    }
}

/// Builds the network from the mesh and the flow solution of fehm-rs.
///
/// `inflow` and `outflow` are names of `boundary__<name>` arrays of the
/// mesh. A node in both is in the in-flow boundary (as DFNTrans).
pub fn build_network(
    flow: &fehm::Input,
    solution: &Solution,
    inflow: &[String],
    outflow: &[String],
) -> Result<Network, DfntransError> {
    let mesh = &flow.mesh;
    let node_count = mesh.points.len();
    check_solution_sizes(mesh, solution)?;

    // Fracture indices in increasing order of the fracture_id values
    let fracture_ids: Vec<i64> = {
        let mut ids = mesh.fracture_ids.clone();
        ids.sort_unstable();
        ids.dedup();
        ids
    };
    let fracture_index: HashMap<i64, usize> = fracture_ids
        .iter()
        .enumerate()
        .map(|(index, &id)| (id, index))
        .collect();

    let triangles: Vec<Triangle> = mesh
        .triangles
        .iter()
        .zip(&mesh.fracture_ids)
        .map(|(&nodes, id)| Triangle {
            nodes,
            fracture: fracture_index[id],
            sectors: [0; 3],
        })
        .collect();

    let mut node_fractures: Vec<Vec<usize>> = vec![Vec::new(); node_count];
    for triangle in &triangles {
        for &node in &triangle.nodes {
            node_fractures[node].push(triangle.fracture);
        }
    }
    for fractures in &mut node_fractures {
        fractures.sort_unstable();
        fractures.dedup();
    }
    if let Some(node) = node_fractures.iter().position(Vec::is_empty) {
        return Err(DfntransError::InvalidNetwork(format!(
            "point {node} is not a point of a triangle"
        )));
    }

    // The aperture of a node is the aperture of its first fracture, also
    // on an intersection: DFNTrans reads the apertures (`ReadAperture`)
    // before it finds the second fracture of a node (`AdjacentCells`)
    let mut apertures = vec![0.0_f64; node_count];
    for (triangle, aperture) in triangles.iter().zip(&mesh.apertures) {
        for &node in &triangle.nodes {
            if triangle.fracture == node_fractures[node][0] {
                apertures[node] = apertures[node].max(*aperture);
            }
        }
    }

    // The triangles of each edge, in increasing order
    let mut edge_triangles: BTreeMap<[usize; 2], Vec<usize>> = BTreeMap::new();
    for (index, triangle) in triangles.iter().enumerate() {
        for (first, second) in [(0, 1), (1, 2), (2, 0)] {
            let (a, b) = (triangle.nodes[first], triangle.nodes[second]);
            edge_triangles
                .entry([a.min(b), a.max(b)])
                .or_default()
                .push(index);
        }
    }

    let mut neighbors: Vec<Vec<Neighbor>> = vec![Vec::new(); node_count];
    for (connection, &mass_flux) in solution
        .control_volumes
        .connections
        .iter()
        .zip(&solution.mass_fluxes)
    {
        let [first, second] = connection.nodes;
        let triangles = edge_triangles
            .get(&connection.nodes)
            .cloned()
            .unwrap_or_default();
        let area_coefficient = connection.aperture_coefficient.abs();
        neighbors[first].push(Neighbor {
            node: second,
            mass_flux,
            area_coefficient,
            triangles: triangles.clone(),
        });
        neighbors[second].push(Neighbor {
            node: first,
            mass_flux: -mass_flux,
            area_coefficient,
            triangles,
        });
    }
    for list in &mut neighbors {
        list.sort_by_key(|neighbor| neighbor.node);
    }

    let mut is_exterior = vec![false; node_count];
    for (edge, triangles) in &edge_triangles {
        if triangles.len() == 1 {
            is_exterior[edge[0]] = true;
            is_exterior[edge[1]] = true;
        }
    }

    let zones = flow_zones(mesh, inflow, outflow)?;

    let fractures = build_fractures(&fracture_ids, &triangles, &node_fractures, &mesh.points);

    let nodes: Vec<Node> = (0..node_count)
        .map(|index| {
            let fractures_of_node = std::mem::take(&mut node_fractures[index]);
            let positions_2d = fractures_of_node
                .iter()
                .map(|&fracture| fractures[fracture].to_plane(mesh.points[index]))
                .collect();
            let point_type = match (is_exterior[index], fractures_of_node.len() > 1) {
                (false, false) => PointType::Interior,
                (true, false) => PointType::Exterior,
                (false, true) => PointType::Interface,
                (true, true) => PointType::ExteriorInterface,
            };
            let pressure = solution.pressures[index];
            Node {
                position: mesh.points[index],
                fractures: fractures_of_node,
                positions_2d,
                point_type,
                zone: zones[index],
                neighbors: std::mem::take(&mut neighbors[index]),
                pressure,
                density: eos::density(pressure, flow.temperature).value,
                volume: solution.control_volumes.volumes[index],
                aperture: apertures[index],
                sectors: Vec::new(),
            }
        })
        .collect();

    let nodes_in = |zone: FlowZone| -> Vec<usize> {
        (0..node_count)
            .filter(|&node| zones[node] == Some(zone))
            .collect()
    };
    let inflow_nodes = nodes_in(FlowZone::Inflow);
    let outflow_nodes = nodes_in(FlowZone::Outflow);
    if inflow_nodes.is_empty() || outflow_nodes.is_empty() {
        return Err(DfntransError::InvalidNetwork(
            "the in-flow and the out-flow boundary must each have at least one point (a point in both is in the in-flow boundary)".to_string(),
        ));
    }

    Ok(Network {
        nodes,
        triangles,
        fractures,
        inflow_nodes,
        outflow_nodes,
    })
}

fn check_solution_sizes(mesh: &fehm::mesh::Mesh, solution: &Solution) -> Result<(), DfntransError> {
    let node_count = mesh.points.len();
    if solution.pressures.len() != node_count
        || solution.control_volumes.volumes.len() != node_count
        || solution.mass_fluxes.len() != solution.control_volumes.connections.len()
    {
        return Err(DfntransError::InvalidNetwork(
            "the flow solution is not a solution of this mesh".to_string(),
        ));
    }
    Ok(())
}

/// The boundary zone of each node, from the `boundary__<name>` arrays
fn flow_zones(
    mesh: &fehm::mesh::Mesh,
    inflow: &[String],
    outflow: &[String],
) -> Result<Vec<Option<FlowZone>>, DfntransError> {
    let mut zones = vec![None; mesh.points.len()];
    // Out-flow first: the in-flow zone replaces it
    for (names, zone) in [(outflow, FlowZone::Outflow), (inflow, FlowZone::Inflow)] {
        for name in names {
            let boundary = mesh
                .boundaries
                .iter()
                .find(|boundary| &boundary.name == name)
                .ok_or_else(|| DfntransError::UnknownBoundary(name.clone()))?;
            for &node in &boundary.nodes {
                zones[node] = Some(zone);
            }
        }
    }
    Ok(zones)
}

/// The plane and the rotations of each fracture
fn build_fractures(
    fracture_ids: &[i64],
    triangles: &[Triangle],
    node_fractures: &[Vec<usize>],
    points: &[[f64; 3]],
) -> Vec<Fracture> {
    let normals = fracture_normals(fracture_ids.len(), triangles, points);
    // DFNTrans `firstnode`: the first node of the fracture (nodes sorted by
    // their first fracture). Without such a node (all nodes of the fracture
    // are also in a fracture with a smaller index): a node of its first triangle.
    let first_nodes = (0..fracture_ids.len()).map(|fracture| {
        node_fractures
            .iter()
            .position(|fractures| fractures.first() == Some(&fracture))
            .or_else(|| {
                triangles
                    .iter()
                    .find(|triangle| triangle.fracture == fracture)
                    .map(|triangle| triangle.nodes[0])
            })
            .expect("each fracture has a triangle")
    });

    fracture_ids
        .iter()
        .zip(normals)
        .zip(first_nodes)
        .map(|((&id, normal), first_node)| {
            let horizontal_length = normal[0].hypot(normal[1]);
            let angle = normal[2].clamp(-1.0, 1.0).acos();
            let (rotation_to_plane, rotation_to_space) = if horizontal_length > 1e-12 {
                let axis = [
                    normal[1] / horizontal_length,
                    -normal[0] / horizontal_length,
                ];
                (
                    Some(rotation_matrix(axis, angle)),
                    Some(rotation_matrix([-axis[0], -axis[1]], angle)),
                )
            } else {
                (None, None)
            };
            let mut fracture = Fracture {
                id,
                normal,
                angle,
                rotation_to_plane,
                rotation_to_space,
                third_coordinate: 0.0,
                first_node,
            };
            fracture.third_coordinate = fracture.rotated_z(points[first_node]);
            fracture
        })
        .collect()
}

/// The area-weighted mean of the triangle normals of each fracture, each
/// normal turned to the side of the normal of the largest triangle; z ≥ 0
fn fracture_normals(
    fracture_count: usize,
    triangles: &[Triangle],
    points: &[[f64; 3]],
) -> Vec<[f64; 3]> {
    let area_normals: Vec<[f64; 3]> = triangles
        .iter()
        .map(|triangle| {
            let [a, b, c] = triangle.nodes.map(|node| points[node]);
            cross(subtract(b, a), subtract(c, a))
        })
        .collect();

    let mut reference: Vec<Option<[f64; 3]>> = vec![None; fracture_count];
    for (triangle, normal) in triangles.iter().zip(&area_normals) {
        let largest = &mut reference[triangle.fracture];
        if largest.is_none_or(|largest| norm(*normal) > norm(largest)) {
            *largest = Some(*normal);
        }
    }

    let mut sums = vec![[0.0; 3]; fracture_count];
    for (triangle, normal) in triangles.iter().zip(&area_normals) {
        let sign = match reference[triangle.fracture] {
            Some(reference) if dot(reference, *normal) < 0.0 => -1.0,
            _ => 1.0,
        };
        for axis in 0..3 {
            sums[triangle.fracture][axis] += sign * normal[axis];
        }
    }

    sums.into_iter()
        .map(|sum| {
            let length = norm(sum);
            let unit = sum.map(|component| component / length);
            if unit[2] < 0.0 {
                unit.map(|component| -component)
            } else {
                unit
            }
        })
        .collect()
}

/// Rotation about the horizontal unit axis (x, y, 0) by `angle` (Rodrigues;
/// the matrices of `Convertto2d` and `Convertto3d`)
fn rotation_matrix(axis: [f64; 2], angle: f64) -> [[f64; 3]; 3] {
    let (sine, cosine) = angle.sin_cos();
    let [x, y] = axis;
    [
        [
            x * x * (1.0 - cosine) + cosine,
            x * y * (1.0 - cosine),
            y * sine,
        ],
        [
            x * y * (1.0 - cosine),
            y * y * (1.0 - cosine) + cosine,
            -x * sine,
        ],
        [-y * sine, x * sine, cosine],
    ]
}

fn multiply(matrix: &[[f64; 3]; 3], vector: [f64; 3]) -> [f64; 3] {
    matrix.map(|row| row[0] * vector[0] + row[1] * vector[1] + row[2] * vector[2])
}

fn normalized_2d(vector: [f64; 2]) -> [f64; 2] {
    let length = vector[0].hypot(vector[1]);
    if length > 0.0 {
        vector.map(|component| component / length)
    } else {
        [0.0, 0.0]
    }
}

pub(crate) fn subtract(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(crate) fn norm(vector: [f64; 3]) -> f64 {
    dot(vector, vector).sqrt()
}
