//! Reconstruction of the Darcy velocity at the nodes from the fluxes of the
//! control volume faces (DFNTrans `VelocityReconstruction.c`; Painter et
//! al. 2012, Eq. 5 and 7), and the time step of each node (`DefineTimeStep`).
//!
//! DFNTrans has one velocity for a node in one fracture, and two velocities
//! in each fracture for a node on an intersection: the intersection line
//! divides its control volume into two halves (`HalfPolygonVelocity`). The
//! port divides the triangles around a node into sectors in the same way. A
//! node on more than two fractures (a crossing of intersection lines, which
//! DFNTrans does not support) gets one sector for each group of triangles
//! between two intersection edges.

use tracing::warn;

use crate::network::{Network, Node, PI, Sector};

/// Calculates the sector velocities and time steps of all nodes and the
/// sector of each triangle corner, then corrects the velocities at the
/// fracture edges (`DarcyVelocity`, `BoundaryCells`, `DefineTimeStep`).
///
/// The velocity is the Darcy velocity / (density · porosity), in m per time
/// unit; `seconds_per_time_unit` converts it.
pub fn reconstruct_velocities(network: &mut Network, porosity: f64, seconds_per_time_unit: f64) {
    let mut nodes_without_boundary_edges = 0;
    let mut singular = SingularSectors::default();
    for node in 0..network.nodes.len() {
        let has_boundary_edges = if network.nodes[node].fractures.len() == 1 {
            single_fracture_velocity(network, node, porosity, seconds_per_time_unit)
        } else {
            // Each fracture of the node, also after one without boundary edges
            let mut has_all = true;
            for slot in 0..network.nodes[node].fractures.len() {
                has_all &= interface_velocities(
                    network,
                    node,
                    slot,
                    porosity,
                    seconds_per_time_unit,
                    &mut singular,
                );
            }
            has_all
        };
        if !has_boundary_edges {
            nodes_without_boundary_edges += 1;
        }
    }
    if nodes_without_boundary_edges > 0 {
        warn!(
            "{nodes_without_boundary_edges} nodes on a fracture edge do not have two edges on it; their velocity has no boundary condition"
        );
    }
    if singular.sector_count > 0 {
        warn!(
            "{} sectors of nodes on an intersection have faces in one direction; they get the velocity of all faces of the node in their fracture ({} of them: 0, these faces also have one direction)",
            singular.sector_count, singular.without_velocity_count
        );
    }

    correct_boundary_velocities(network);
    define_time_steps(network);
}

/// One row of the least squares fit: the face normal times the face area
/// (`normxarea`), and the mass flux through the face
#[derive(Clone, Copy, Debug)]
struct FaceFlux {
    normal_area: [f64; 2],
    mass_flux: f64,
}

/// The counts of the sectors whose faces have one direction (`is_singular`)
#[derive(Default)]
struct SingularSectors {
    sector_count: usize,
    /// Also all faces of the node in the fracture have one direction: velocity 0
    without_velocity_count: usize,
}

/// The no-flow condition on the fracture edge at an exterior node (`struct lb`)
#[derive(Clone, Copy, Debug)]
struct BoundaryCondition {
    length: f64,
    normal: [f64; 2],
}

/// The fit of the velocity of a sector
#[derive(Clone, Copy, Debug)]
enum Fit {
    /// `VelocityInteriorNode`
    Interior,
    /// `VelocityExteriorNode`
    Exterior(BoundaryCondition),
    /// DFNTrans leaves the velocity 0
    None,
}

/// The face flux of a neighbor of a node, in the 2D coordinates of `fracture`
fn face_flux(network: &Network, node: usize, neighbor_index: usize, fracture: usize) -> FaceFlux {
    let own = network.nodes[node].position_in(fracture);
    let neighbor = &network.nodes[node].neighbors[neighbor_index];
    let other = network.nodes[neighbor.node].position_in(fracture);
    FaceFlux {
        normal_area: [
            -(own[0] - other[0]) * neighbor.area_coefficient,
            -(own[1] - other[1]) * neighbor.area_coefficient,
        ],
        mass_flux: neighbor.mass_flux,
    }
}

/// Velocity of a node in one fracture: one sector with all its faces
/// (the interior and exterior parts of `DarcyVelocity`). Returns false if
/// the node is on a fracture edge without two edges on it.
fn single_fracture_velocity(
    network: &mut Network,
    node: usize,
    porosity: f64,
    seconds_per_time_unit: f64,
) -> bool {
    let fracture = network.nodes[node].fractures[0];
    let faces: Vec<FaceFlux> = (0..network.nodes[node].neighbors.len())
        .map(|index| face_flux(network, node, index, fracture))
        .collect();

    let mut has_boundary_edges = true;
    // Exterior nodes in the in-flow or out-flow zone are fitted as interior nodes
    let fit = if network.nodes[node].type_number() == 10 {
        let boundary_edges: Vec<usize> = network.nodes[node]
            .neighbors
            .iter()
            .enumerate()
            .filter(|(_, neighbor)| neighbor.triangles.len() == 1)
            .map(|(index, _)| index)
            .collect();
        match boundary_edges.as_slice() {
            // `edge01` is the first and `edge02` the last of them
            [first, .., last] => {
                Fit::Exterior(boundary_condition(network, node, *first, *last, fracture))
            }
            _ => {
                has_boundary_edges = false;
                Fit::Interior
            }
        }
    } else {
        Fit::Interior
    };

    let velocity = fit_velocity(
        &network.nodes[node],
        &faces,
        fit,
        porosity,
        seconds_per_time_unit,
    );
    network.nodes[node].sectors = vec![Sector {
        fracture,
        velocity,
        time_step: 0.0,
    }];
    has_boundary_edges
}

/// The division of the faces of a node in one fracture into sectors
struct Sectoring {
    /// The faces (neighbor indices) of each sector
    faces: Vec<Vec<usize>>,
    /// (triangle, corner, sector): the sector of each corner; a later entry
    /// replaces an earlier one
    corners: Vec<(usize, usize, usize)>,
}

/// Velocities of a node on an intersection in the fracture of `slot`
/// (`HalfPolygonVelocity`). Returns false if the node is on a fracture edge
/// and the fit has no boundary condition. Unlike DFNTrans, a sector whose
/// faces have one direction gets the fit of all faces of the node in the
/// fracture (DFNTrans: a determinant of 1e-75, a velocity of about 1e40).
fn interface_velocities(
    network: &mut Network,
    node: usize,
    slot: usize,
    porosity: f64,
    seconds_per_time_unit: f64,
    singular: &mut SingularSectors,
) -> bool {
    let fracture = network.nodes[node].fractures[slot];
    let triangles_in_fracture = |edge: usize| {
        network.nodes[node].neighbors[edge]
            .triangles
            .iter()
            .filter(|&&triangle| network.triangles[triangle].fracture == fracture)
            .count()
    };

    // The faces of the node in this fracture (`fract_j1`, `fract_j2`)
    let edges: Vec<usize> = (0..network.nodes[node].neighbors.len())
        .filter(|&edge| triangles_in_fracture(edge) > 0)
        .collect();
    // Faces with one triangle in this fracture: on the fracture edge (`ss == 1`)
    let boundary_edges: Vec<usize> = edges
        .iter()
        .copied()
        .filter(|&edge| triangles_in_fracture(edge) == 1)
        .collect();
    // Unlike DFNTrans (a face to a node on an intersection), a face with
    // triangles in two fractures: the neighbor of a node on the fracture
    // edge can be on another intersection
    let intersection_edges: Vec<usize> = edges
        .iter()
        .copied()
        .filter(|&edge| is_intersection_edge(network, node, edge))
        .collect();

    let sectoring = if network.nodes[node].fractures.len() == 2 {
        half_planes(network, node, fracture, &edges, &intersection_edges)
    } else {
        triangle_groups(network, node, fracture, &edges)
    };

    let mut has_boundary_edges = true;
    let fit = if boundary_edges.is_empty() || network.nodes[node].zone.is_some() {
        Fit::Interior
    } else if let [first, second] = boundary_edges[..] {
        Fit::Exterior(boundary_condition(network, node, first, second, fracture))
    } else {
        has_boundary_edges = false;
        Fit::None
    };

    // The sectors of this fracture follow those of the fractures before it
    if slot == 0 {
        network.nodes[node].sectors.clear();
    }
    let first_sector = network.nodes[node].sectors.len();

    let all_face_fluxes: Vec<FaceFlux> = edges
        .iter()
        .map(|&edge| face_flux(network, node, edge, fracture))
        .collect();
    for faces in &sectoring.faces {
        let face_fluxes: Vec<FaceFlux> = faces
            .iter()
            .map(|&edge| face_flux(network, node, edge, fracture))
            .collect();
        let fitted_faces = if face_fluxes.is_empty() {
            None
        } else if !is_singular(&face_fluxes) {
            Some(&face_fluxes)
        } else {
            singular.sector_count += 1;
            if is_singular(&all_face_fluxes) {
                singular.without_velocity_count += 1;
                None
            } else {
                Some(&all_face_fluxes)
            }
        };
        let velocity = fitted_faces.map_or([0.0, 0.0], |faces| {
            fit_velocity(
                &network.nodes[node],
                faces,
                fit,
                porosity,
                seconds_per_time_unit,
            )
        });
        network.nodes[node].sectors.push(Sector {
            fracture,
            velocity,
            time_step: 0.0,
        });
    }

    // Default sector of a corner which no face gives (DFNTrans: velocity 0
    // of the first fracture; here the first sector of this fracture)
    for &edge in &edges {
        for triangle in network.nodes[node].neighbors[edge].triangles.clone() {
            if network.triangles[triangle].fracture == fracture
                && let Some(corner) = corner_of(network, triangle, node)
            {
                network.triangles[triangle].sectors[corner] = first_sector;
            }
        }
    }
    for (triangle, corner, sector) in sectoring.corners {
        network.triangles[triangle].sectors[corner] = first_sector + sector;
    }
    has_boundary_edges
}

/// The sine of the largest angle between a face and the intersection line
/// for which the face is on the line
const LINE_TOLERANCE: f64 = 1e-10;

/// The relative difference below which two edge lengths are equal
const LENGTH_TOLERANCE: f64 = 1e-10;

/// A face whose triangles are in more than one fracture: an edge of an
/// intersection line
fn is_intersection_edge(network: &Network, node: usize, edge: usize) -> bool {
    let triangles = &network.nodes[node].neighbors[edge].triangles;
    triangles.iter().any(|&triangle| {
        network.triangles[triangle].fracture != network.triangles[triangles[0]].fracture
    })
}

/// The corner of `node` in a triangle
fn corner_of(network: &Network, triangle: usize, node: usize) -> Option<usize> {
    network.triangles[triangle]
        .nodes
        .iter()
        .position(|&corner| corner == node)
}

/// The two halves of the control volume on the two sides of the line of
/// the last intersection edge (`HalfPolygonVelocity`). A face on the line
/// is in both halves. The corners of the triangles of a face to a node not
/// on an intersection get the sector of the face, in the order of the faces.
fn half_planes(
    network: &Network,
    node: usize,
    fracture: usize,
    edges: &[usize],
    intersection_edges: &[usize],
) -> Sectoring {
    let own = network.nodes[node].position_in(fracture);
    let line = intersection_edges.last().map_or([0.0, 0.0], |&edge| {
        let other = network.nodes[network.nodes[node].neighbors[edge].node].position_in(fracture);
        [own[0] - other[0], own[1] - other[1]]
    });

    let mut sectoring = Sectoring {
        faces: vec![Vec::new(), Vec::new()],
        corners: Vec::new(),
    };
    for &edge in edges {
        let neighbor = &network.nodes[node].neighbors[edge];
        let other = network.nodes[neighbor.node].position_in(fracture);
        let offset = [own[0] - other[0], own[1] - other[1]];
        // A face on the line is in both halves. Unlike DFNTrans, with a
        // tolerance: the other edge of a straight intersection line is on
        // the line, and the rounding must not put it into one half only.
        let mut product = offset[0] * line[1] - offset[1] * line[0];
        if product.abs() <= LINE_TOLERANCE * offset[0].hypot(offset[1]) * line[0].hypot(line[1]) {
            product = 0.0;
        }
        let assigns_corners = !network.nodes[neighbor.node].is_on_interface();

        for (sector, is_on_side) in [product <= 0.0, product >= 0.0].into_iter().enumerate() {
            if !is_on_side {
                continue;
            }
            sectoring.faces[sector].push(edge);
            if assigns_corners {
                for &triangle in &neighbor.triangles {
                    if let Some(corner) = corner_of(network, triangle, node) {
                        sectoring.corners.push((triangle, corner, sector));
                    }
                }
            }
        }
    }
    sectoring
}

/// Groups of the triangles around the node in one fracture which are
/// connected through faces that are not intersection edges (edges of
/// triangles of more than one fracture). For a node on more than two
/// fractures.
fn triangle_groups(network: &Network, node: usize, fracture: usize, edges: &[usize]) -> Sectoring {
    let neighbors = &network.nodes[node].neighbors;
    let in_fracture = |edge: usize| -> Vec<usize> {
        neighbors[edge]
            .triangles
            .iter()
            .copied()
            .filter(|&triangle| network.triangles[triangle].fracture == fracture)
            .collect()
    };

    // The group of each triangle: start with one group for each triangle,
    // then join the groups of the triangles of each face which is not an
    // intersection edge
    let mut triangles: Vec<usize> = edges.iter().flat_map(|&edge| in_fracture(edge)).collect();
    triangles.sort_unstable();
    triangles.dedup();
    let mut group_of: Vec<usize> = (0..triangles.len()).collect();
    let index_of = |triangle: usize| {
        triangles
            .binary_search(&triangle)
            .expect("triangle of a face")
    };
    for &edge in edges {
        if is_intersection_edge(network, node, edge) {
            continue;
        }
        let members: Vec<usize> = in_fracture(edge).into_iter().map(index_of).collect();
        if let [first, second] = members[..] {
            let (kept, replaced) = (group_of[first], group_of[second]);
            for group in &mut group_of {
                if *group == replaced {
                    *group = kept;
                }
            }
        }
    }

    // Number the groups from 0 in the order of their first triangle
    let mut group_numbers: Vec<usize> = Vec::new();
    let sector_of = |group: usize, numbers: &mut Vec<usize>| {
        numbers
            .iter()
            .position(|&known| known == group)
            .unwrap_or_else(|| {
                numbers.push(group);
                numbers.len() - 1
            })
    };
    let mut sectoring = Sectoring {
        faces: Vec::new(),
        corners: Vec::new(),
    };
    for (index, &triangle) in triangles.iter().enumerate() {
        let sector = sector_of(group_of[index], &mut group_numbers);
        if let Some(corner) = corner_of(network, triangle, node) {
            sectoring.corners.push((triangle, corner, sector));
        }
    }
    sectoring.faces = vec![Vec::new(); group_numbers.len()];
    for &edge in edges {
        for triangle in in_fracture(edge) {
            let sector = sector_of(group_of[index_of(triangle)], &mut group_numbers);
            if !sectoring.faces[sector].contains(&edge) {
                sectoring.faces[sector].push(edge);
            }
        }
    }
    sectoring
}

/// The no-flow condition from two edges on the fracture edge (`DefineBoundaryAngle`)
fn boundary_condition(
    network: &Network,
    node: usize,
    first_edge: usize,
    second_edge: usize,
    fracture: usize,
) -> BoundaryCondition {
    let own = network.nodes[node].position_in(fracture);
    let edge_vector = |edge: usize| {
        let other = network.nodes[network.nodes[node].neighbors[edge].node].position_in(fracture);
        [other[0] - own[0], other[1] - own[1]]
    };
    let u = edge_vector(first_edge);
    let v = edge_vector(second_edge);
    let length_u = u[0].hypot(u[1]);
    let length_v = v[0].hypot(v[1]);
    // The normal of the longer edge. Unlike DFNTrans, edges of equal length
    // within the rounding take the second edge (as DFNTrans in exact
    // arithmetic); a mesh has many edges of equal length on its outline.
    let normal = if length_u > length_v * (1.0 + LENGTH_TOLERANCE) {
        [-u[1], u[0]]
    } else {
        [-v[1], v[0]]
    };
    BoundaryCondition {
        length: length_u + length_v,
        normal,
    }
}

/// GᵀG of the faces (`MatrixProducts`)
fn normal_products(faces: &[FaceFlux]) -> [[f64; 2]; 2] {
    let mut product = [[0.0; 2]; 2];
    for face in faces {
        for (row, &first) in product.iter_mut().zip(&face.normal_area) {
            for (entry, &second) in row.iter_mut().zip(&face.normal_area) {
                *entry += first * second;
            }
        }
    }
    product
}

/// The relative size of the determinant of GᵀG (to the product of its
/// diagonal) at or below which the faces have one direction
const SINGULAR_TOLERANCE: f64 = 1e-10;

/// The faces have one direction within the rounding, as two faces on a
/// straight fracture edge: GᵀG is singular and the fit has no solution
fn is_singular(faces: &[FaceFlux]) -> bool {
    let product = normal_products(faces);
    let determinant = product[0][0] * product[1][1] - product[0][1] * product[1][0];
    determinant.abs() <= SINGULAR_TOLERANCE * product[0][0] * product[1][1]
}

/// Least squares fit of the velocity to the face fluxes, divided by the
/// density and the porosity (`MatrixProducts`, `VelocityInteriorNode`,
/// `VelocityExteriorNode`)
fn fit_velocity(
    node: &Node,
    faces: &[FaceFlux],
    fit: Fit,
    porosity: f64,
    seconds_per_time_unit: f64,
) -> [f64; 2] {
    const SMALLEST_DETERMINANT: f64 = 1e-75;

    let product = normal_products(faces);
    let mut determinant = product[0][0] * product[1][1] - product[0][1] * product[1][0];
    if determinant.abs() < SMALLEST_DETERMINANT {
        determinant = SMALLEST_DETERMINANT;
    }
    let inverse = [
        [product[1][1] / determinant, -product[0][1] / determinant],
        [-product[1][0] / determinant, product[0][0] / determinant],
    ];

    // (GᵀG)⁻¹ Gᵀ F
    let mut flux = [0.0; 2];
    for (m, component) in flux.iter_mut().enumerate() {
        for face in faces {
            let weight = inverse[m][0] * face.normal_area[0] + inverse[m][1] * face.normal_area[1];
            *component += weight * face.mass_flux;
        }
    }

    let flux = match fit {
        Fit::Interior => flux,
        Fit::Exterior(condition) => {
            let b = condition
                .normal
                .map(|component| component * condition.length);
            let b_flux = b[0] * flux[0] + b[1] * flux[1];
            let inverse_b = [
                b[0] * inverse[0][0] + b[1] * inverse[1][0],
                b[0] * inverse[0][1] + b[1] * inverse[1][1],
            ];
            let scale = 1.0 / (inverse_b[0] * b[0] + inverse_b[1] * b[1]);
            let scaled_b = b.map(|component| component * scale);
            let correction = [
                inverse[0][0] * scaled_b[0] + inverse[0][1] * scaled_b[1],
                inverse[1][0] * scaled_b[0] + inverse[1][1] * scaled_b[1],
            ];
            [
                flux[0] - correction[0] * b_flux,
                flux[1] - correction[1] * b_flux,
            ]
        }
        Fit::None => return [0.0, 0.0],
    };
    flux.map(|component| component / (node.density * porosity) * seconds_per_time_unit)
}

/// Angle between two 2D vectors (`DefineAngle`; single precision `acosf` as DFNTrans)
pub fn angle_between(u: [f64; 2], v: [f64; 2]) -> f64 {
    let length_u = (u[0] * u[0] + u[1] * u[1]).sqrt();
    let length_v = (v[0] * v[0] + v[1] * v[1]).sqrt();
    let cosine = (u[0] * v[0] + u[1] * v[1]) / (length_u * length_v);
    f64::from((cosine as f32).acos())
}

/// Turns a velocity by `angle`, first counterclockwise, else clockwise, so
/// that it is parallel to `direction`; `None` if neither is (the rotation
/// steps of `CornerVelocity` and `BoundaryCells`)
fn turn_along(
    velocity: [f64; 2],
    angle: f64,
    direction: [f64; 2],
    tolerance: f64,
) -> Option<[f64; 2]> {
    let (sine, cosine) = angle.sin_cos();
    let mut turned = [
        velocity[0] * cosine - velocity[1] * sine,
        velocity[0] * sine + velocity[1] * cosine,
    ];
    let mut remaining = angle_between(turned, direction);
    if remaining > tolerance {
        turned = [
            velocity[0] * cosine + velocity[1] * sine,
            -velocity[0] * sine + velocity[1] * cosine,
        ];
        remaining = angle_between(turned, direction);
    }
    (remaining < tolerance).then_some(turned)
}

fn sector_velocity(network: &Network, node: usize, sector: usize) -> [f64; 2] {
    network.nodes[node].sectors[sector].velocity
}

fn set_sector_velocity(network: &mut Network, node: usize, sector: usize, velocity: [f64; 2]) {
    network.nodes[node].sectors[sector].velocity = velocity;
}

/// Corrections of the velocities at the fracture edges (`BoundaryCells`):
/// at the corners of a fracture, a velocity which points out of the
/// fracture is turned along its edge; of two velocities on one edge which
/// point against each other, one is changed, else a particle stays between
/// them.
fn correct_boundary_velocities(network: &mut Network) {
    for cell in 0..network.triangles.len() {
        let triangle = &network.triangles[cell];
        let cell_fracture = triangle.fracture;

        // k1, k2: the first two nodes of the cell on the fracture edge; kk3: a third one
        let mut on_edge = triangle
            .nodes
            .into_iter()
            .zip(triangle.sectors)
            .filter(|&(node, _)| matches!(network.nodes[node].type_number(), 10 | 12 | 310 | 312));
        let (Some((k1, v1)), Some((k2, v2))) = (on_edge.next(), on_edge.next()) else {
            continue;
        };
        let third = on_edge.next();

        // The node k3 on the fracture edge next to the center node, and its
        // sector in the cell `celln`. The search around k2 keeps k3 and
        // `celln` of the search around k1 if it finds no node (as DFNTrans).
        let mut k3: Option<usize> = None;
        let mut v3 = 0;
        let mut corner_cell = 0;
        let mut is_straight_at_k1 = false;
        let mut is_straight_at_k2 = false;

        for (center, other) in [(k1, k2), (k2, k1)] {
            let mut third_index = 0;
            let mut other_index = 0;
            for (index, neighbor) in network.nodes[center].neighbors.iter().enumerate() {
                let neighbor_node = &network.nodes[neighbor.node];
                if (matches!(neighbor_node.base_type_number(), 10 | 12)
                    || neighbor_node.type_number() == 310)
                    && neighbor.node != other
                    && neighbor.triangles.len() == 1
                    && network.triangles[neighbor.triangles[0]].fracture == cell_fracture
                {
                    k3 = Some(neighbor.node);
                    third_index = index;
                    corner_cell = neighbor.triangles[0];
                }
                if neighbor.node == other {
                    other_index = index;
                }
            }
            let has_one_triangle = |index: usize| {
                network.nodes[center]
                    .neighbors
                    .get(index)
                    .is_none_or(|neighbor| neighbor.triangles.len() <= 1)
            };
            let Some(k3_node) = k3 else {
                continue;
            };
            if !(has_one_triangle(other_index) && has_one_triangle(third_index)) {
                continue;
            }
            if let Some(corner) = corner_of(network, corner_cell, k3_node) {
                v3 = network.triangles[corner_cell].sectors[corner];
            }
            if center == k1 {
                is_straight_at_k1 = corner_velocity(network, cell, [k1, k2, k3_node], [v1, v2, v3]);
            } else {
                is_straight_at_k2 = corner_velocity(network, cell, [k2, k1, k3_node], [v2, v1, v3]);
            }
        }

        if is_straight_at_k1 || is_straight_at_k2 {
            correct_opposite_velocities(network, cell, (k1, v1), (k2, v2));
        }

        if let Some((kk3, vv3)) = third
            && network.nodes[kk3].neighbors.len() == 2
        {
            corner_velocity(network, cell, [kk3, k2, k1], [vv3, v2, v1]);
        }
    }
}

/// Two velocities on one fracture edge which point against each other: on
/// an edge without intersection, the velocity of k1 is turned to its
/// neighbor with the lowest pressure; else the velocity of the node on the
/// intersection is reversed (part of `BoundaryCells`)
fn correct_opposite_velocities(
    network: &mut Network,
    cell: usize,
    (k1, v1): (usize, usize),
    (k2, v2): (usize, usize),
) {
    const TOLERANCE: f64 = 0.0001;
    let velocity_1 = sector_velocity(network, k1, v1);
    let velocity_2 = sector_velocity(network, k2, v2);
    let are_opposite = velocity_1[0] * velocity_2[0] < 0.0 && velocity_1[1] * velocity_2[1] < 0.0;
    let type_sum = network.nodes[k1].type_number() + network.nodes[k2].type_number();
    let cell_fracture = network.triangles[cell].fracture;

    let angle = angle_between(velocity_1, velocity_2);
    if are_opposite && angle < PI + TOLERANCE && angle > PI - TOLERANCE && type_sum == 20 {
        let mut lowest_pressure = 1000.0;
        let mut direction = [0.0, 0.0];
        let own = network.nodes[k1].positions_2d[0];
        for neighbor in &network.nodes[k1].neighbors {
            let neighbor_node = &network.nodes[neighbor.node];
            if neighbor_node.pressure < lowest_pressure && neighbor.node != k2 {
                lowest_pressure = neighbor_node.pressure;
                let other = neighbor_node.position_in(cell_fracture);
                direction = [other[0] - own[0], other[1] - own[1]];
            }
        }
        let turn = angle_between(velocity_1, direction);
        if let Some(turned) = turn_along(velocity_1, turn, direction, TOLERANCE) {
            set_sector_velocity(network, k1, v1, turned);
        }
    }

    if are_opposite && type_sum != 20 {
        if network.nodes[k1].type_number() == 12 {
            set_sector_velocity(network, k1, v1, velocity_1.map(|component| -component));
        } else if network.nodes[k2].type_number() == 12 {
            set_sector_velocity(network, k2, v2, velocity_2.map(|component| -component));
        }
    }
}

/// At a corner of a fracture (two edges on the fracture edge with an angle
/// between them): a velocity of the corner node m1 which points along one
/// edge, away from the other, is turned along the other edge; a velocity
/// at an in-flow node is turned along the edge to a node which is not in
/// the in-flow zone (`CornerVelocity`). Returns true if the two edges are
/// on one line (no corner).
fn corner_velocity(
    network: &mut Network,
    cell: usize,
    [m1, m2, m3]: [usize; 3],
    [s1, s2, s3]: [usize; 3],
) -> bool {
    const TOLERANCE: f64 = 0.001;
    let fracture = network.triangles[cell].fracture;
    let position = |node: usize| network.nodes[node].position_in(fracture);
    let center = position(m1);
    let u = [position(m2)[0] - center[0], position(m2)[1] - center[1]];
    let v = [position(m3)[0] - center[0], position(m3)[1] - center[1]];
    let edge_angle = angle_between(u, v);

    if !(edge_angle < PI - TOLERANCE && edge_angle > TOLERANCE) {
        return (PI - TOLERANCE..=PI + TOLERANCE).contains(&edge_angle)
            || (-TOLERANCE..=TOLERANCE).contains(&edge_angle);
    }

    let velocity_1 = sector_velocity(network, m1, s1);
    let velocity_2 = sector_velocity(network, m2, s2);
    let velocity_3 = sector_velocity(network, m3, s3);
    let angle_v = angle_between(velocity_1, v);
    let angle_u = angle_between(velocity_1, u);
    let angle_12 = angle_between(velocity_1, velocity_2);
    let angle_13 = angle_between(velocity_1, velocity_3);
    let [type_1, type_2, type_3] = [m1, m2, m3].map(|node| network.nodes[node].type_number());
    let is_opposite = |angle: f64| angle > PI - TOLERANCE && angle < PI + TOLERANCE;
    let is_parallel = |angle: f64| angle < TOLERANCE || angle > PI - TOLERANCE;

    if is_opposite(angle_u) && is_parallel(angle_12) {
        if type_1 < 300
            && type_3 < 300
            && let Some(turned) = turn_along(velocity_1, angle_v, v, TOLERANCE)
        {
            set_sector_velocity(network, m1, s1, turned);
        }
    } else if is_opposite(angle_v)
        && is_parallel(angle_13)
        && type_1 < 300
        && type_2 < 300
        && let Some(turned) = turn_along(velocity_1, angle_u, u, TOLERANCE)
    {
        set_sector_velocity(network, m1, s1, turned);
    }

    // A node in the in-flow zone: along the edge to a node out of the zone
    let is_inflow_edge_node = |type_number: u32| matches!(type_number, 310 | 312);
    let current = sector_velocity(network, m1, s1);
    let edge = if is_inflow_edge_node(type_1) && !is_inflow_edge_node(type_2) {
        Some(u)
    } else if is_inflow_edge_node(type_1) && !is_inflow_edge_node(type_3) {
        Some(v)
    } else {
        None
    };
    if let Some(edge) = edge {
        let turn = angle_between(velocity_1, edge);
        if let Some(turned) = turn_along(current, turn, edge, TOLERANCE) {
            set_sector_velocity(network, m1, s1, turned);
        }
    }
    false
}

/// Time step of each sector: 0.005 · √(Voronoi area / speed²) (`DefineTimeStep`)
fn define_time_steps(network: &mut Network) {
    const SMALLEST_SQUARED_SPEED: f64 = 1e-10;
    for node in &mut network.nodes {
        let area = node.volume / node.aperture;
        for sector in &mut node.sectors {
            let squared_speed =
                sector.velocity[0] * sector.velocity[0] + sector.velocity[1] * sector.velocity[1];
            let divisor = if squared_speed > SMALLEST_SQUARED_SPEED {
                squared_speed
            } else {
                SMALLEST_SQUARED_SPEED
            };
            sector.time_step = 0.005 * (area / divisor).sqrt();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn face(normal_area: [f64; 2]) -> FaceFlux {
        FaceFlux {
            normal_area,
            mass_flux: 1.0,
        }
    }

    #[test]
    fn faces_on_one_line_are_singular() {
        assert!(is_singular(&[face([1.0, 1.0]), face([-12.0, -12.0])]));
    }

    #[test]
    fn faces_almost_on_one_line_are_singular() {
        assert!(is_singular(&[face([1.0, 1.0]), face([-12.0, -12.0 * (1.0 + 1e-7)])]));
    }

    #[test]
    fn faces_in_two_directions_are_not_singular() {
        assert!(!is_singular(&[face([1.0, 1.0]), face([-12.0, 3.0])]));
    }
}
