# dfntrans-rs

Particle tracking in a discrete fracture network (DFN), in Rust. It is a
port of DFNTrans, the particle tracking of
[dfnWorks](https://github.com/lanl/dfnWorks) (LANL, C), for the work flow

    dfngen-rs (dfn_mesh.vtu) → fehm-rs (pressures, mass fluxes) → dfntrans-rs (particles)

The flow solution comes from [fehm-rs](../fehm-rs), on a triangle mesh of
[dfngen-rs](../dfngen-rs). The velocity at each node is reconstructed from
the fluxes of its control volume; the particles move with a
predictor-corrector scheme in the plane of each fracture, and an
intersection rule selects the next fracture.

## Usage

```bash
cargo run --release -- <input_file.toml> <output_folder>
```

Example: `cargo run --release -- tests/data/crossing_fractures.toml output`.

The binary reads the fehm-rs input file which the input file names, solves
the flow with fehm-rs, and tracks the particles.

## Input file

One TOML file. An unknown key is an error.

```toml
flow = "flow.toml"              # fehm-rs input file; a relative path is relative to the folder of this file
seed = 12345                    # random numbers (srand48 of DFNTrans); 0 = from the clock
time_unit = "days"              # seconds, minutes, hours, days or years: unit of the times and velocities
porosity = 1.0                  # fracture porosity, above 0, at most 1
max_time_steps = 1000000        # time steps of a particle, at most
inflow = ["west"]               # boundary__<name> point arrays of the mesh: the in-flow boundary
outflow = ["east"]              # the out-flow boundary (a point in both is in the in-flow boundary)
routing = "complete_mixing"     # or "streamline_routing": the rule at an intersection

[seeding]                       # the initial positions; one method
method = "flux_weighted"
particle_count = 1000

[tdrw]                          # optional: matrix diffusion (time domain random walk)
matrix_porosity = 0.02
matrix_diffusivity = 1e-13      # m²/s

[control_planes]                # optional: record the particles where they cross planes
shape = "plane"
direction = "x"                 # the planes are normal to x, y or z
spacing = 1.0                   # m, from the in-flow boundary to the out-flow boundary

[output]                        # optional section
trajectories = false            # optional, default false: trajectories.vtp
all_particles = false           # optional, default false: also the particles which do not exit
velocities = false              # optional, default false: velocities.vtp
```

`inflow` and `outflow` name `boundary__<name>` point arrays of the mesh;
they need not be pressure boundaries of the flow input.

### Seeding methods

The in-flow lines are the chains of mesh edges between two in-flow points
in one fracture.

| `method` | Keys | Particles | Flux weight | DFNTrans |
|---|---|---|---|---|
| `uniform_per_line` | `particles_per_line`, `weight_by_flux` (optional, default true) | the same number on each in-flow line, at equal distances | in-flow flux at the position, or 0 | `init_nf` |
| `equidistant` | `particles_per_line`, `weight_by_flux` (optional, default true) | on all in-flow lines at one distance: total length / (`particles_per_line` · number of lines) | in-flow flux at the position, or 0 | `init_eqd` |
| `random` | `particle_count` | at the centers of random triangles without an out-flow point | aperture | `init_random` |
| `flux_weighted` | `particle_count` | on the in-flow lines: at least one for each point, and one for each part `total flux / particle_count` of its flux | `total flux / particle_count` (m³/s) | `init_fluxw` |
| `at_nodes` | `particles_per_node` | at each in-flow point | aperture | `init_well` |

A particle which starts exactly at a point of a fracture edge can move
along the edge and leave the fracture (also in DFNTrans: at an in-flow
point on a corner of the edge, DFNTrans turns the velocity along the
edge). On a dfngen-rs mesh, use `flux_weighted`, `uniform_per_line` or
`equidistant`.

### Control planes and cylinders

```toml
[control_planes]
shape = "cylinder"
axis = "z"                      # cylinders about the line through the origin along z
spacing = 1.0                   # m
length_to_well = 10.0           # m, radius of the outer cylinder
well_thickness = 0.5            # m, the cylinders stop at radius well_thickness / 2
```

The planes start at the coordinate of the first in-flow point, in the
direction of the first out-flow point. A particle is recorded at its start
(planes), at the first time step at or behind each plane or cylinder, and at
its exit.

## Output files

All files are in `<output_folder>`. Times are in `time_unit`.

| File | Content | DFNTrans |
|---|---|---|
| `travel_times.dat` | One line for each particle which exited (with `all_particles`, for each particle): particle number, last time step, flux weight, travel time, final x y z, β, path length. With `[tdrw]`: advective time, advective + diffusion time and diffusion time in place of the travel time. | `partime.dat` |
| `inflow_flux.dat` | The volume flux into the network through the in-flow boundary (m³/s) | `inputflux_m3s` |
| `trajectories.vtp` | With `output.trajectories`: one polyline for each particle, a point for each time step and each intersection crossing. Point data: `time`, `velocity`, `fracture_id`, `aperture`, `pressure`. Cell data: `particle`, `exited`. | `traject_*`, `part_*.inp` |
| `control_planes.vtp` | With `[control_planes]`: one vertex for each crossing. Point data as `trajectories.vtp`; cell data `particle`, `length` (and `total_time`, `diffusion_time` with `[tdrw]`). | `part_control_*` |
| `velocities.vtp` | With `output.velocities`: one vertex for each velocity of a point (a point on an intersection has more than one). Point data: `velocity` (3D), `node`, `fracture_id`, `time_step`. | `Velocity3D` |

β = Σ 2 · step length / (speed · aperture) along the path.

## Library

The crate is a library (`dfntrans`) and a binary. The library compiles for
`wasm32-unknown-unknown`: `build_input_from_text` takes the text of the
input file, `track_particles` takes the fehm-rs input and solution, and
`write_output` gets a writer for each output file. `read_input` and
`write_output_to_folder` are the file-system wrappers.

```rust
let input = dfntrans::build_input_from_text(&toml_text)?;
let transport = dfntrans::track_particles(&input, &flow_input, &flow_solution)?;
```

## Differences from DFNTrans

The data which DFNTrans reads from the LaGriT, FEHM and dfnGen files come
from the mesh and the flow solution: the node types (interior, fracture
edge, intersection) from the triangles of each edge, the rotation of each
fracture from the mean normal of its triangles, and the density of each
node from the water properties of fehm-rs (DFNTrans has one density).

Deliberate changes:

- A node on three or more fractures (a crossing of two intersection lines,
  possible in a dfngen-rs mesh) gets one velocity for each group of
  triangles between two intersection edges. DFNTrans supports two
  fractures for each node.
- The two halves of the control volume of a node on an intersection: a face
  within 1e-10 (relative) of the intersection line is in both halves. In
  DFNTrans, the other edge of a straight line is in one or both halves by
  the rounding of its single precision rotation and of fused
  multiply-add; in exact arithmetic it is in both.
- The line of the halves is the last face whose triangles are in two
  fractures (an intersection edge). DFNTrans takes the last face to a node
  on an intersection, which is an edge of the fracture outline if the
  intersection ends at the outline and the outline neighbor is on another
  intersection (frequent in a network without FRAM). Then one half has only
  the two outline faces, on one line.
- A sector whose faces have one direction (the determinant of the least
  squares fit is at most 1e-10 of the product of its diagonal) gets the fit
  of all faces of the node in its fracture, or the velocity 0 if these also
  have one direction, with a warning. DFNTrans sets the determinant to
  1e-75 and gets a velocity of about 1e40.
- The line of the halves uses the coordinates of the neighbor in the
  fracture of the half. DFNTrans uses the coordinates of the slot of the
  node (`coord_xy[ic]`), which is another fracture if the neighbor is on
  another intersection line.
- At a fracture edge, two edges of equal length within 1e-10 give the
  normal of the second edge (DFNTrans: the longer edge; for equal lengths,
  the rounding decides).
- The search of the first cell of a particle uses only the cells of its fracture.
- Streamline routing compares the fractures of the two out-going cells;
  DFNTrans compares those of the next cells (`cell[c]` for `cell[c-1]`).
- No limit of the trajectory buffer: DFNTrans stops a particle with output
  files after 0.8 · `max_time_steps` steps between two intersections.
- The seeding methods on the in-flow boundary use the in-flow lines of
  any shape; DFNTrans uses the segment between the first and the last
  in-flow point (by number) of each fracture.
- Bad input gives an error; DFNTrans prints "Error" and continues.

Kept as in DFNTrans (they change the results):

- The aperture of a node is the aperture of its first fracture: DFNTrans
  reads the apertures before it finds the second fracture of a node.
- `DefineAngle` uses the single precision `acosf`; π is 3.14159265359.
- The random numbers of `drand48`; the same seed gives the same numbers.

Not ported: PFLOTRAN input, `init_oneregion`, `init_matrix`, dispersivity
(`out_disp`), MARFA and PLUMECALC files, tortuosity and fracture list
files, AVS files, `out_curv`, `out_filetemp`, `out_grid`.

## Comparison with DFNTrans

`tests/dfntrans_reference.rs` compares the paths of 90 particles in the
crossing fractures (`tests/data/reference`) with the output of DFNTrans:
complete mixing, streamline routing and matrix diffusion. Each particle
makes the same random decisions; the travel times differ by less than
4e-5 (the density).

To repeat a DFNTrans run, build DFNTrans from the dfnWorks sources (the
build without fused multiply-add is that of gcc on x86-64):

```bash
make CFLAGS="-lm -w -O3 -std=gnu89 -Wno-error=return-mismatch -ffp-contract=off"
scripts/compare_with_dfntrans.sh tests/data/reference/complete_mixing.toml <work_folder> <path to DFNTrans>
```

`examples/export_dfntrans.rs` writes the DFNTrans input files of a case
(`full_mesh.inp`, `full_mesh.stor`, `run.fin`, `boundary.zone`,
`params.txt`, `aperture.dat`, `control.dat`), with the mesh sorted by
fracture as DFNTrans needs, and the dfntrans-rs output for the same mesh.

On a dfngen-rs mesh (31 129 nodes, 88 fractures), the node velocities of
dfntrans-rs and of DFNTrans agree within 1e-3 at all 23 787 interior nodes
and at all but 9 of 37 753 velocities, if DFNTrans gets the changes
above of the halves and of the edge lengths (single precision tolerance
1e-6). Without them, about 3000 velocities at intersection nodes differ.

Found in DFNTrans, not a part of the port: `init_random` calls
`Moving2Center` and `InsideCell` with the global particle index `np`, thus
its particles have no interpolation weights; `DarcyVelocity` reads past the
four cells of an edge (`while (fracts[j][l] != 0)`); a `.fin` token longer
than 15 characters overflows a buffer.

## Test data

`tests/data`: the meshes `single_fracture.vtu` and `crossing_fractures.vtu`
of fehm-rs (`scripts/make_test_meshes.py` there). `single_fracture_flow.toml`
has no gravity, thus a uniform flow. `tests/data/reference`:
`crossing_fractures.vtu` with the array `boundary__west_inner`, the west
points without the two corners (`scripts/add_inner_boundary.py`), and the
`partime.dat` files of DFNTrans.

`tests/dfngen_chain.rs` makes a network with the dfngen-rs library (the case
`mesh` of dfngen-rs with 300 fractures, the largest cluster, seed 3), solves
the flow from the face x = -5 to the face x = 5 and tracks 665 particles;
all exit. Run the tests with `--release`.

`scripts/add_face_boundary.py` adds a boundary array for one face of the
domain to a dfngen-rs mesh, from its array `boundary__mesh_surface`.

## License

LGPL-3.0, as dfnWorks (`LICENSE.md`, `LICENSE_DOE_LANL.md`). DFNTrans:
LANL LA-CC-14-091.
