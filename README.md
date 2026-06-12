# Tensor/Graph Electric Simulation

This project represents cells and junctions with physical tensors instead
of exposing one manually indexed node vector.

```text
cell potential       phi_cell(Ncell)
ID potential         phi_ID(M, 2, Njunction)
boundary potential   phi_boundary(Nface, Ncell)
cleft potential      phi_cleft(M, Njunction)
cleft concentration  S_cleft(4, M, Njunction)
```

`M` is the number of FEM patches in one junction mesh. The second ID
dimension contains the two cells facing the same cleft.

The sparse global vector still exists inside
`matrix/assemble_graph_network.m`, because electrically connected cells
must be solved together. All flattening is centralized there. The main
simulation only uses the tensor layouts and does not construct `P`, `Q`,
or `C` matrices.

## Run

```matlab
cd("electric_sim_tensor_graph");
main
```

The example in `main.m` runs a three-cell 1-D cable and plots:

1. Cell membrane-voltage traces.
2. The graph topology colored by final voltage.
3. Mean cleft potential at every junction.

## Topology Input

Occupied entries are cells. Orthogonally adjacent entries create
junctions.

### 1-D cable

```matlab
cfg.topology = [1, 1, 1];
```

This creates three cells and two junctions.

### 2-D cross

```matlab
cfg.topology = [
    0, 1, 0;
    1, 1, 1;
    0, 1, 0
];
```

This creates five cells and four junctions.

### 3-D topology

```matlab
mask = zeros(3, 3, 2);
mask(2, 2, 1) = 1;
mask(1, 2, 1) = 1;
mask(3, 2, 1) = 1;
mask(2, 2, 2) = 1;
cfg.topology = mask;
```

The same setup and simulation files are used for 1-D, 2-D, and 3-D
topologies.

## Junction Ordering and Meshes

`build_topology` creates:

```matlab
topology.junction_cells
topology.junction_faces
```

Each column is:

```text
[first_cell; second_cell]
```

Assign one mesh index to each column:

```matlab
cfg.junction_mesh = [1, 2, 1, 2];
```

The indices refer to:

```matlab
cfg.mesh_files
```

A scalar applies the same mesh to every junction:

```matlab
cfg.junction_mesh = 1;
```

All selected meshes must currently contain the same number `M` of FEM
patches so they can share one dense junction tensor.

The default mesh is `FEMDATA_baseline.mat`, which is the same baseline
mesh used by the original 1-D source code.

## Main Files

```text
config/default_config.m
    User parameters and topology input.

topology/build_topology.m
    Converts the occupancy array into cells, coordinates, and junctions.

matrix/assemble_graph_network.m
    Loads each junction mesh, creates tensor layouts, and assembles the
    hidden global conductance and capacitance operators.

simulation/run_time_loop.m
    Contains one direct while loop for ionic currents, cleft
    concentrations, and the global voltage solve.

run_graph_simulation.m
    Runs setup, simulation, and plotting.
```

## Concentration Modes

For dynamic cleft concentrations:

```matlab
cfg.clamp_flag = [false; false; false; false];
```

For a faster electrical demonstration with fixed bulk concentrations:

```matlab
cfg.clamp_flag = true(4, 1);
```

## Accessing Tensor Results

```matlab
result.final.phi_cell
result.final.phi_ID
result.final.phi_boundary
result.final.phi_cleft
result.final.S_cleft
```

Adding cells changes only the cell and junction dimensions. Local cleft
indices remain `1:M` for every junction.

## Cell Faces and Boundary Membrane

The graph uses `2*dimension` possible membrane faces per cell:

```text
1-D: 2 faces
2-D: 4 faces
3-D: 6 faces
```

A connected face contains the detailed `M`-patch ID mesh. An unconnected
boundary face contains one lumped membrane node, matching the terminal
disc construction in the original 1-D code. It is connected to the cell
node by one myoplasmic conductance and is not represented by `M` nodes.

Face order is `[-axis1, +axis1, -axis2, +axis2, ...]`. A zero in
`network.layout.boundary(face, cell)` means that the face is connected to
another cell. Channel localization is divided by the number of possible
faces, not by the number of connected neighbors. In a 1-D cable, an end
cell therefore keeps half of its ID-localized current on its external
terminal disc.
