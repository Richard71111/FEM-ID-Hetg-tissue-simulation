# Tensor/Graph Electric Simulation

This project represents cardiac cells as graph nodes and cell-cell
junctions as graph edges. Each junction keeps its own FEM cleft mesh.

```text
cell voltage         Vm_cell(Ncell, time)
ID potential         phi_ID(M, 2, Njunction)
boundary potential   phi_boundary(Nport, Ncell)
cleft potential      phi_cleft(M, Njunction)
cleft concentration  S_cleft(4, M, Njunction)
```

`M` is the number of FEM patches in one junction mesh. The second ID
dimension contains the two cells connected by that junction.

## Run

```matlab
cd("electric_sim_tensor_graph");
main
```

## Adjacency-Matrix Input

The topology input is a symmetric `Ncell`-by-`Ncell` matrix. A value of
one means that two cells share one junction.

### Three-cell cable

```matlab
cfg.adjacency_matrix = [
    0, 1, 0;
    1, 0, 1;
    0, 1, 0
];
cfg.cell_port_count = 2;
cfg.cell_coordinates = (1:3)';
```

This creates:

```text
cell 1 -- cell 2 -- cell 3
```

### Five-cell star

```matlab
cfg.adjacency_matrix = [
    0, 1, 1, 1, 1;
    1, 0, 0, 0, 0;
    1, 0, 0, 0, 0;
    1, 0, 0, 0, 0;
    1, 0, 0, 0, 0
];
cfg.cell_port_count = [4, 2, 2, 2, 2];
cfg.cell_coordinates = [
     0,  0;
     1,  0;
     0,  1;
    -1,  0;
     0, -1
];
```

The adjacency matrix defines electrical connectivity. Coordinates are
optional and are used only for plotting. Therefore arbitrary graphs,
including branches, cycles, and non-grid networks, can be simulated.

## Cell Ports

Every junction occupies one port on each connected cell. Remaining ports
become lumped boundary-disc nodes.

```matlab
cfg.cell_port_count = 2;   % Typical 1-D cable
cfg.cell_port_count = 4;   % Typical 2-D grid interpretation
cfg.cell_port_count = 6;   % Typical 3-D grid interpretation
```

The value may also be one number per cell. It must not be smaller than
the degree of that cell. If it is empty, the code uses:

```matlab
max(2, cell_degree)
```

ID-localized ionic current is divided by the total port count of its
cell. This preserves the original 1-D behavior: each terminal cell in a
cable has one connected junction port and one boundary port.

## Junction Ordering and Meshes

`build_topology` stores the edge list in:

```matlab
result.topology.junction_cells
```

Each column contains:

```text
[first_cell; second_cell]
```

Edges are generated from the upper triangle of the adjacency matrix.
Assign one mesh index to each edge:

```matlab
cfg.junction_mesh = [1, 2, 1, 2];
```

A scalar applies the same mesh to every junction:

```matlab
cfg.junction_mesh = 1;
```

All selected meshes must currently have the same number `M` of FEM
patches.

## Main Files

```text
config/default_config.m
    Model parameters, adjacency matrix, coordinates, and port counts.

topology/build_topology.m
    Validates the adjacency matrix and creates cells, edges, and ports.

matrix/assemble_graph_network.m
    Loads junction meshes and assembles sparse global operators.

simulation/run_time_loop.m
    Advances ionic states, cleft concentrations, and voltages.

run_graph_simulation.m
    Runs topology setup, matrix assembly, simulation, and plotting.
```

## Concentration Modes

Dynamic cleft concentrations:

```matlab
cfg.clamp_flag = false(4, 1);
```

Fixed cleft concentrations:

```matlab
cfg.clamp_flag = true(4, 1);
```

## Results

```matlab
result.topology.adjacency
result.topology.junction_cells
result.topology.junction_ports
result.final.phi_cell
result.final.phi_ID
result.final.phi_boundary
result.final.phi_cleft
result.final.S_cleft
```

The sparse flat vector remains internal to
`matrix/assemble_graph_network.m`; cell and junction results are exposed
with stable tensor dimensions.
