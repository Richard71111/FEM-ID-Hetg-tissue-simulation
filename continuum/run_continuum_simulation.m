function result = run_continuum_simulation(cfg)
%RUN_CONTINUUM_SIMULATION Build graph/FEM data and run note-style solver.
% The state is advanced in separated variables named after the continuum
% derivation, while sparse vectors are used only at solver boundaries.

continuum_folder = fileparts(mfilename("fullpath"));
project_folder = fileparts(continuum_folder);
addpath(fullfile(project_folder, "models"));
addpath(continuum_folder);

if nargin < 1
    cfg = default_continuum_config();
end
if isempty(cfg.T)
    cfg.T = cfg.BCL * cfg.nbeats;
end

cell_coordinates = [];
cell_port_count = [];
if isfield(cfg, "cell_coordinates")
    cell_coordinates = cfg.cell_coordinates;
end
if isfield(cfg, "cell_port_count")
    cell_port_count = cfg.cell_port_count;
end

topology = build_continuum_topology( ...
    cfg.adjacency_matrix, cell_coordinates, cell_port_count);
if cfg.stim_cell < 1 || cfg.stim_cell > topology.Ncell
    error("stim_cell must be between 1 and Ncell.");
end

model = setup_ionic_model(cfg);
network = assemble_continuum_network(cfg, topology, model);
result = run_time_loop_continuum(cfg, topology, model, network);
result.cfg = cfg;
result.topology = topology;
result.network = network;
end
