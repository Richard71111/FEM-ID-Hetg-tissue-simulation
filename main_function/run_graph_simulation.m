function result = run_graph_simulation(cfg)
%RUN_GRAPH_SIMULATION Build topology, assemble operators, and run the model.
% Input: optional cfg structure; default_config is used when omitted.
% Output: result with graph metadata and tensor-shaped physical variables.

function_folder = fileparts(mfilename("fullpath"));
project_folder  = fileparts(function_folder);
addpath(fullfile(project_folder, "config"));
addpath(fullfile(project_folder, "topology"));
addpath(fullfile(project_folder, "models"));
addpath(fullfile(project_folder, "matrix"));
addpath(fullfile(project_folder, "simulation"));
addpath(fullfile(project_folder, "plot_function"));

if nargin < 1
    cfg = default_config();
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
topology = build_topology( ...
    cfg.adjacency_matrix, cell_coordinates, cell_port_count);
if cfg.stim_cell < 1 || cfg.stim_cell > topology.Ncell
    error("stim_cell must be between 1 and Ncell.");
end

model = setup_ionic_model(cfg);
network = assemble_graph_network(cfg, topology, model);
result = run_time_loop(cfg, topology, model, network);
result.cfg = cfg;
result.topology = topology;
result.network = network;

if cfg.make_plots
    plot_graph_result(result);
end
end
