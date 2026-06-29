function result = run_graph_voltage_clamp_dynamic_simulation(cfg)
%RUN_GRAPH_VOLTAGE_CLAMP_DYNAMIC_SIMULATION Build and run a waveform clamp.
% Mirrors run_graph_voltage_clamp_simulation but drives the dynamic time
% loop (run_time_loop_voltage_clamp_dynamic) that accepts an arbitrary
% sampled voltage waveform in cfg.vclamp.command_voltage.

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
if ~isfield(cfg, "vclamp") || ~isstruct(cfg.vclamp)
    error("cfg.vclamp must define the voltage-clamp protocol.");
end

cell_coordinates = [];
cell_port_count  = [];
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

model   = setup_ionic_model(cfg);
network = assemble_graph_network(cfg, topology, model);
result  = run_time_loop_voltage_clamp_dynamic(cfg, topology, model, network);
result.cfg      = cfg;
result.topology = topology;
result.network  = network;

if cfg.make_plots
    if isfield(cfg, "vclamp") && isfield(cfg.vclamp, "figure_name") && ...
            strlength(string(cfg.vclamp.figure_name)) > 0
        figure_name = string(cfg.vclamp.figure_name);
    else
        figure_name = "Voltage clamp dynamic waveform";
    end
    plot_voltage_clamp_result(result, figure_name);
end
end
