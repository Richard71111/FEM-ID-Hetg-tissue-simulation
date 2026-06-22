function result = run_graph_voltage_clamp_simulation(cfg)
%RUN_GRAPH_VOLTAGE_CLAMP_SIMULATION Build and run a voltage-clamp protocol.
% This runner mirrors run_graph_simulation but calls the voltage-clamp time
% loop so the original simulation code remains unchanged.

project_folder = fileparts(mfilename("fullpath"));
addpath(fullfile(project_folder, "config"));
addpath(fullfile(project_folder, "topology"));
addpath(fullfile(project_folder, "models"));
addpath(fullfile(project_folder, "matrix"));
addpath(fullfile(project_folder, "simulation"));

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
result = run_time_loop_voltage_clamp(cfg, topology, model, network);
result.cfg = cfg;
result.topology = topology;
result.network = network;

if cfg.make_plots
    if isfield(cfg, "vclamp") && isfield(cfg.vclamp, "figure_name") && ...
            strlength(string(cfg.vclamp.figure_name)) > 0
        figure_name = string(cfg.vclamp.figure_name);
    else
        figure_name = default_voltage_clamp_figure_name(cfg, topology.Ncell);
    end
    plot_voltage_clamp_result(result, figure_name);
end
end

function figure_name = default_voltage_clamp_figure_name(cfg, Ncell)
protocol = cfg.vclamp;
if isfield(protocol, "cells") && ~isempty(protocol.cells)
    clamped_cells = unique(protocol.cells(:));
elseif isfield(protocol, "mode") && any(lower(string(protocol.mode)) == ["two_cell", "both_cells"])
    clamped_cells = (1:min(2, Ncell)).';
else
    clamped_cells = 1;
end

if numel(clamped_cells) == 1
    cell_text = sprintf("cell %d clamped", clamped_cells);
elseif numel(clamped_cells) == Ncell
    cell_text = "all cells clamped";
else
    cell_text = "cells " + strjoin(string(clamped_cells.'), "-") + " clamped";
end

if isfield(protocol, "release_after_step") && logical(protocol.release_after_step)
    protocol_text = "release after pulse";
else
    protocol_text = "return to rest";
end

figure_name = "Voltage clamp " + protocol_text + ": " + cell_text;
end
