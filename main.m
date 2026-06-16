%% Run a graph/tensor cardiac simulation
% Edit the adjacency matrix and junction mesh list below, then run this file.

clear;
clc;
close all;
project_folder = fileparts(mfilename("fullpath"));

addpath(fullfile(project_folder, "config"));

cfg = default_config();

%% A three-cell 1-D cable
cfg.adjacency_matrix = [
    0, 1, 0;
    1, 0, 1;
    0, 1, 0
];
cfg.cell_coordinates = (1:3)';  % Used only to draw the graph.
cfg.cell_port_count = 2;  % Two membrane ports per cell for a 1-D cable.

% Junction ordering is stored in result.topology.junction_cells.
cfg.junction_mesh = 1;

cfg.model = "ORd11";
cfg.stim_cell = 1;
cfg.nbeats = 2;
cfg.stim_amp = 50;
cfg.T = cfg.BCL * cfg.nbeats;  % Plot one complete action potential without a second stimulus.
cfg.clamp_flag = false(4, 1);
cfg.show_progress = true;
cfg.make_plots = true;  % run_graph_simulation calls plot_graph_result.

tic
result = run_graph_simulation(cfg);
toc

%% Save data
% Only the retained variables are stored:
%   phi_axial : axial (intracellular) potential, Ncell-by-Nt.
%   Gstate    : ionic state vector incl. intracellular concentration,
%               (Nstate*Npatches)-by-Nt.
%   Icleft    : cleft (axial-to-ID) current per junction side, 2-by-Njunction-by-Nt.
%   S_cleft   : cleft (extracellular) concentration, 4-by-M-by-Njunction-by-Nt.
% time = result.time;
% phi_axial = result.phi_axial;
% Gstate = result.Gstate;
% Icleft = result.Icleft;
% S_cleft = result.S_cleft;
% 
% save(fullfile(cfg.save_data_path, "new_model.mat"), ...
%     "time", "phi_axial", "Gstate", "Icleft", "S_cleft");
