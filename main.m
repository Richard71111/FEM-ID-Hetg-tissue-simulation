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
% cfg.clamp_flag = false(4, 1);
cfg.show_progress = true;
cfg.make_plots = false;

tic
result = run_graph_simulation(cfg);
toc

%% Plots
time = result.time;
phi_axial = result.Vm_cell;

fig1 = figure;
plot(time,phi_axial)
xlabel("Time(ms)")
ylabel("$\phi_{axial}$","Interpreter", "latex")
title("reconstructed code")

time = time(:);              % Nt x 1
Ncell = size(phi_axial, 1);

[Cell, Time] = meshgrid(1:Ncell, time);

fig2 = figure;
surf(Cell, Time, phi_axial.');

xlabel("Cell number", "Interpreter", "latex");
ylabel("Time (ms)", "Interpreter", "latex");
zlabel("$\phi_{axial}$ (mV)", "Interpreter", "latex");

title("Reconstructed code", "Interpreter", "latex");
shading interp;
colorbar;
view(45, 30);

file_name_general = sprintf("axial_phi_%s_BCL_%d_nbeats_%d_Ncell_%d_general.png", ...
    cfg.mesh_files(1),cfg.BCL,cfg.nbeats,result.topology.Ncell);

file_name_3D = sprintf("axial_phi_%s_BCL_%d_nbeats_%d_Ncell_%d_3D.png", ...
    cfg.mesh_files(1),cfg.BCL,cfg.nbeats,result.topology.Ncell);

saveas(fig1,fullfile(cfg.save_plot_path,file_name_general))
saveas(fig2,fullfile(cfg.save_plot_path,file_name_3D))

% Save data
save(fullfile(cfg.save_data_path, "new_model.mat"), ...
    "time", "phi_axial");
%% compare 2 mode
old_model = load(fullfile(cfg.save_data_path, "old_model.mat"));
new_model = load(fullfile(cfg.save_data_path, "new_model.mat"));
