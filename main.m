%% Run a graph/tensor cardiac simulation
% Edit the occupancy array and junction mesh list below, then run this file.

clear;
clc;
close all;
project_folder = fileparts(mfilename("fullpath"));

addpath(fullfile(project_folder, "config"));

cfg = default_config();

%% A three-cell 1-D cable
cfg.topology = [1,1,1];

% Junction ordering is stored in result.topology.junction_cells.
cfg.junction_mesh = 1;

cfg.model = "ORd11";
cfg.stim_cell = 1;
cfg.nbeats = 2;
cfg.stim_amp = 50;
cfg.T = cfg.BCL * cfg.nbeats;  % Plot one complete action potential without a second stimulus.
% cfg.clamp_flag = true(4, 1);
cfg.show_progress = true;
cfg.make_plots = false;

tic
result = run_graph_simulation(cfg);
toc

%% Plots
cfg = default_config();

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
    cfg.mesh_files(1),cfg.BCL,cfg.nbeats,numel(cfg.topology));

file_name_3D = sprintf("axial_phi_%s_BCL_%d_nbeats_%d_Ncell_%d_3D.png", ...
    cfg.mesh_files(1),cfg.BCL,cfg.nbeats,numel(cfg.topology));

saveas(fig1,fullfile(cfg.save_plot_path,file_name_general))
saveas(fig2,fullfile(cfg.save_plot_path,file_name_3D))