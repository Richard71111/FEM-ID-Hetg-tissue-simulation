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

%% Reduce to the summed quantities (order-independent, matches old runner)
Np = result.network.Npatches;
Nt = numel(result.time);
Nstate = size(result.Gstate, 1) / Np;

% sum each state over all membrane patches -> Nstate x Nt
Gsum = squeeze(sum(reshape(result.Gstate, Np, Nstate, Nt), 1));
if Nstate == 1
    Gsum = reshape(Gsum, 1, Nt);
end

% sum each cleft ion over all cleft nodes (M x Njunction) -> 4 x Nt
Ssum = squeeze(sum(sum(result.S_cleft, 2), 3));
Ssum = reshape(Ssum, 4, Nt);

new = struct();
new.time      = result.time;        % 1 x Nt
new.phi_axial = result.phi_axial;   % Ncell x Nt
new.Gsum      = Gsum;               % Nstate x Nt
new.Ssum      = Ssum;               % 4 x Nt
new.Ncell     = result.topology.Ncell;
new.Nstate    = Nstate;
new.Npatches  = Np;
new.M         = result.network.M;
new.Njunc     = result.topology.Njunction;
new.dt        = cfg.dt;
new.save_every = cfg.save_every;

save_folder = fullfile(project_folder, "Output", "Save data");
if ~exist(save_folder, "dir"); mkdir(save_folder); end
save(fullfile(save_folder, "new_model.mat"), "new");
fprintf("Saved new_model.mat (%d samples) to %s\n", Nt, save_folder);

%% quick look
fig = figure;
plot(new.time, new.phi_axial');
xlabel("Time (ms)"); ylabel("\phi_{axial} (mV)");
title("New code - axial potential"); grid on;
saveas(fig, fullfile(save_folder, "new_phi_axial.png"));
