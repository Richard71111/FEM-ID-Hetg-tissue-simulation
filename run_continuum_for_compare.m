%% RUN_CONTINUUM_FOR_COMPARE
% Run the note-style continuum implementation on the tensor benchmark.
%
% This is the main entry point for producing continuum_model.mat, which is
% then consumed by compare_continuum_with_tensor.m.

clear;
clc;
close all;

project_folder = fileparts(mfilename("fullpath"));
addpath(fullfile(project_folder, "continuum"));

cfg = default_continuum_config();

%% Match the working tensor graph benchmark exactly.
cfg.adjacency_matrix = [
    0, 1, 0;
    1, 0, 1;
    0, 1, 0
];
cfg.cell_coordinates = (1:3)';
cfg.cell_port_count = 2;
cfg.junction_mesh = 1;

cfg.model = "ORd11";
cfg.stim_cell = 1;
cfg.nbeats = 2;
cfg.stim_amp = 50;
cfg.mesh_files = ["FEMDATA_2.mat"];
cfg.D = 1;
cfg.T = cfg.BCL * cfg.nbeats;
cfg.clamp_flag = false(4, 1);
cfg.show_progress = true;
cfg.make_plots = false;

tic;
result = run_continuum_simulation(cfg);
toc;

%% Reduce to order-independent quantities used by compare_old_new.m.
Np = result.network.Npatches;
Nt = numel(result.time);
Nstate = size(result.Gstate, 1) / Np;

Gsum = squeeze(sum(reshape(result.Gstate, Np, Nstate, Nt), 1));
if Nstate == 1
    Gsum = reshape(Gsum, 1, Nt);
end

Ssum = squeeze(sum(sum(result.S_cleft, 2), 3));
Ssum = reshape(Ssum, 4, Nt);

continuum = struct();
continuum.time = result.time;
continuum.phi_axial = result.phi_axial;
continuum.Gsum = Gsum;
continuum.Ssum = Ssum;
continuum.Ncell = result.topology.Ncell;
continuum.Nstate = Nstate;
continuum.Npatches = Np;
continuum.M = result.network.M;
continuum.Njunc = result.topology.Njunction;
continuum.dt = cfg.dt;
continuum.save_every = cfg.save_every;

save_folder = fullfile(project_folder, "Output", "Save data");
if ~exist(save_folder, "dir")
    mkdir(save_folder);
end
save(fullfile(save_folder, "continuum_model.mat"), "continuum");
fprintf("Saved continuum_model.mat (%d samples) to %s\n", Nt, save_folder);
