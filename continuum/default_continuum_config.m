function cfg = default_continuum_config()
%DEFAULT_CONTINUUM_CONFIG Define graph, model, mesh, and time parameters.
% This config mirrors the tensor graph benchmark, but the continuum runner
% consumes adjacency_matrix directly and keeps note-style variables exposed.

continuum_folder = fileparts(mfilename("fullpath"));
project_folder = fileparts(continuum_folder);

%% Topology
cfg.adjacency_matrix = [ ...
    0, 1, 0; ...
    1, 0, 1; ...
    0, 1, 0];
cfg.cell_coordinates = [];
cfg.cell_port_count = 2;
cfg.junction_mesh = 1;

%% FEM mesh files
cfg.mesh_folder = fullfile(project_folder, "data", "384");
cfg.mesh_files = ["FEMDATA_1.mat", "FEMDATA_1.mat"];
cfg.scale_gj_loc = 1;
cfg.scale_chan_loc = 1;

%% Ionic model and cell geometry
cfg.model = "ORd11";
cfg.L = 100;
cfg.r = 11;
cfg.Cm = 1e-8;

%% Current localization
cfg.locINa = 0.7;
cfg.locIK1 = 0.2;
cfg.locICa = 0.2;
cfg.locINaK = 0.2;

%% Electrical properties
cfg.rho_myo = 1500;
cfg.rho_ext = 1500;
cfg.rho_ie = 1;
cfg.ggap = [];
cfg.D = 0.1;
cfg.f_cleft = 1;
cfg.f_bulk = 1;
cfg.fVol = 1;

%% Bulk concentrations and physical constants
cfg.Na_b = 140;
cfg.K_b = 5.4;
cfg.Ca_b = 1.8;
cfg.A_b = cfg.Na_b + cfg.K_b + 2 * cfg.Ca_b;
cfg.clamp_flag = true(4, 1);
cfg.F = 96.5;
cfg.R = 8.314;
cfg.Temp = 310;

%% Time and stimulus
cfg.BCL = 1000;
cfg.nbeats = 1;
cfg.T = [];
cfg.adaptive_dt = true;
cfg.twin = 50;
cfg.dt = 0.01;
cfg.dt2 = 0.1;
cfg.dtS = cfg.dt / 5;
cfg.dtS2 = cfg.dt2 / 10;
cfg.save_every = 10;
cfg.electric_integrator = "backward_euler";  % "backward_euler", "forward_euler", or "crank_nicolson".
cfg.cleft_integrator = "forward_euler_split";
cfg.stim_cell = 1;
cfg.stim_dur = [];
cfg.stim_amp = [];

%% Output
cfg.show_progress = true;
cfg.progress_interval = 500;
cfg.make_plots = true;

cfg.save_plot_path = fullfile(project_folder, "Output", "Save plot");
cfg.save_data_path = fullfile(project_folder, "Output", "Save data");

folders_to_make = {cfg.save_plot_path, cfg.save_data_path};
for k = 1:numel(folders_to_make)
    if ~exist(folders_to_make{k}, "dir")
        mkdir(folders_to_make{k});
    end
end
end
