function cfg = default_config()
%DEFAULT_CONFIG Define topology, model, mesh, and time parameters.
% Input: none.
% Output: cfg structure used by the graph/tensor simulation.

config_folder = fileparts(mfilename("fullpath"));
project_folder = fileparts(config_folder);

%% Topology
cfg.adjacency_matrix = [ ...
    0, 1, 0; ...
    1, 0, 1; ...
    0, 1, 0];  % Symmetric cell-cell connectivity matrix.
cfg.cell_coordinates = [];  % Optional Ncell-by-1, 2, or 3 plotting coordinates.
cfg.cell_port_count = 2;  % Scalar/vector total ports; empty uses max(2, degree).
cfg.junction_mesh = 1;  % Scalar or one mesh index for each junction.

%% FEM mesh files
% The first mesh is the baseline mesh used by the original 1-D code.
cfg.mesh_folder = fullfile(project_folder, "data","384");
cfg.mesh_files = ["FEMDATA_1.mat", "FEMDATA_1.mat"];
cfg.scale_gj_loc = 1;  % Scalar or one GJ localization scale per junction.
cfg.scale_chan_loc = 1;  % Scalar or one channel scale per junction.

%% Ionic model and cell geometry
cfg.model = "ORd11";  % "LR1", "ORd11", or "Court98".
cfg.L = 100;  % Reference cell length, um.
cfg.r = 11;  % Reference cell radius, um.
cfg.Cm = 1e-8;  % Membrane capacitance density, uF/um^2.

%% Current localization
%Fraction of the total channel density localized to the ID membrane.
cfg.locINa = 0.7; 
cfg.locIK1 = 0.2;
cfg.locICa = 0.2;
cfg.locINaK = 0.2;

%% Electrical properties
cfg.rho_myo = 1500;  % Myoplasmic resistivity, kOhm*um.
cfg.rho_ext = 1500;  % Extracellular resistivity, kOhm*um.
cfg.rho_ie = 1;  % Ratio of ID resistivity to cleft resistivity.
cfg.ggap = [];  % Total GJ conductance, mS; empty uses 7.35e-4*D.
cfg.D = 0.1;  % Effective diffusion coefficient, cm^2/s.
cfg.f_cleft = 1;
cfg.f_bulk = 1;
cfg.fVol = 1;

%% Bulk concentrations and physical constants
cfg.Na_b = 140;
cfg.K_b = 5.4;
cfg.Ca_b = 1.8;
cfg.A_b = cfg.Na_b + cfg.K_b + 2 * cfg.Ca_b;
cfg.clamp_flag = true(4, 1);  % Na, K, Ca, and anion.
cfg.F = 96.5;  % C/mmol.
cfg.R = 8.314;  % J/(mol*K).
cfg.Temp = 310;  % K.

%% Time and stimulus
cfg.BCL = 1000;  % Basic cycle length, ms.
cfg.nbeats = 1;
cfg.T = [];  % Total time, ms; empty uses BCL*nbeats.

% Adaptive (dual) time step, matching the original 1-D source
% (run_generic_RCL_multi_mesh.m). Within the first `twin` ms after each
% beat onset (mod(t, BCL) < twin) the solver uses the fine step dt; for the
% rest of each cycle it uses the coarse step dt2. Set adaptive_dt = false to
% use the fine step cfg.dt everywhere.
cfg.adaptive_dt = true;  % Enable dual fine/coarse time stepping.
cfg.twin = 50;  % Fine-step window after each beat onset, ms.
cfg.dt = 0.01;  % Fine voltage time step (used within twin), ms.
cfg.dt2 = 0.1;  % Coarse voltage time step (used outside twin), ms.
cfg.dtS = cfg.dt/5;  % Fine cleft concentration time step, ms.
cfg.dtS2 = cfg.dt2/10;  % Coarse cleft concentration time step, ms.
cfg.save_every = 10;  % Cross-step saving: store outputs every N accepted steps.
cfg.stim_cell = 1;  % Cell receiving the stimulus.
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
