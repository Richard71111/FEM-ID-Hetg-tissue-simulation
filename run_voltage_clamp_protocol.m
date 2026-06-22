%% RUN_VOLTAGE_CLAMP_PROTOCOL
% Two-cell voltage-clamp examples without modifying the original solvers.
%
% Protocol family 1: hold at Vrest, pulse to V_step, then hold at Vrest.
% Protocol family 2: hold at Vrest, pulse to V_step, then release the clamp.
% Each family is run for cell 1 only and for both cells.

clear;
clc;
close all;

project_folder = fileparts(mfilename("fullpath"));
addpath(project_folder);
addpath(fullfile(project_folder, "config"));
addpath(fullfile(project_folder, "topology"));
addpath(fullfile(project_folder, "models"));
addpath(fullfile(project_folder, "matrix"));
addpath(fullfile(project_folder, "simulation"));

base_cfg = default_config();

%% Two-cell cable setup
base_cfg.adjacency_matrix = [
    0, 1;
    1, 0
];
base_cfg.cell_coordinates = (1:2)';
base_cfg.cell_port_count = 2;
base_cfg.junction_mesh = 1;

base_cfg.model = "ORd11";
base_cfg.stim_cell = 1;
base_cfg.stim_amp = 0;       % No stimulus current during voltage clamp.
base_cfg.stim_dur = 0;
base_cfg.nbeats = 1;
base_cfg.clamp_flag = false(4, 1);
base_cfg.show_progress = true;
base_cfg.make_plots = false;
base_cfg.save_every = 10;

% Voltage-clamp protocols are defined by their own command waveform rather
% than by beat/stimulus timing. With adaptive_dt=false, the time loop uses
% cfg.dt everywhere; cfg.dt2 is ignored except that we keep it equal to cfg.dt
% to make the intent obvious.
base_cfg.adaptive_dt = false;
base_cfg.dt = 0.01;
base_cfg.dt2 = base_cfg.dt;
base_cfg.dtS = base_cfg.dt / 5;
base_cfg.dtS2 = base_cfg.dt2 / 10;

%% Shared protocol timing
Vrest = -0.879989146999539e2;                 % mV.
t_rest = 100;                % ms at Vrest before the pulse.
t_tail = 100;                % ms after the longest pulse.

%% Protocol family 1: force the voltage back to Vrest after the pulse.
release_after_step = false;

cfg_return_one = make_one_cell_cfg( ...
    base_cfg, Vrest, t_rest, 20, 50, t_tail, release_after_step);
run_voltage_clamp_case( ...
    cfg_return_one, ...
    "return-to-rest, cell 1 clamped", ...
    "voltage_clamp_return_one_cell", ...
    "Voltage clamp return to rest: cell 1 clamped");

cfg_return_two = make_two_cell_cfg( ...
    base_cfg, Vrest, t_rest, [20; 30], [50; 60], t_tail, ...
    release_after_step);
run_voltage_clamp_case( ...
    cfg_return_two, ...
    "return-to-rest, both cells clamped", ...
    "voltage_clamp_return_two_cell", ...
    "Voltage clamp return to rest: both cells clamped");

%% Protocol family 2: release the voltage clamp after the pulse.
release_after_step = true;

cfg_release_one = make_one_cell_cfg( ...
    base_cfg, Vrest, t_rest, 20, 50, t_tail, release_after_step);
run_voltage_clamp_case( ...
    cfg_release_one, ...
    "release after pulse, cell 1 clamped", ...
    "voltage_clamp_release_one_cell", ...
    "Voltage clamp release after pulse: cell 1 clamped");

cfg_release_two = make_two_cell_cfg( ...
    base_cfg, Vrest, t_rest, [20; 30], [50; 60], t_tail, ...
    release_after_step);
run_voltage_clamp_case( ...
    cfg_release_two, ...
    "release after pulse, both cells clamped", ...
    "voltage_clamp_release_two_cell", ...
    "Voltage clamp release after pulse: both cells clamped");

function result = run_voltage_clamp_case(cfg, case_label, save_name, figure_name)
fprintf("\nRunning voltage clamp case: %s.\n", case_label);
tic
result = run_graph_voltage_clamp_simulation(cfg);
toc
save_voltage_clamp_result(result, save_name);
plot_voltage_clamp_result(result, figure_name);
end

function cfg = make_one_cell_cfg( ...
    base_cfg, Vrest, t_rest, V1, t1, t_tail, release_after_step)
cfg = base_cfg;
cfg.vclamp = struct();
cfg.vclamp.mode = "one_cell";
cfg.vclamp.cells = 1;
cfg.vclamp.Vrest = Vrest;
cfg.vclamp.t_rest = t_rest;
cfg.vclamp.V1 = V1;
cfg.vclamp.t1 = t1;
cfg.vclamp.t_tail = t_tail;
cfg.vclamp.release_after_step = release_after_step;
cfg.T = t_rest + t1 + t_tail;
end

function cfg = make_two_cell_cfg( ...
    base_cfg, Vrest, t_rest, V_step, t_step, t_tail, release_after_step)
cfg = base_cfg;
cfg.vclamp = struct();
cfg.vclamp.mode = "two_cell";
cfg.vclamp.cells = [1, 2];
cfg.vclamp.Vrest = Vrest;
cfg.vclamp.t_rest = t_rest;
cfg.vclamp.V_step = V_step(:);
cfg.vclamp.t_step = t_step(:);
cfg.vclamp.t_tail = t_tail;
cfg.vclamp.release_after_step = release_after_step;
cfg.T = t_rest + max(t_step) + t_tail;
end

function save_voltage_clamp_result(result, name)
cfg = result.cfg;
if ~exist(cfg.save_data_path, "dir")
    mkdir(cfg.save_data_path);
end
save_file = fullfile(cfg.save_data_path, name + ".mat");
save(save_file, "result", "-v7.3");
fprintf("Saved %s\n", save_file);
end
