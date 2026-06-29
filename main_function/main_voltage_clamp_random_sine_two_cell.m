function main_voltage_clamp_random_sine_two_cell(num_cases)
%MAIN_VOLTAGE_CLAMP_RANDOM_SINE_TWO_CELL Generate random sine-clamp data.
%
% Usage:
%   main_voltage_clamp_random_sine_two_cell        % one case + plot
%   main_voltage_clamp_random_sine_two_cell(10000) % batch, no plots

if nargin < 1
    num_cases = 1;
end
validateattributes(num_cases, {'numeric'}, ...
    {'scalar', 'integer', 'positive'});

function_folder = fileparts(mfilename("fullpath"));
project_folder  = fileparts(function_folder);
addpath(fullfile(project_folder, "main_function"));
addpath(fullfile(project_folder, "plot_function"));
addpath(fullfile(project_folder, "config"));
addpath(fullfile(project_folder, "topology"));
addpath(fullfile(project_folder, "models"));
addpath(fullfile(project_folder, "matrix"));
addpath(fullfile(project_folder, "simulation"));

%% Fixed generation settings
dt = 0.01;                 % ms
T = 50;                    % ms
Vrest = -0.879989146999539e2; % mV
ramp_time = 2;             % ms, avoids an artificial t=0 current spike
center_range = [-60, 50]; % mV
amplitude_range = [5, 20]; % mV
period_range = [0.5, 10];  % ms

cfg_base = make_base_config(dt);
command_time = 0:dt:T;
output_dir = fullfile(cfg_base.save_data_path, ...
    "random_sine_vclamp_dataset");
if ~exist(output_dir, "dir")
    mkdir(output_dir);
end

rng("shuffle");
batch_tag = "random_sine_vclamp_" ...
    + string(datetime("now", "Format", "yyyyMMdd_HHmmss"));

for case_idx = 1:num_cases
    case_seed = randi(2^31 - 1);
    stream = RandStream("mt19937ar", "Seed", case_seed);

    center = sample_uniform(stream, center_range, 2);
    amplitude = sample_uniform(stream, amplitude_range, 2);
    period = sample_log_uniform(stream, period_range, 2);
    phase = 2 * pi * rand(stream, 2, 1);

    sine_voltage = center ...
        + amplitude .* sin(2 * pi * command_time ./ period + phase);
    ramp_fraction = min(command_time / ramp_time, 1);
    envelope = 0.5 * (1 - cos(pi * ramp_fraction));
    command_voltage = Vrest + envelope .* (sine_voltage - Vrest);

    cfg = cfg_base;
    cfg.T = T;
    cfg.vclamp = struct( ...
        "mode", "waveform", ...
        "cells", [1; 2], ...
        "command_time", command_time, ...
        "command_voltage", command_voltage, ...
        "interpolation", "linear", ...
        "figure_name", sprintf("Random sine clamp case %03d", case_idx));

    fprintf("\nCase %d/%d\n", case_idx, num_cases);
    fprintf("  center    = [%.3f, %.3f] mV\n", center);
    fprintf("  amplitude = [%.3f, %.3f] mV\n", amplitude);
    fprintf("  period    = [%.3f, %.3f] ms\n", period);
    tic
    result = run_graph_voltage_clamp_dynamic_simulation(cfg);
    toc
    t = result.time;
    phi_axial = result.phi_axial;
    Icleft = result.Icleft;
    vclamp_command = result.vclamp.command;
    protocol_metadata = struct( ...
        "seed", case_seed, ...
        "center_mV", center, ...
        "amplitude_mV", amplitude, ...
        "period_ms", period, ...
        "phase_rad", phase);

    % save_file = fullfile(output_dir, sprintf( ...
    %     "%s_case_%05d.mat", batch_tag, case_idx - 1));
    % save(save_file, "t", "phi_axial", "Icleft", ...
    %     "vclamp_command", "protocol_metadata", "-v7");
    % fprintf("  Saved %s\n", save_file);

    if num_cases == 1
        plot_voltage_clamp_result(result, cfg.vclamp.figure_name);
    end
end
end

function cfg = make_base_config(dt)
cfg = default_config();
cfg.adjacency_matrix = [0, 1; 1, 0];
cfg.cell_coordinates = (1:2)';
cfg.cell_port_count = 2;
cfg.junction_mesh = 1;
cfg.model = "ORd11";
cfg.stim_cell = 1;
cfg.stim_amp = 0;
cfg.stim_dur = 0;
cfg.nbeats = 1;
cfg.clamp_flag = false(4, 1);
cfg.show_progress = false;
cfg.make_plots = false;
cfg.save_every = 1;
cfg.adaptive_dt = false;
cfg.dt = dt;
cfg.dt2 = dt;
cfg.dtS = dt;
cfg.dtS2 = dt;
end

function values = sample_uniform(stream, range, n)
values = range(1) + (range(2) - range(1)) * rand(stream, n, 1);
end

function values = sample_log_uniform(stream, range, n)
values = exp(log(range(1)) ...
    + (log(range(2)) - log(range(1))) * rand(stream, n, 1));
end
