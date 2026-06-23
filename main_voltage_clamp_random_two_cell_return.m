%% MAIN_VOLTAGE_CLAMP_RANDOM_TWO_CELL_RETURN
% Protocol family 1: force the voltage back to Vrest after the pulse.
% Both cells are voltage clamped. t1/t2 are sampled log-uniformly from
% [1, 20] ms, and v1/v2 are sampled linearly from [-70, 30] mV.

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
base_cfg.save_every = 1;

% Use a fixed time step so the random pulse durations are resolved exactly.
base_cfg.adaptive_dt = false;
base_cfg.dt = 0.01;
base_cfg.dt2 = base_cfg.dt;
base_cfg.dtS = base_cfg.dt;
base_cfg.dtS2 = base_cfg.dt2;

%% Random protocol settings
rng('shuffle');

num_cases = 1;               % Increase this for a random batch.
make_plots = true;

Vrest = -0.879989146999539e2;    % mV.
t_rest = 5;                     % ms at Vrest before the pulse.
t_tail = 5;                     % ms held at Vrest after the longest pulse.

t_range = [1, 20];               % ms, sampled on log scale.
v_range = [-70, 30];             % mV, sampled on linear scale.

case_id = (1:num_cases)';
t1_ms = nan(num_cases, 1);
t2_ms = nan(num_cases, 1);
v1_mV = nan(num_cases, 1);
v2_mV = nan(num_cases, 1);
save_name = strings(num_cases, 1);

batch_tag = "family1_random_two_cell_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss"));

for k = 1:num_cases
    t_step = round(sample_log_uniform(t_range, 2),2);
    V_step = sample_linear_uniform(v_range, 2);

    cfg = make_two_cell_return_cfg(base_cfg, Vrest, t_rest, V_step, t_step, t_tail);

    t1_ms(k) = t_step(1);
    t2_ms(k) = t_step(2);
    v1_mV(k) = V_step(1);
    v2_mV(k) = V_step(2);
    save_name(k) = sprintf("%s_case_%03d", batch_tag, k);

    fprintf("\nCase %d/%d\n", k, num_cases);
    fprintf("Cell 1: t1 = %.6g ms, v1 = %.6g mV\n", t_step(1), V_step(1));
    fprintf("Cell 2: t2 = %.6g ms, v2 = %.6g mV\n", t_step(2), V_step(2));

    tic
    result = run_graph_voltage_clamp_simulation(cfg);
    toc

    result.random_protocol = struct( ...
        "t1_ms", t_step(1), ...
        "t2_ms", t_step(2), ...
        "v1_mV", V_step(1), ...
        "v2_mV", V_step(2), ...
        "t_sampling", "log-uniform", ...
        "v_sampling", "linear-uniform");

    % save_voltage_clamp_result(result, save_name(k));

    if make_plots
        figure_name = sprintf( ...
            "Voltage clamp return to rest random two cells case %03d", k);
        plot_voltage_clamp_result(result, figure_name);
    end
end

summary = table(case_id, t1_ms, t2_ms, v1_mV, v2_mV, save_name);
summary_file = fullfile(base_cfg.save_data_path, batch_tag + "_summary.csv");
% writetable(summary, summary_file);
fprintf("Saved %s\n", summary_file);

function cfg = make_two_cell_return_cfg( ...
    base_cfg, Vrest, t_rest, V_step, t_step, t_tail)
cfg = base_cfg;
cfg.vclamp = struct();
cfg.vclamp.mode = "two_cell";
cfg.vclamp.cells = [1; 2];
cfg.vclamp.Vrest = Vrest;
cfg.vclamp.t_rest = t_rest;
cfg.vclamp.V_step = V_step(:);
cfg.vclamp.t_step = t_step(:);
cfg.vclamp.t_tail = t_tail;
cfg.vclamp.release_after_step = false;
cfg.T = t_rest + max(t_step) + t_tail;
end

function values = sample_log_uniform(range, n)
lo = range(1);
hi = range(2);
values = exp(log(lo) + (log(hi) - log(lo)) * rand(n, 1));
end

function values = sample_linear_uniform(range, n)
lo = range(1);
hi = range(2);
values = lo + (hi - lo) * rand(n, 1);
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