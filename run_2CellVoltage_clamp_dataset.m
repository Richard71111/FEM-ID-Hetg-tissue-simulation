%% RUN_2CELLVOLTAGE_CLAMP_DATASET
% Slurm-array version.
% One Slurm array task runs one ordinary MATLAB for-loop.
% No parpool, no parfor.

clear;
clc;
close all;

%% Slurm/task settings
GJ_coupling = string(getenv("GJ_coupling"));
if strlength(GJ_coupling) == 0
    GJ_coupling = "strong";
end
task_id = read_env_number("SLURM_ARRAY_TASK_ID", 1);
cases_per_task = read_env_number("CASES_PER_TASK", 1);
case_offset    = read_env_number("CASE_ID_OFFSET", 0); 
base_seed = 20260623;

case_ids = case_offset + (((task_id - 1) * cases_per_task + 1):(task_id * cases_per_task));

fprintf("[%s] MATLAB task started. SLURM_ARRAY_TASK_ID=%d, cases_per_task=%d, case range=%d-%d\n", ...
    string(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss")), ...
    task_id, cases_per_task, case_ids(1), case_ids(end));

%% Paths
project_folder = fileparts(mfilename("fullpath"));
addpath(fullfile(project_folder, "main_function"));
addpath(fullfile(project_folder, "config"));
addpath(fullfile(project_folder, "topology"));
addpath(fullfile(project_folder, "models"));
addpath(fullfile(project_folder, "matrix"));
addpath(fullfile(project_folder, "simulation"));

%% Base config
base_cfg = default_config();

base_cfg.adjacency_matrix = [
    0, 1;
    1, 0
];
base_cfg.cell_coordinates = (1:2)';
base_cfg.cell_port_count = 2;
base_cfg.junction_mesh = 1;

base_cfg.model = "ORd11";
base_cfg.stim_cell = 1;
base_cfg.stim_amp = 0;
base_cfg.stim_dur = 0;
base_cfg.nbeats = 1;
base_cfg.clamp_flag = false(4, 1);
base_cfg.show_progress = false;
base_cfg.make_plots = false;
base_cfg.save_every = 10;

base_cfg.adaptive_dt = false;
base_cfg.dt = 0.01;
base_cfg.dt2 = base_cfg.dt;
base_cfg.dtS = base_cfg.dt;
base_cfg.dtS2 = base_cfg.dt2;
switch GJ_coupling
    case "strong", base_cfg.D = 1.0;
    case "weak",   base_cfg.D = 0.1;
    otherwise, error("Unknown GJ_coupling: '%s'", GJ_coupling);
end
%% Protocol settings
Vrest = -0.879989146999539e2;
t_rest = 5;
t_tail = 5;
t_range = [1, 20];
v_range = [-70, 30];
release_after_step = false;

%% Save folder
save_folder_path = "/fs/ess/PAS1622/RichardSui/FML_ID_dataset/voltage_clamp_ds";
mesh_name = get_mesh_name(base_cfg.mesh_files);
target_folder = fullfile(save_folder_path, mesh_name, GJ_coupling);

ensure_folder(target_folder);

fprintf("[%s] Target folder: %s\n", ...
    string(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss")), target_folder);

%% Local temporary folder
local_dir = string(getenv("TMPDIR"));
if strlength(local_dir) == 0 || ~exist(local_dir, "dir")
    local_dir = string(tempdir);
end

fprintf("[%s] Local tmp folder: %s\n", ...
    string(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss")), local_dir);
fprintf("[%s] GJ_coupling=%s, D=%.3g\n", ...
    string(datetime("now","Format","yyyy-MM-dd HH:mm:ss")), GJ_coupling, base_cfg.D);

%% Main loop for this Slurm task
for local_case_idx = 1:cases_per_task
    case_id = case_ids(local_case_idx);
    rng(base_seed + case_id, "twister");

    t_step = round(sample_log_uniform(t_range, 2), 2);
    V_step = sample_linear_uniform(v_range, 2);

    file_name = generate_case_filename(base_cfg, case_id, t_step, V_step);
    local_save_path = fullfile(local_dir, file_name);
    final_save_path = fullfile(target_folder, file_name);

    fprintf("[%s] START task=%d local_case=%d/%d case=%d | t=[%.2f %.2f], V=[%.3f %.3f]\n", ...
        string(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss")), ...
        task_id, local_case_idx, cases_per_task, case_id, ...
        t_step(1), t_step(2), V_step(1), V_step(2));

    try
        cfg = make_two_cell_cfg( ...
            base_cfg, Vrest, t_rest, V_step, t_step, t_tail, release_after_step);

        case_tic = tic;
        result = run_graph_voltage_clamp_simulation(cfg);
        elapsed_sec = toc(case_tic);

        t = result.time;
        phi_axial = result.phi_axial;
        Icleft = result.Icleft;
        S_cleft = result.S_cleft;
        G_state = result.Gstate;

        save_data(local_save_path, ...
            t, phi_axial, Icleft, S_cleft, G_state, ...
            case_id, task_id, local_case_idx, t_step, V_step);

        if exist(final_save_path, "file")
            fprintf("[%s] SKIP  task=%d local_case=%d/%d case=%d | file exists: %s\n", ...
                string(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss")), ...
                task_id, local_case_idx, cases_per_task, case_id, final_save_path);
            if exist(local_save_path, "file"), delete(local_save_path); end
            continue;
        end

        [status, msg] = copyfile(local_save_path, final_save_path, 'f');
        if ~status
            error("Failed to copy %s to %s. Reason: %s", ...
                local_save_path, final_save_path, msg);
        end

        if exist(local_save_path, "file")
            delete(local_save_path);
        end

        fprintf("[%s] DONE  task=%d local_case=%d/%d case=%d | %.2f sec | %s\n", ...
            string(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss")), ...
            task_id, local_case_idx, cases_per_task, case_id, ...
            elapsed_sec, final_save_path);
        clear result t phi_axial Icleft S_cleft G_state t_step V_step
    catch ME
        fprintf("[%s] FAIL  task=%d local_case=%d/%d case=%d | %s\n", ...
            string(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss")), ...
            task_id, local_case_idx, cases_per_task, case_id, ME.message);
        continue;
    end
end

fprintf("[%s] MATLAB task finished. task=%d\n", ...
    string(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss")), task_id);

%% Helper functions
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

function filename = generate_case_filename(cfg, case_id, t_step, V_step)
mesh_name = get_mesh_name(cfg.mesh_files);

filename = sprintf("%s_case_%06d_t%s_%s_v%s_%s.mat", ...
    char(mesh_name), ...
    case_id, ...
    number_tag(t_step(1)), ...
    number_tag(t_step(2)), ...
    number_tag(V_step(1)), ...
    number_tag(V_step(2)));
end

function tag = number_tag(value)
tag = sprintf("%.6g", value);
tag = strrep(tag, "-", "m");
tag = strrep(tag, "+", "");
tag = strrep(tag, ".", "p");
end

function mesh_name = get_mesh_name(mesh_files)
mesh_files = string(mesh_files);
mesh_file = mesh_files(1);

[~, stem] = fileparts(char(mesh_file));
token = regexp(stem, '^FEMDATA_(\d+)$', 'tokens', 'once');

if ~isempty(token)
    mesh_name = "FEM" + string(token{1});
else
    mesh_name = string(regexprep(stem, '[^A-Za-z0-9_]+', '_'));
end
end

function save_data(save_path, ...
    t, phi_axial, Icleft, S_cleft, G_state, ...
    case_id, task_id, local_case_idx, t_step, V_step)

save(save_path, ...
    "t", "phi_axial", "Icleft", "S_cleft", "G_state", ...
    "case_id", "task_id", "local_case_idx", "t_step", "V_step", ...
    "-v7");

fprintf("[%s] Saved local file: %s\n", ...
    string(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss")), save_path);
end

function value = read_env_number(name, default_value)
text = getenv(char(name));
if isempty(text)
    value = default_value;
    return;
end

value = str2double(text);
if isnan(value)
    value = default_value;
end
value = floor(value);
end

function ensure_folder(folder_path)
if ~exist(folder_path, "dir")
    [status, msg, msgID] = mkdir(folder_path);
    if ~status
        error(msgID, "Failed to create folder: %s\nReason: %s", folder_path, msg);
    end
end
end