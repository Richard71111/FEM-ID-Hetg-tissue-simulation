function result = run_time_loop_voltage_clamp(cfg, topology, model, network)
%RUN_TIME_LOOP_VOLTAGE_CLAMP Advance the model with axial-cell voltage clamp.
% Inputs: cfg, topology, ionic model, and assembled graph network.
% Output: sampled time histories of the retained physical variables.
%
% Adaptive (dual) time step (cfg.adaptive_dt): the fine step cfg.dt is used
% within the first cfg.twin ms after each beat onset and the coarse step
% cfg.dt2 for the rest of the cycle, matching the original 1-D source.
%
% Cross-step saving (cfg.save_every): outputs are stored every N accepted
% time steps (plus the initial state and the final step), instead of on a
% fixed time interval. This mirrors the step-based saving of the original
% 1-D source code.
%
% Saved time histories (everything else is dropped):
%   time        1-by-Nt sample times, ms.
%   phi_axial   Ncell-by-Nt axial (intracellular) node potential, mV.
%   Gstate      Nstate-by-Ncell-by-Nt axial ionic state tensor. The ORd11/LR1/
%               Court98 intracellular concentrations live inside this state,
%               so intracellular concentration is NOT saved separately.
%   Icleft      2-by-Njunction-by-Nt cleft (axial-to-ID) current per side.
%   S_cleft     4-by-M-by-Njunction-by-Nt cleft (extracellular) concentration.

Ncell = topology.Ncell;
Njunction = topology.Njunction;
M = network.M;
Nnodes = network.Nnodes;
Npatches = network.Npatches;
layout = network.layout;
patch = network.patch;
mesh = network.mesh;
protocol = normalize_vclamp_protocol(cfg.vclamp, Ncell, cfg.T);
all_nodes = (1:Nnodes)';

p = model.p;
p.BCL = cfg.BCL;
p.Npatches = Npatches;
p.f_I = network.f_I;
p.indstim = zeros(Npatches, 1);
if cfg.stim_amp ~= 0 && cfg.stim_dur > 0
    p.indstim(patch.axial(cfg.stim_cell)) = 1;
end
if lower(cfg.model) == "ord11"
    p.celltype = zeros(Npatches, 1);
end

%% Initial state tensors
phi = model.x0(1) * ones(Nnodes, 1);
phi(layout.cleft(:)) = 0;
[active_cells0, command0] = voltage_clamp_command(0, protocol);
phi(layout.cell(active_cells0)) = command0;

bulk = [cfg.Na_b; cfg.K_b; cfg.Ca_b; cfg.A_b];
S_cleft = repmat(reshape(bulk, 4, 1, 1), 1, M, Njunction);
z = [1; 1; 2; -1];

Gstate = zeros(model.Nstate * Npatches, 1);
for state_number = 1:model.Nstate
    rows = (state_number - 1) * Npatches + (1:Npatches);
    Gstate(rows) = model.x0(state_number + 1);
end

RTF = cfg.R * cfg.Temp / cfg.F;

% Time-step setup (adaptive dual step, matching the original 1-D source).
% Within the first cfg.twin ms after each beat onset the solver uses the fine
% step cfg.dt; for the rest of the cycle it uses the coarse step cfg.dt2.
% The coefficient factorization for each step is precomputed once and reused,
% so switching steps does not refactor the system every iteration.
use_adaptive = isfield(cfg, "adaptive_dt") && cfg.adaptive_dt;
dt_fine = cfg.dt;
Nsplit_fine = max(1, round(cfg.dt / cfg.dtS));
if use_adaptive
    dt_coarse = cfg.dt2;
    Nsplit_coarse = max(1, round(cfg.dt2 / cfg.dtS2));
else
    dt_coarse = dt_fine;
    Nsplit_coarse = Nsplit_fine;
end
solver_cache = struct("key", {}, "solver", {});

%% Axial-to-ID coupling used for the cleft current Icleft
% Each junction has two disc faces (pre- and post-junctional), each carrying
% its own axial coupling current into the ID patches of the connected cell:
%   Icleft(side, junction) = (gmyo/M) * sum_{j=1..M} (phi_cell - phi_ID_{j,side})
% Columns of icleft_cell / icleft_ID_idx run junction-major, side-minor so
% that reshape(.,2,Njunction) maps row = side, column = junction.
gmyo = network.gmyo;
cell_index = layout.cell(:);            % Ncell-by-1 axial node indices.
side_count = 2 * Njunction;
icleft_cell = zeros(1, side_count);     % Owner cell of each junction side.
icleft_ID_idx = zeros(M, max(side_count, 1));  % ID node indices per side.
col = 0;
for j = 1:Njunction
    for side = 1:2
        col = col + 1;
        icleft_cell(col) = topology.junction_cells(side, j);
        icleft_ID_idx(:, col) = layout.ID(:, side, j);
    end
end

%% Sampled output buffers (preallocated, grown on demand)
% Step-based saving: store every cfg.save_every accepted steps. With a fixed
% step there are about ceil(T/dt) steps; allocate from that and grow the
% buffers if an adaptive step ever produces more samples than expected.
save_every = max(1, round(cfg.save_every));
if use_adaptive
    % Fine steps within twin plus coarse steps for the rest of each cycle.
    per_beat = cfg.twin / dt_fine + max(cfg.BCL - cfg.twin, 0) / dt_coarse;
    nsteps_guess = max(1, ceil(per_beat * (cfg.T / cfg.BCL)));
else
    nsteps_guess = max(1, ceil(cfg.T / max(dt_fine, eps)));
end
max_samples = ceil(nsteps_guess / save_every) + 2;

Nstate = model.Nstate;
axial_patches = patch.axial(:);
time = nan(1, max_samples);
phi_axial = nan(Ncell, max_samples);
Gstate_hist = nan(Nstate, Ncell, max_samples);
Icleft = nan(2, Njunction, max_samples);          % Row 1 = pre-side, row 2 = post-side.
S_cleft_hist = nan(4, M, Njunction, max_samples);
vclamp_command = nan(Ncell, max_samples);

%% Initial sample (count = 1, t = 0)
count = 1;
time(count) = 0;
phi_axial(:, count) = phi(cell_index);
Gstate_hist(:, :, count) = extract_axial_Gstate( ...
    Gstate, Npatches, Nstate, axial_patches);
vclamp_command(active_cells0, count) = command0;
Icleft(:, :, count) = compute_icleft(phi, ...
    icleft_cell, icleft_ID_idx, gmyo, M, Njunction);
if Njunction > 0
    S_cleft_hist(:, :, :, count) = S_cleft;
end

ti = 0;
step = 0;

%% Main simulation loop
while ti < cfg.T
    % Adaptive step selection: fine step inside the post-beat window
    % (mod(t, BCL) < twin), coarse step otherwise.
    if use_adaptive && mod(ti, cfg.BCL) >= cfg.twin
        dt = dt_coarse;
        Nsplit = Nsplit_coarse;
    else
        dt = dt_fine;
        Nsplit = Nsplit_fine;
    end
    % Clip steps to land exactly on final time and voltage-command edges.
    next_event = next_voltage_clamp_event(ti, protocol, cfg.T);
    if dt > next_event - ti
        dt = next_event - ti;
    end
    p.dt = dt;

    Vm = network.Am * phi;
    Sp = repmat(bulk(1:3)', Npatches, 1);
    ID_patches = patch.ID(:);
    cleft_index = patch.cleft_linear(ID_patches);
    for ion = 1:3
        values = reshape(S_cleft(ion, :, :), M * Njunction, 1);
        Sp(ID_patches, ion) = values(cleft_index);
    end

    [Gnew, Iion, Ivec, ~, ~] = ...
        model.ionic_fun(ti, [Vm; Gstate], p, Sp(:));
    Ivec = reshape(Ivec, Npatches, 3);

    Esource = zeros(M, Njunction);
    if ~all(cfg.clamp_flag)
        Idisc = zeros(4, M, Njunction);
        for j = 1:Njunction
            for side = 1:2
                ID_patches = patch.ID(:, side, j);
                for ion = 1:3
                    Idisc(ion, :, j) = Idisc(ion, :, j) + ...
                        reshape(Ivec(ID_patches, ion), 1, M);
                end
            end
        end

        dtS = dt / Nsplit;
        for split_step = 1:Nsplit
            Esource(:) = 0;
            for j = 1:Njunction
                phi_cleft = phi(layout.cleft(:, j));
                voltage_difference = phi_cleft - phi_cleft';
                Gc_species = mesh.Gc(:, :, j) / 4;
                Gb_species = mesh.Gb(:, j) / 4;

                for ion = 1:4
                    concentration = reshape(S_cleft(ion, :, j), M, 1);
                    Erev = RTF / z(ion) * ...
                        log(concentration' ./ concentration);
                    Ebulk = RTF / z(ion) * ...
                        log(bulk(ion) ./ concentration);

                    Iout = sum(Gc_species .* ...
                        (voltage_difference - Erev), 2) + ...
                        Gb_species .* (phi_cleft - Ebulk);
                    Esource(:, j) = Esource(:, j) + ...
                        sum(Gc_species .* Erev, 2) + ...
                        Gb_species .* Ebulk;

                    if ~cfg.clamp_flag(ion)
                        dS = dtS * ...
                            (reshape(Idisc(ion, :, j), M, 1) - Iout) * ...
                            1e6 ./ (z(ion) * cfg.F * mesh.volume(:, j));
                        S_cleft(ion, :, j) = ...
                            reshape(max(concentration + dS, 1e-9), 1, M);
                    end
                end
            end
        end
    end

    source = zeros(Nnodes, 1);
    source(layout.cleft(:)) = -Esource(:);
    rhs = network.Cm / dt * phi - network.Am' * Iion + source;
    [active_cells, command] = voltage_clamp_command(ti + dt, protocol);
    active_nodes = layout.cell(active_cells);
    [solver, solver_cache] = get_vclamp_solver( ...
        solver_cache, network, dt, all_nodes, active_nodes);
    phi_new = solve_vclamp_step(solver, rhs, command);

    if any(~isfinite(phi_new)) || any(~isfinite(Gnew))
        error("Non-finite state at t = %.6f ms.", ti);
    end

    phi = phi_new;
    Gstate = Gnew;
    ti = ti + dt;
    step = step + 1;
    if step > 4 * max_samples * save_every
        error("Step runaway at t=%.6f (T=%.6f).", ti, cfg.T);
    end

    % Cross-step saving: store every save_every accepted steps and the final.
    if mod(step, save_every) == 0 || ti >= cfg.T
        count = count + 1;
        if count > size(time, 2)
            % Grow buffers (robust to adaptive-step runs with extra samples).
            new_size = 2 * size(time, 2);
            time(1, new_size) = nan;
            phi_axial(Ncell, new_size) = nan;
            Gstate_hist(Nstate, Ncell, new_size) = nan;
            Icleft(2, max(Njunction, 1), new_size) = nan;
            S_cleft_hist(4, M, max(Njunction, 1), new_size) = nan;
            vclamp_command(Ncell, new_size) = nan;
        end
        time(count) = ti;
        phi_axial(:, count) = phi(cell_index);
        Gstate_hist(:, :, count) = extract_axial_Gstate( ...
            Gstate, Npatches, Nstate, axial_patches);
        vclamp_command(active_cells, count) = command;
        Icleft(:, :, count) = compute_icleft(phi, ...
            icleft_cell, icleft_ID_idx, gmyo, M, Njunction);
        if Njunction > 0
            S_cleft_hist(:, :, :, count) = S_cleft;
        end
    end

    if cfg.show_progress && ...
            abs(ti / cfg.progress_interval - ...
            round(ti / cfg.progress_interval)) < 1e-8
        fprintf("t = %.3f / %.3f ms\n", ti, cfg.T);
    end
end

%% Retained time histories
result.time = time(1:count);
result.phi_axial = phi_axial(:, 1:count);          % Axial node potential, Ncell-by-Nt.
result.Gstate = Gstate_hist(:, :, 1:count);       % Axial ionic states, Nstate-by-Ncell-by-Nt.
result.Icleft = Icleft(:, :, 1:count);             % Cleft current per junction side.
result.S_cleft = S_cleft_hist(:, :, :, 1:count);   % Cleft concentration, 4-by-M-by-Njunction-by-Nt.
result.vclamp = struct();
result.vclamp.protocol = protocol;
result.vclamp.command = vclamp_command(:, 1:count);
end

function protocol = normalize_vclamp_protocol(vclamp, Ncell, T)
%NORMALIZE_VCLAMP_PROTOCOL Convert user-facing protocol fields to vectors.
if isfield(vclamp, "mode")
    mode = lower(string(vclamp.mode));
else
    mode = "one_cell";
end

if isfield(vclamp, "cells") && ~isempty(vclamp.cells)
    cells = vclamp.cells(:);
elseif mode == "two_cell" || mode == "both_cells"
    cells = [1; 2];
else
    cells = 1;
end

if any(cells < 1) || any(cells > Ncell) || any(cells ~= round(cells))
    error("cfg.vclamp.cells must contain valid cell indices.");
end
if numel(unique(cells)) ~= numel(cells)
    error("cfg.vclamp.cells must not contain duplicate cells.");
end

Vrest = get_protocol_field(vclamp, "Vrest", -87);
t_rest = get_protocol_field(vclamp, "t_rest", 100);

if isfield(vclamp, "V_step")
    V_step = vclamp.V_step(:);
else
    V1 = get_protocol_field(vclamp, "V1", -20);
    V2 = get_protocol_field(vclamp, "V2", V1);
    V_step = [V1; V2];
end

if isfield(vclamp, "t_step")
    t_step = vclamp.t_step(:);
else
    t1 = get_protocol_field(vclamp, "t1", 200);
    t2 = get_protocol_field(vclamp, "t2", t1);
    t_step = [t1; t2];
end

if isscalar(V_step)
    V_step = repmat(V_step, numel(cells), 1);
end
if isscalar(t_step)
    t_step = repmat(t_step, numel(cells), 1);
end
if numel(V_step) < numel(cells) || numel(t_step) < numel(cells)
    error("Voltage-clamp V_step and t_step must match cfg.vclamp.cells.");
end
V_step = V_step(1:numel(cells));
t_step = t_step(1:numel(cells));

if t_rest < 0 || any(t_step < 0)
    error("Voltage-clamp durations must be nonnegative.");
end

protocol = struct();
protocol.mode = mode;
protocol.cells = cells;
protocol.Vrest = Vrest;
protocol.t_rest = t_rest;
protocol.V_step = V_step;
protocol.t_step = t_step;
protocol.release_after_step = isfield(vclamp, "release_after_step") && ...
    logical(vclamp.release_after_step);
protocol.events = unique([t_rest; t_rest + t_step(:); T]);
protocol.events = protocol.events(protocol.events >= 0 & protocol.events <= T);
end

function value = get_protocol_field(s, name, default_value)
if isfield(s, name)
    value = s.(name);
else
    value = default_value;
end
end

function [active_cells, command] = voltage_clamp_command(t, protocol)
%VOLTAGE_CLAMP_COMMAND Return the cells actively held at the current time.
active_cells = zeros(0, 1);
command = zeros(0, 1);
for k = 1:numel(protocol.cells)
    step_end = protocol.t_rest + protocol.t_step(k);
    if t < protocol.t_rest
        active = true;
        voltage = protocol.Vrest;
    elseif t <= step_end
        active = true;
        voltage = protocol.V_step(k);
    elseif protocol.release_after_step
        active = false;
        voltage = nan;
    else
        active = true;
        voltage = protocol.Vrest;
    end

    if active
        active_cells(end + 1, 1) = protocol.cells(k); %#ok<AGROW>
        command(end + 1, 1) = voltage; %#ok<AGROW>
    end
end
end

function next_event = next_voltage_clamp_event(t, protocol, T)
tol = 1e-10;
events = protocol.events(protocol.events > t + tol);
if isempty(events)
    next_event = T;
else
    next_event = min(events);
end
if next_event <= t + tol
    next_event = T;
end
end

function [solver, solver_cache] = get_vclamp_solver( ...
    solver_cache, network, dt, all_nodes, active_nodes)
key = solver_cache_key(dt, active_nodes);
if isempty(solver_cache)
    cache_index = [];
else
    cache_index = find(strcmp({solver_cache.key}, key), 1);
end

if isempty(cache_index)
    solver = make_vclamp_solver(network, dt, all_nodes, active_nodes);
    solver_cache(end + 1).key = key;
    solver_cache(end).solver = solver;
else
    solver = solver_cache(cache_index).solver;
end
end

function key = solver_cache_key(dt, active_nodes)
node_text = strjoin(string(active_nodes(:).'), ",");
key = char("dt=" + string(sprintf("%.12g", dt)) + ";nodes=" + node_text);
end

function solver = make_vclamp_solver(network, dt, all_nodes, active_nodes)
A = network.G + network.Cm / dt;
solver.A = A;
solver.clamped_nodes = active_nodes(:);
solver.free_nodes = setdiff(all_nodes(:), solver.clamped_nodes);
solver.Afc = A(solver.free_nodes, solver.clamped_nodes);
solver.free_solver = decomposition(A(solver.free_nodes, solver.free_nodes), "lu");
end

function phi_new = solve_vclamp_step(solver, rhs, command)
%SOLVE_VCLAMP_STEP Apply Dirichlet axial-cell voltages.
phi_new = zeros(size(rhs));
if isempty(solver.clamped_nodes)
    phi_new(solver.free_nodes) = solver.free_solver \ rhs(solver.free_nodes);
else
    phi_new(solver.clamped_nodes) = command(:);
    phi_new(solver.free_nodes) = solver.free_solver \ ...
        (rhs(solver.free_nodes) - solver.Afc * command(:));
end
end

function Ic = compute_icleft(phi, icleft_cell, icleft_ID_idx, gmyo, M, Njunction)
%COMPUTE_ICLEFT Axial-to-ID (cleft) current for each junction side.
% Every junction has two disc faces, so two currents:
%   Ic(1, junction) = (gmyo/M) * sum_{j=1..M} (phi_preCell  - phi_ID_pre,j )
%   Ic(2, junction) = (gmyo/M) * sum_{j=1..M} (phi_postCell - phi_ID_post,j)
% Row 1 is the pre-junctional side (first cell of the edge), row 2 the
% post-junctional side (second cell). icleft_cell / icleft_ID_idx columns
% are ordered (junction-major, side-minor) to match reshape(.,2,Njunction).
Ic = zeros(2, Njunction);
if Njunction == 0
    return;
end
pc = phi(icleft_cell(:)).';                      % 1-by-side_count axial potential of each side.
sum_phi_ID = sum(phi(icleft_ID_idx), 1);         % 1-by-side_count sum of ID potentials.
contrib = (gmyo / M) * (M * pc - sum_phi_ID);    % 1-by-side_count.
Ic = reshape(contrib, 2, Njunction);             % Row = side, column = junction.
end

function axial_Gstate = extract_axial_Gstate(Gstate, Npatches, Nstate, axial_patches)
%EXTRACT_AXIAL_GSTATE Return axial ionic states as Nstate-by-Ncell.
Gstate_by_patch = reshape(Gstate, Npatches, Nstate);
axial_Gstate = Gstate_by_patch(axial_patches, :).';
end
