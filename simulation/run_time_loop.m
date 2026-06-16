function result = run_time_loop(cfg, topology, model, network)
%RUN_TIME_LOOP Advance ionic states, cleft concentrations, and potentials.
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
%   Gstate      (Nstate*Npatches)-by-Nt ionic state vector. The ORd11/LR1/
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

p = model.p;
p.BCL = cfg.BCL;
p.Npatches = Npatches;
p.f_I = network.f_I;
p.indstim = zeros(Npatches, 1);
p.indstim(patch.axial(cfg.stim_cell)) = 1;
if lower(cfg.model) == "ord11"
    p.celltype = zeros(Npatches, 1);
end

%% Initial state tensors
phi = model.x0(1) * ones(Nnodes, 1);
phi(layout.cleft(:)) = 0;

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
solver_fine = decomposition(network.G + network.Cm / dt_fine, "lu");
if use_adaptive
    dt_coarse = cfg.dt2;
    Nsplit_coarse = max(1, round(cfg.dt2 / cfg.dtS2));
    solver_coarse = decomposition(network.G + network.Cm / dt_coarse, "lu");
else
    dt_coarse = dt_fine;
    Nsplit_coarse = Nsplit_fine;
    solver_coarse = solver_fine;
end

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

NG = model.Nstate * Npatches;
time = nan(1, max_samples);
phi_axial = nan(Ncell, max_samples);
Gstate_hist = nan(NG, max_samples);
Icleft = nan(2, Njunction, max_samples);          % Row 1 = pre-side, row 2 = post-side.
S_cleft_hist = nan(4, M, Njunction, max_samples);

%% Initial sample (count = 1, t = 0)
count = 1;
time(count) = 0;
phi_axial(:, count) = phi(cell_index);
Gstate_hist(:, count) = Gstate;
Icleft(:, :, count) = compute_icleft(phi, ...
    icleft_cell, icleft_ID_idx, gmyo, M, Njunction);
if Njunction > 0
    S_cleft_hist(:, :, :, count) = S_cleft;
end

ti = 0;
step = 0;
Iall = [];

%% Main simulation loop
while ti < cfg.T
    % Adaptive step selection: fine step inside the post-beat window
    % (mod(t, BCL) < twin), coarse step otherwise.
    if use_adaptive && mod(ti, cfg.BCL) >= cfg.twin
        dt = dt_coarse;
        solver = solver_coarse;
        Nsplit = Nsplit_coarse;
    else
        dt = dt_fine;
        solver = solver_fine;
        Nsplit = Nsplit_fine;
    end
    % Clip the final step to land exactly on T (refactor once for this step).
    if dt > cfg.T - ti
        dt = cfg.T - ti;
        solver = decomposition(network.G + network.Cm / dt, "lu");
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
    phi_new = solver \ rhs;

    if any(~isfinite(phi_new)) || any(~isfinite(Gnew))
        error("Non-finite state at t = %.6f ms.", ti);
    end

    phi = phi_new;
    Gstate = Gnew;
    ti = round(ti + dt, 5);
    step = step + 1;

    % Cross-step saving: store every save_every accepted steps and the final.
    if mod(step, save_every) == 0 || ti >= cfg.T
        count = count + 1;
        if count > size(time, 2)
            % Grow buffers (robust to adaptive-step runs with extra samples).
            new_size = 2 * size(time, 2);
            time(1, new_size) = nan;
            phi_axial(Ncell, new_size) = nan;
            Gstate_hist(NG, new_size) = nan;
            Icleft(2, max(Njunction, 1), new_size) = nan;
            S_cleft_hist(4, M, max(Njunction, 1), new_size) = nan;
        end
        time(count) = ti;
        phi_axial(:, count) = phi(cell_index);
        Gstate_hist(:, count) = Gstate;
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
result.Gstate = Gstate_hist(:, 1:count);           % Ionic states (incl. intracellular conc.).
result.Icleft = Icleft(:, :, 1:count);             % Cleft current per junction side.
result.S_cleft = S_cleft_hist(:, :, :, 1:count);   % Cleft concentration, 4-by-M-by-Njunction-by-Nt.
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
