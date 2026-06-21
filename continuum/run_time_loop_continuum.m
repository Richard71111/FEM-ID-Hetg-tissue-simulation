function result = run_time_loop_continuum(cfg, topology, model, network)
%RUN_TIME_LOOP_CONTINUUM Advance the ID-cleft continuum notation variables.
%
% This file is the "equations and time stepping" layer. It deliberately uses
% the same names as the continuum note for the main physical fields:
%
% Primary variables match the derivation:
%   phi_m      Ncell-by-1 axial/myoplasmic potential.
%   phi_pre    M-by-Njunction pre-junctional ID potential.
%   phi_post   M-by-Njunction post-junctional ID potential.
%   phi_c      M-by-Njunction cleft potential.
%   s_c        4-by-M-by-Njunction cleft concentrations.
%
% The helper functions at the bottom do the unavoidable conversions between
% these readable fields and the sparse operators assembled in network:
%   stack_phi_continuum      readable fields -> global node vector
%   unstack_phi_continuum    global node vector -> readable fields
%   stack_membrane_voltage   note variables -> ionic-model Vm vector
%
% The electrical equation is kept in mass-matrix form:
%   Cm * dphi/dt = -G*phi - Am'*Iion + source.
% That is why the code builds rhs_phi first, then passes it to a selectable
% integrator in advance_phi_continuum.

Ncell = topology.Ncell;
Njunction = topology.Njunction;
M = network.M;
Nnodes = network.Nnodes;
Npatches = network.Npatches;
Nports = network.Nports;
patch = network.patch;
mesh = network.mesh;

p = model.p;
p.BCL = cfg.BCL;
p.Npatches = Npatches;
p.f_I = network.f_I;

% Stimulus and current-localization data are passed into the ionic model in
% the same patch order used by network.patch.
p.indstim = zeros(Npatches, 1);
p.indstim(patch.axial(cfg.stim_cell)) = 1;
if lower(cfg.model) == "ord11"
    p.celltype = zeros(Npatches, 1);
end

%% Initial state in continuum-note notation
phi_m = model.x0(1) * ones(Ncell, 1);
phi_pre = model.x0(1) * ones(M, Njunction);
phi_post = model.x0(1) * ones(M, Njunction);
phi_boundary = model.x0(1) * ones(Nports, Ncell);
phi_c = zeros(M, Njunction);

% s_b is the fixed bulk reservoir in the note. s_c is the cleft
% concentration field for Na, K, Ca, and the anion A.
s_b = [cfg.Na_b; cfg.K_b; cfg.Ca_b; cfg.A_b];
s_c = repmat(reshape(s_b, 4, 1, 1), 1, M, Njunction);
z_s = [1; 1; 2; -1];
alpha_s = cfg.R * cfg.Temp ./ (z_s * cfg.F);

% Gstate stores all non-voltage ionic gating/concentration states. It is not
% a continuum-note electrical field, but it belongs to the ionic model RHS.
Gstate = zeros(model.Nstate * Npatches, 1);
for state_number = 1:model.Nstate
    rows = (state_number - 1) * Npatches + (1:Npatches);
    Gstate(rows) = model.x0(state_number + 1);
end

%% Time-step setup
% Only the step size policy is chosen here. The electrical update method is
% selected later by cfg.electric_integrator inside advance_phi_continuum.
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

%% Sampled output buffers
save_every = max(1, round(cfg.save_every));
if use_adaptive
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
Icleft = nan(2, Njunction, max_samples);
S_cleft_hist = nan(4, M, Njunction, max_samples);

count = 1;
time(count) = 0;
phi_axial(:, count) = phi_m;
Gstate_hist(:, count) = Gstate;
Icleft(:, :, count) = compute_icleft_continuum( ...
    phi_m, phi_pre, phi_post, topology, network.gmyo, M);
if Njunction > 0
    S_cleft_hist(:, :, :, count) = s_c;
end

ti = 0;
step = 0;

%% Main simulation loop
while ti < cfg.T
    % Adaptive step choice: small voltage/concentration steps near stimulus,
    % larger steps later in the beat. This matches the validated tensor code.
    if use_adaptive && mod(ti, cfg.BCL) >= cfg.twin
        dt = dt_coarse;
        Nsplit = Nsplit_coarse;
    else
        dt = dt_fine;
        Nsplit = Nsplit_fine;
    end
    if dt > cfg.T - ti
        dt = cfg.T - ti;
    end
    p.dt = dt;

    % Note mapping:
    %   V_m    = phi_m - phi_bulk, with bulk reference 0.
    %   V_pre  = phi_pre  - phi_c.
    %   V_post = phi_post - phi_c.
    % The ionic model still expects one patch-ordered Vm vector, so this is
    % the point where the readable fields are temporarily stacked.
    V_m = phi_m;
    V_pre = phi_pre - phi_c;
    V_post = phi_post - phi_c;
    V_membrane = stack_membrane_voltage( ...
        V_m, V_pre, V_post, phi_boundary, network);

    % ID patch ionic currents depend on the local cleft concentrations.
    % Axial/boundary patches see the fixed bulk concentration.
    S_patch = build_patch_concentration(s_c, s_b, network);

    % Ionic model evaluation:
    %   Iion is the total membrane current per patch, used in the electrical
    %   potential equation through -Am'*Iion.
    %   Ivec splits the ionic current by species for the cleft s_c equation.
    [Gnew, Iion, Ivec, ~, ~] = ...
        model.ionic_fun(ti, [V_membrane; Gstate], p, S_patch(:));
    Ivec = reshape(Ivec, Npatches, 3);

    % Convert patch-order species currents into note-style quantities:
    %   i_pre_s(junction patch)
    %   i_post_s(junction patch)
    % These are the i_pre^(s) and i_post^(s) terms in ds_c/dt.
    [i_pre_s, i_post_s] = split_ID_ionic_current(Ivec, topology, network);
    Esource = zeros(M, Njunction);
    if ~all(cfg.clamp_flag)
        % Cleft concentration equation from the note:
        %   ds_c^s/dt = 1/(z_s F delta_c) *
        %     [i_pre^s + i_post^s - I_b^s - I_c^s].
        % Esource is the electro-diffusive source term needed by the cleft
        % potential equation after moving Nernst terms to the RHS.
        [s_c, Esource] = update_cleft_concentration_continuum( ...
            phi_c, s_c, s_b, z_s, alpha_s, i_pre_s, i_post_s, ...
            mesh, cfg, dt, Nsplit);
    end

    % Electrical RHS in mass-matrix form:
    %   Cm*dphi/dt = -G*phi - Am'*Iion + cleft_source.
    %
    % The variables are stacked only here because G, Cm, and Am are sparse
    % global operators. The model-facing state remains phi_m/phi_pre/etc.
    phi_old = stack_phi_continuum( ...
        phi_m, phi_pre, phi_post, phi_boundary, phi_c, network);
    forcing_phi = potential_forcing_continuum(Iion, Esource, Nnodes, network);
    rhs_phi = potential_rhs_continuum(phi_old, forcing_phi, network);

    % Numerical integration is isolated here. Changing from backward Euler
    % to Crank-Nicolson, or later to a DAE solver, should happen in this
    % helper rather than by changing the physical RHS construction above.
    phi_new = advance_phi_continuum( ...
        phi_old, rhs_phi, forcing_phi, dt, cfg, network);

    if any(~isfinite(phi_new)) || any(~isfinite(Gnew))
        error("Non-finite state at t = %.6f ms.", ti);
    end

    % Return from sparse solver ordering back to note-style fields.
    [phi_m, phi_pre, phi_post, phi_boundary, phi_c] = ...
        unstack_phi_continuum(phi_new, topology, network);
    Gstate = Gnew;
    ti = round(ti + dt, 5);
    step = step + 1;

    if mod(step, save_every) == 0 || ti >= cfg.T
        count = count + 1;
        if count > size(time, 2)
            new_size = 2 * size(time, 2);
            time(1, new_size) = nan;
            phi_axial(Ncell, new_size) = nan;
            Gstate_hist(NG, new_size) = nan;
            Icleft(2, max(Njunction, 1), new_size) = nan;
            S_cleft_hist(4, M, max(Njunction, 1), new_size) = nan;
        end
        time(count) = ti;
        phi_axial(:, count) = phi_m;
        Gstate_hist(:, count) = Gstate;
        Icleft(:, :, count) = compute_icleft_continuum( ...
            phi_m, phi_pre, phi_post, topology, network.gmyo, M);
        if Njunction > 0
            S_cleft_hist(:, :, :, count) = s_c;
        end
    end

    if cfg.show_progress && ...
            abs(ti / cfg.progress_interval - ...
            round(ti / cfg.progress_interval)) < 1e-8
        fprintf("t = %.3f / %.3f ms\n", ti, cfg.T);
    end
end

result.time = time(1:count);
result.phi_axial = phi_axial(:, 1:count);
result.Gstate = Gstate_hist(:, 1:count);
result.Icleft = Icleft(:, :, 1:count);
result.S_cleft = S_cleft_hist(:, :, :, 1:count);
end

function phi = stack_phi_continuum( ...
    phi_m, phi_pre, phi_post, phi_boundary, phi_c, network)
%STACK_PHI_CONTINUUM Pack readable note fields into global node order.
% This is bookkeeping for sparse operators only. It is the inverse of
% unstack_phi_continuum and should not be treated as the conceptual state.
layout = network.layout;
M = network.M;
Njunction = size(phi_c, 2);

phi = zeros(network.Nnodes, 1);
phi(layout.cell(:)) = phi_m(:);
if Njunction > 0
    phi_ID = zeros(M, 2, Njunction);
    phi_ID(:, 1, :) = reshape(phi_pre, M, 1, Njunction);
    phi_ID(:, 2, :) = reshape(phi_post, M, 1, Njunction);
    phi(layout.ID(:)) = phi_ID(:);
    phi(layout.cleft(:)) = phi_c(:);
end
boundary_mask = network.boundary_mask;
phi(layout.boundary(boundary_mask)) = phi_boundary(boundary_mask);
end

function [phi_m, phi_pre, phi_post, phi_boundary, phi_c] = ...
    unstack_phi_continuum(phi, topology, network)
%UNSTACK_PHI_CONTINUUM Restore note-style field names after a node solve.
layout = network.layout;
M = network.M;
Njunction = topology.Njunction;

phi_m = reshape(phi(layout.cell), topology.Ncell, 1);
phi_pre = zeros(M, Njunction);
phi_post = zeros(M, Njunction);
phi_c = zeros(M, Njunction);
if Njunction > 0
    phi_ID = reshape(phi(layout.ID), M, 2, Njunction);
    phi_pre = reshape(phi_ID(:, 1, :), M, Njunction);
    phi_post = reshape(phi_ID(:, 2, :), M, Njunction);
    phi_c = reshape(phi(layout.cleft), M, Njunction);
end
phi_boundary = zeros(network.Nports, topology.Ncell);
boundary_mask = network.boundary_mask;
phi_boundary(boundary_mask) = phi(layout.boundary(boundary_mask));
end

function V_membrane = stack_membrane_voltage( ...
    V_m, V_pre, V_post, V_boundary, network)
%STACK_MEMBRANE_VOLTAGE Build the patch-ordered Vm expected by ionic_fun.
% Conceptually this is Am*phi, but using already separated note variables is
% more readable here than forming Am*phi just to unpack it again.
M = network.M;
Njunction = size(V_pre, 2);
V_ID = zeros(M, 2, Njunction);
if Njunction > 0
    V_ID(:, 1, :) = reshape(V_pre, M, 1, Njunction);
    V_ID(:, 2, :) = reshape(V_post, M, 1, Njunction);
end
V_membrane = [ ...
    V_m(:); ...
    V_ID(:); ...
    V_boundary(network.boundary_mask)];
end

function S_patch = build_patch_concentration(s_c, s_b, network)
%BUILD_PATCH_CONCENTRATION Give each membrane patch its extracellular ions.
% Axial and boundary membrane patches use the bulk reservoir s_b. ID patches
% look up the local cleft concentration s_c at the facing cleft patch.
patch = network.patch;
M = network.M;
Njunction = size(s_c, 3);
S_patch = repmat(s_b(1:3)', network.Npatches, 1);
if Njunction == 0
    return;
end
ID_patches = patch.ID(:);
cleft_index = patch.cleft_linear(ID_patches);
for ion = 1:3
    values = reshape(s_c(ion, :, :), M * Njunction, 1);
    S_patch(ID_patches, ion) = values(cleft_index);
end
end

function [i_pre_s, i_post_s] = split_ID_ionic_current(Ivec, topology, network)
%SPLIT_ID_IONIC_CURRENT Convert ionic_fun output into note notation.
% Ivec is patch-by-ion over all membrane patches. The cleft equation needs
% species currents from the two facing ID membranes separately.
M = network.M;
Njunction = topology.Njunction;
i_pre_s = zeros(4, M, Njunction);
i_post_s = zeros(4, M, Njunction);
for junction = 1:Njunction
    pre_patch = network.patch.ID(:, 1, junction);
    post_patch = network.patch.ID(:, 2, junction);
    for ion = 1:3
        i_pre_s(ion, :, junction) = reshape(Ivec(pre_patch, ion), 1, M);
        i_post_s(ion, :, junction) = reshape(Ivec(post_patch, ion), 1, M);
    end
end
end

function [s_c, Esource] = update_cleft_concentration_continuum( ...
    phi_c, s_c, s_b, z_s, alpha_s, i_pre_s, i_post_s, ...
    mesh, cfg, dt, Nsplit)
%UPDATE_CLEFT_CONCENTRATION_CONTINUUM Evaluate the note's cleft equations.
% sigma_c_s and bar_g_b are the FEM-derived conductance matrices/vectors.
%
% For species s:
%   eta_c^s = phi_c + alpha_s log(s_c^s)
%   eta_b^s = alpha_s log(s_b^s)
%
% The lateral term I_c_s and bulk term I_b_s are exactly the discrete
% patch-graph versions of
%   div(sigma_c^s grad eta_c^s)
% and
%   bar_g_b(eta_c^s - eta_b^s).
%
% Esource collects the Nernst/electro-diffusive source terms that enter the
% cleft potential equation when the electrical solve is written in phi.

M = size(phi_c, 1);
Njunction = size(phi_c, 2);
dtS = dt / Nsplit;
Esource = zeros(M, Njunction);

for split_step = 1:Nsplit 
    Esource(:) = 0;
    for junction = 1:Njunction
        phi_c_j = phi_c(:, junction);
        voltage_difference = phi_c_j - phi_c_j';
        sigma_c_s = mesh.Gc(:, :, junction) / 4;
        bar_g_b = mesh.Gb(:, junction) / 4;

        for ion = 1:4
            concentration = reshape(s_c(ion, :, junction), M, 1);
            eta_c = phi_c_j + alpha_s(ion) * log(concentration);
            eta_b = alpha_s(ion) * log(s_b(ion)); 

            Erev = alpha_s(ion) * log(concentration' ./ concentration);
            Ebulk = alpha_s(ion) * log(s_b(ion) ./ concentration);
            I_c_s = sum(sigma_c_s .* (voltage_difference - Erev), 2);
            I_b_s = bar_g_b .* (phi_c_j - Ebulk);
            Iout = I_c_s + I_b_s;

            Esource(:, junction) = Esource(:, junction) + ...
                sum(sigma_c_s .* Erev, 2) + bar_g_b .* Ebulk;

            if ~cfg.clamp_flag(ion)
                i_pre_post_s = reshape( ...
                    i_pre_s(ion, :, junction) + i_post_s(ion, :, junction), ...
                    M, 1);
                ds_c_dt = (i_pre_post_s - Iout) * ...
                    1e6 ./ (z_s(ion) * cfg.F * mesh.volume(:, junction));
                s_c(ion, :, junction) = ...
                    reshape(max(concentration + dtS * ds_c_dt, 1e-9), 1, M);
            end

            if any(~isfinite(eta_c))
                error("Non-finite eta_c for ion %d at junction %d.", ion, junction);
            end
        end
    end
end
end

function Icleft = compute_icleft_continuum( ...
    phi_m, phi_pre, phi_post, topology, gmyo, M)
%COMPUTE_ICLEFT Report axial-to-ID current on each side of each junction.
% This is the discrete version of
%   g_myo(phi_m - <phi_pre>) and g_myo(phi_m - <phi_post>).
Njunction = topology.Njunction;
Icleft = zeros(2, Njunction);
for junction = 1:Njunction
    pre_cell = topology.junction_cells(1, junction);
    post_cell = topology.junction_cells(2, junction);
    Icleft(1, junction) = (gmyo / M) * ...
        (M * phi_m(pre_cell) - sum(phi_pre(:, junction)));
    Icleft(2, junction) = (gmyo / M) * ...
        (M * phi_m(post_cell) - sum(phi_post(:, junction)));
end
end

function forcing_phi = potential_forcing_continuum( ...
    Iion, Esource, Nnodes, network)
%POTENTIAL_FORCING_CONTINUUM Terms in the electrical equation not in -G*phi.
% -Am'*Iion projects membrane ionic currents into nodal KCL equations.
% Esource is nonzero only on cleft nodes and comes from cleft ion gradients.
forcing_phi = -network.Am' * Iion;
cleft_source = zeros(Nnodes, 1);
cleft_source(network.layout.cleft(:)) = -Esource(:);
forcing_phi = forcing_phi + cleft_source;
end

function rhs_phi = potential_rhs_continuum(phi, forcing_phi, network)
%POTENTIAL_RHS_CONTINUUM Semi-discrete electrical equation:
%   C_m * dphi/dt = -G * phi - A_m' * Iion + source.
rhs_phi = -network.G * phi + forcing_phi;
end

function phi_new = advance_phi_continuum( ...
    phi_old, rhs_phi, forcing_phi, dt, cfg, network)
%ADVANCE_PHI_CONTINUUM Choose the numerical update for the electrical part.
% The continuum RHS is assembled before this function. Only this function
% knows whether the update is explicit, implicit, or midpoint-like.

if isfield(cfg, "electric_integrator")
    method = lower(string(cfg.electric_integrator));
else
    method = "backward_euler";
end

switch method
    case "backward_euler"
        lhs = network.Cm / dt + network.G;
        rhs = network.Cm / dt * phi_old + forcing_phi;
        phi_new = lhs \ rhs;

    case "forward_euler"
        error([
            "forward_euler requires solving dphi/dt = C_m \\ RHS, " + ...
            "but this ID-cleft potential system has a singular mass " + ...
            "matrix and is a mass-matrix DAE. Use backward_euler, " + ...
            "crank_nicolson, or a DAE/mass-matrix ODE solver."
        ]);

    case "crank_nicolson"
        lhs = network.Cm / dt + 0.5 * network.G;
        rhs = (network.Cm / dt - 0.5 * network.G) * phi_old + forcing_phi;
        phi_new = lhs \ rhs;

    otherwise
        error("Unknown electric_integrator: %s", method);
end
end
