function result = run_time_loop(cfg, topology, model, network)
%RUN_TIME_LOOP Advance ionic states, cleft concentrations, and potentials.
% Inputs: cfg, topology, ionic model, and assembled graph network.
% Output: sampled histories and final values in physical tensor form.

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

Nsplit = max(1, round(cfg.dt / cfg.dtS));
RTF = cfg.R * cfg.Temp / cfg.F;
system_matrix = network.G + network.Cm / cfg.dt;
solver = decomposition(system_matrix, "lu");

%% Sampled output
max_samples = ceil(cfg.T / cfg.sample_dt) + 2;
time = nan(1, max_samples);
Vm_cell = nan(Ncell, max_samples);
Vm_ID_mean = nan(2, Njunction, max_samples);
phi_cleft_mean = nan(Njunction, max_samples);
S_cleft_mean = nan(4, Njunction, max_samples);

Vm = network.Am * phi;
count = 1;
time(count) = 0;
Vm_cell(:, count) = Vm(patch.axial);
if Njunction > 0
    Vm_ID = reshape(Vm(patch.ID(:)), M, 2, Njunction);
    Vm_ID_mean(:, :, count) = reshape(mean(Vm_ID, 1), 2, Njunction);
    phi_cleft_mean(:, count) = mean(phi(layout.cleft), 1)';
    S_cleft_mean(:, :, count) = reshape(mean(S_cleft, 2), 4, Njunction);
end
next_sample = cfg.sample_dt;
ti = 0;
Iall = [];

%% Main simulation loop
while ti < cfg.T
    dt = min(cfg.dt, cfg.T - ti);
    if abs(dt - cfg.dt) > eps(cfg.dt)
        system_matrix = network.G + network.Cm / dt;
        solver = decomposition(system_matrix, "lu");
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

    [Gnew, Iion, Ivec, ~, Iall] = ...
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
    source(layout.cleft(:)) = Esource(:);
    rhs = network.Cm / dt * phi - network.Am' * Iion + source;
    phi_new = solver \ rhs;

    if any(~isfinite(phi_new)) || any(~isfinite(Gnew))
        error("Non-finite state at t = %.6f ms.", ti);
    end

    phi = phi_new;
    Gstate = Gnew;
    ti = round(ti + dt, 10);

    if ti + 1e-10 >= next_sample || ti >= cfg.T
        count = count + 1;
        Vm = network.Am * phi;
        time(count) = ti;
        Vm_cell(:, count) = Vm(patch.axial);
        if Njunction > 0
            Vm_ID = reshape(Vm(patch.ID(:)), M, 2, Njunction);
            Vm_ID_mean(:, :, count) = ...
                reshape(mean(Vm_ID, 1), 2, Njunction);
            phi_cleft_mean(:, count) = mean(phi(layout.cleft), 1)';
            S_cleft_mean(:, :, count) = ...
                reshape(mean(S_cleft, 2), 4, Njunction);
        end
        next_sample = next_sample + cfg.sample_dt;
    end

    if cfg.show_progress && ...
            abs(ti / cfg.progress_interval - ...
            round(ti / cfg.progress_interval)) < 1e-8
        fprintf("t = %.3f / %.3f ms\n", ti, cfg.T);
    end
end

result.time = time(1:count);
result.Vm_cell = Vm_cell(:, 1:count);
result.Vm_ID_mean = Vm_ID_mean(:, :, 1:count);
result.phi_cleft_mean = phi_cleft_mean(:, 1:count);
result.S_cleft_mean = S_cleft_mean(:, :, 1:count);
result.final.phi_cell = reshape(phi(layout.cell), 1, Ncell);
result.final.phi_ID = reshape(phi(layout.ID), M, 2, Njunction);
result.final.phi_boundary = nan(network.Nfaces, Ncell);
result.final.phi_boundary(network.boundary_mask) = ...
    phi(layout.boundary(network.boundary_mask));
result.final.phi_cleft = reshape(phi(layout.cleft), M, Njunction);
result.final.S_cleft = S_cleft;
result.final.Gstate = Gstate;
result.final.Iall = Iall;
end
