%% COMPARE_CONTINUUM_WITH_TENSOR
% Compare the note-style continuum implementation against the working
% tensor/graph output. The comparison uses the same reduced quantities as
% compare_old_new.m: phi_axial, Gsum, and Ssum.

clear;
clc;
close all;

project_folder = fileparts(mfilename("fullpath"));
data_folder = fullfile(project_folder, "Output", "Save data");

continuum_file = fullfile(data_folder, "continuum_model.mat");
assert(isfile(continuum_file), ...
    "Missing %s. Run run_continuum_for_compare.m first.", continuum_file);

reference_candidates = [
    fullfile(data_folder, "new_model.mat")
    "C:\Users\Administrator\Desktop\electric_sim_tensor_graph\Output\Save data\new_model.mat"
];
reference_file = "";
for k = 1:numel(reference_candidates)
    if isfile(reference_candidates(k))
        reference_file = reference_candidates(k);
        break;
    end
end
assert(reference_file ~= "", ...
    "Missing tensor reference new_model.mat. Run the working tensor graph benchmark first.");

S = load(reference_file);
if isfield(S, "new")
    old = S.new;
elseif isfield(S, "reference")
    old = S.reference;
else
    error("Reference file must contain variable new or reference.");
end
S = load(continuum_file);
new = S.continuum;

fprintf("\nReference file:\n  %s\n", reference_file);
fprintf("Continuum file:\n  %s\n", continuum_file);

fprintf("\n========== SETUP ==========\n");
fprintf("Ncell    old=%d  new=%d\n", old.Ncell, new.Ncell);
fprintf("Nstate   old=%d  new=%d\n", old.Nstate, new.Nstate);
fprintf("Npatches old=%d  new=%d\n", old.Npatches, new.Npatches);
fprintf("samples  old=%d  new=%d\n", numel(old.time), numel(new.time));

t0 = max(old.time(1), new.time(1));
t1 = min(old.time(end), new.time(end));
tc = old.time(old.time >= t0 & old.time <= t1);
tc = tc(:)';

interp_rows = @(t, Y, tq) interp1(t(:), Y.', tq(:), "linear").';

phi_o = interp_rows(old.time, old.phi_axial, tc);
phi_n = interp_rows(new.time, new.phi_axial, tc);
G_o = interp_rows(old.time, old.Gsum, tc);
G_n = interp_rows(new.time, new.Gsum, tc);
S_o = interp_rows(old.time, old.Ssum, tc);
S_n = interp_rows(new.time, new.Ssum, tc);

maxabs = @(a, b) max(abs(b(:) - a(:)));
rmse = @(a, b) sqrt(mean((b(:) - a(:)).^2));
relL2 = @(a, b) norm(b(:) - a(:)) / max(norm(a(:)), eps);

report = @(name, a, b) fprintf("%-12s  maxabs=%.3e  rmse=%.3e  relL2=%.3e\n", ...
    name, maxabs(a, b), rmse(a, b), relL2(a, b));

fprintf("\n========== OVERALL ERRORS (continuum vs tensor) ==========\n");
report("phi_axial", phi_o, phi_n);
report("Gsum", G_o, G_n);
report("Ssum", S_o, S_n);

fprintf("\n----- phi_axial per cell -----\n");
for c = 1:size(phi_o, 1)
    report(sprintf("cell %d", c), phi_o(c, :), phi_n(c, :));
end

ion_names = ["Na", "K", "Ca", "A"];
fprintf("\n----- Ssum per ion (cleft) -----\n");
for i = 1:size(S_o, 1)
    report(ion_names(i), S_o(i, :), S_n(i, :));
end

relG = zeros(size(G_o, 1), 1);
for s = 1:size(G_o, 1)
    relG(s) = relL2(G_o(s, :), G_n(s, :));
end
[~, order] = sort(relG, "descend");
fprintf("\n----- Gsum: 5 states with largest relL2 -----\n");
for k = 1:min(5, numel(order))
    state_number = order(k);
    fprintf("state %2d   relL2=%.3e\n", state_number, relG(state_number));
end

summary = struct();
summary.tc = tc;
summary.reference_file = reference_file;
summary.continuum_file = continuum_file;
summary.phi_axial = struct("maxabs", maxabs(phi_o, phi_n), ...
    "rmse", rmse(phi_o, phi_n), "relL2", relL2(phi_o, phi_n));
summary.Gsum = struct("maxabs", maxabs(G_o, G_n), ...
    "rmse", rmse(G_o, G_n), "relL2", relL2(G_o, G_n), ...
    "relL2_per_state", relG);
summary.Ssum = struct("maxabs", maxabs(S_o, S_n), ...
    "rmse", rmse(S_o, S_n), "relL2", relL2(S_o, S_n));
save(fullfile(data_folder, "continuum_compare_summary.mat"), "summary");

fprintf("\ncontinuum_compare_summary.mat saved to:\n  %s\n", data_folder);
