%% COMPARE_OLD_NEW
% Numerical comparison of the OLD 1-D code and the NEW tensor/graph code.
% This script does ONLY analysis: it loads the two saved result files and
% reports errors. Generate the inputs first by running:
%     (old code)  run_old_for_compare.m   -> Output/Save data/old_model.mat
%     (new code)  run_new_for_compare.m   -> Output/Save data/new_model.mat
%
% Compared quantities (all as time histories on a common time grid):
%     phi_axial : axial potential per cell        (Ncell x Nt)
%     Gsum      : sum of each ionic state over all membrane patches (Nstate x Nt)
%     Ssum      : sum of each cleft ion over all cleft nodes         (4 x Nt)
% G and S are compared as sums so the result does not depend on the
% (different) internal patch/node ordering of the two codes.

clear; clc; close all;
project_folder = fileparts(mfilename("fullpath"));
data_folder = fullfile(project_folder, "Output", "Save data");

old_file = fullfile(data_folder, "old_model.mat");
new_file = fullfile(data_folder, "new_model.mat");
assert(isfile(old_file), "Missing %s. Run run_old_for_compare.m first.", old_file);
assert(isfile(new_file), "Missing %s. Run run_new_for_compare.m first.", new_file);

S = load(old_file); old = S.old;
S = load(new_file); new = S.new;

%% Basic consistency checks
fprintf("\n========== SETUP ==========\n");
fprintf("Ncell    old=%d  new=%d\n", old.Ncell, new.Ncell);
fprintf("Nstate   old=%d  new=%d\n", old.Nstate, new.Nstate);
fprintf("Npatches old=%d  new=%d\n", old.Npatches, new.Npatches);
fprintf("samples  old=%d  new=%d\n", numel(old.time), numel(new.time));
if old.Ncell ~= new.Ncell || old.Nstate ~= new.Nstate
    warning("Ncell/Nstate differ - comparison may be invalid.");
end

%% Common time grid (interpolate onto the overlap; grids should already match)
t0 = max(old.time(1),  new.time(1));
t1 = min(old.time(end), new.time(end));
tc = old.time(old.time >= t0 & old.time <= t1);
tc = tc(:)';

interp_rows = @(t, Y, tq) interp1(t(:), Y.', tq(:), "linear").';  % rows x numel(tq)

phi_o = interp_rows(old.time, old.phi_axial, tc);
phi_n = interp_rows(new.time, new.phi_axial, tc);
G_o   = interp_rows(old.time, old.Gsum,      tc);
G_n   = interp_rows(new.time, new.Gsum,      tc);
S_o   = interp_rows(old.time, old.Ssum,      tc);
S_n   = interp_rows(new.time, new.Ssum,      tc);

%% Error metrics
maxabs = @(a, b) max(abs(b(:) - a(:)));
rmse   = @(a, b) sqrt(mean((b(:) - a(:)).^2));
relL2  = @(a, b) norm(b(:) - a(:)) / max(norm(a(:)), eps);

report = @(name, a, b) fprintf("%-12s  maxabs=%.3e  rmse=%.3e  relL2=%.3e\n", ...
    name, maxabs(a, b), rmse(a, b), relL2(a, b));

fprintf("\n========== OVERALL ERRORS (new vs old) ==========\n");
report("phi_axial", phi_o, phi_n);
report("Gsum",      G_o,   G_n);
report("Ssum",      S_o,   S_n);

fprintf("\n----- phi_axial per cell -----\n");
for c = 1:size(phi_o, 1)
    report(sprintf("cell %d", c), phi_o(c, :), phi_n(c, :));
end

ion_names = ["Na", "K", "Ca", "A"];
fprintf("\n----- Ssum per ion (cleft) -----\n");
for i = 1:size(S_o, 1)
    report(ion_names(i), S_o(i, :), S_n(i, :));
end

% states with the largest relative error
relG = zeros(size(G_o, 1), 1);
for s = 1:size(G_o, 1)
    relG(s) = relL2(G_o(s, :), G_n(s, :));
end
[~, order] = sort(relG, "descend");
fprintf("\n----- Gsum: 5 states with largest relL2 -----\n");
for k = 1:min(5, numel(order))
    s = order(k);
    fprintf("state %2d   relL2=%.3e\n", s, relG(s));
end

%% Plots
% 1) phi_axial overlay
f1 = figure("Color", "w", "Name", "phi_axial old vs new");
hold on;
co = lines(size(phi_o, 1));
for c = 1:size(phi_o, 1)
    plot(tc, phi_o(c, :), "-",  "Color", co(c, :), "LineWidth", 1.4);
    plot(tc, phi_n(c, :), "--", "Color", co(c, :), "LineWidth", 1.0);
end
xlabel("Time (ms)"); ylabel("\phi_{axial} (mV)");
title("phi_{axial}: old (solid) vs new (dashed)"); grid on; box on;
saveas(f1, fullfile(data_folder, "cmp_phi_axial.png"));

% 2) phi_axial error
f2 = figure("Color", "w", "Name", "phi_axial error");
plot(tc, (phi_n - phi_o)', "LineWidth", 1.0);
xlabel("Time (ms)"); ylabel("\phi_{new} - \phi_{old} (mV)");
title("Axial potential error per cell"); grid on; box on;
saveas(f2, fullfile(data_folder, "cmp_phi_axial_error.png"));

% 3) Ssum overlay
f3 = figure("Color", "w", "Name", "Ssum old vs new");
hold on;
ci = lines(4);
for i = 1:4
    plot(tc, S_o(i, :), "-",  "Color", ci(i, :), "LineWidth", 1.4);
    plot(tc, S_n(i, :), "--", "Color", ci(i, :), "LineWidth", 1.0);
end
xlabel("Time (ms)"); ylabel("\Sigma cleft concentration");
title("Ssum per ion: old (solid) vs new (dashed)");
legend(ion_names, "Location", "best"); grid on; box on;
saveas(f3, fullfile(data_folder, "cmp_Ssum.png"));

%% Save a summary
summary = struct();
summary.tc = tc;
summary.phi_axial = struct("maxabs", maxabs(phi_o, phi_n), ...
    "rmse", rmse(phi_o, phi_n), "relL2", relL2(phi_o, phi_n));
summary.Gsum = struct("maxabs", maxabs(G_o, G_n), ...
    "rmse", rmse(G_o, G_n), "relL2", relL2(G_o, G_n), "relL2_per_state", relG);
summary.Ssum = struct("maxabs", maxabs(S_o, S_n), ...
    "rmse", rmse(S_o, S_n), "relL2", relL2(S_o, S_n));
save(fullfile(data_folder, "compare_summary.mat"), "summary");

fprintf("\nFigures and compare_summary.mat saved to:\n  %s\n", data_folder);
