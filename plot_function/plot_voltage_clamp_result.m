function plot_voltage_clamp_result(result, figure_name)
%PLOT_VOLTAGE_CLAMP_RESULT Plot voltage-clamp voltages and cleft currents.
% The figure includes axial voltage, command voltage, and left/right cleft
% current with explicit legends for cell and cleft side.

cfg = result.cfg;
if ~exist(cfg.save_plot_path, "dir")
    mkdir(cfg.save_plot_path);
end

time = result.time(:);
phi_axial = result.phi_axial;
command = result.vclamp.command;
Icleft = result.Icleft;

Ncell = size(phi_axial, 1);
cell_labels = compose("Cell %d", 1:Ncell);
clamp_labels = compose("Cell %d command", 1:Ncell);

fig = figure("Color", "w", "Name", figure_name);
tiledlayout(fig, 3, 1, "TileSpacing", "compact");

nexttile;
plot(time, phi_axial', "LineWidth", 1.2);
ylabel("\phi_{axial} (mV)");
title(figure_name);
legend(cell_labels, "Location", "best");
grid on;

nexttile;
command_rows = any(isfinite(command), 2);
plot(time, command(command_rows, :)', "--", "LineWidth", 1.2);
ylabel("V command (mV)");
legend(clamp_labels(command_rows), "Location", "best");
grid on;

nexttile;
plot_cleft_current_panel(time, Icleft);
xlabel("Time (ms)");
ylabel("I cleft (uA)");
grid on;

file_tag = regexprep(char(figure_name), '[<>:"/\\|?*]', '_');
save_file = fullfile(cfg.save_plot_path, file_tag + ".png");
saveas(fig, save_file);
fprintf("Saved %s\n", save_file);
end

function plot_cleft_current_panel(time, Icleft)
Njunction = size(Icleft, 2);
if Njunction == 0
    plot(time, zeros(size(time)), "LineWidth", 1.2);
    legend("No cleft junction", "Location", "best");
    return;
end

hold on;
labels = strings(1, 2 * Njunction);
plot_index = 0;
for junction = 1:Njunction
    left_current = squeeze(Icleft(1, junction, :));
    right_current = squeeze(Icleft(2, junction, :));

    plot_index = plot_index + 1;
    plot(time, left_current, "LineWidth", 1.2);
    labels(plot_index) = sprintf("Junction %d left cleft current", junction);

    plot_index = plot_index + 1;
    plot(time, right_current, "--", "LineWidth", 1.2);
    labels(plot_index) = sprintf("Junction %d right cleft current", junction);
end
hold off;
legend(labels, "Location", "best");
end
