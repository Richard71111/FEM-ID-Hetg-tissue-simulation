function plot_graph_result(result)
%PLOT_GRAPH_RESULT Plot axial potential and a dimension-aware overview.
% Input: result returned by run_graph_simulation (must contain phi_axial).
% Output:
%   - Figure 1: axial potential phi_axial versus time (all cells).
%   - 1-D networks: a 3-D surface over (cell number, time, phi).
%   - 2-D / 3-D networks: a time-lapse video of the node potentials.
% Figures and the video are saved to cfg.save_plot_path.

topology = result.topology;
cfg = result.cfg;
coordinates = topology.coordinates;
edges = topology.junction_cells;
dimension = topology.dimension;
Ncell = topology.Ncell;
Njunction = topology.Njunction;

time = result.time(:);                 % Nt-by-1.
phi_axial = result.phi_axial;          % Ncell-by-Nt.

if isfield(cfg, "save_plot_path") && ~isempty(cfg.save_plot_path)
    save_path = cfg.save_plot_path;
    if ~exist(save_path, "dir")
        mkdir(save_path);
    end
else
    save_path = pwd;
end

if isfield(cfg, "mesh_files") && ~isempty(cfg.mesh_files)
    mesh_tag = char(cfg.mesh_files(1));
else
    mesh_tag = "mesh";
end
tag = sprintf("%s_BCL_%d_nbeats_%d_Ncell_%d_%dD", ...
    mesh_tag, cfg.BCL, cfg.nbeats, Ncell, dimension);

%% Figure 1: axial potential versus time
fig1 = figure("Color", "w", "Name", "Axial potential");
plot(time, phi_axial', "LineWidth", 1.2);
xlabel("Time (ms)");
ylabel("$\phi_{axial}$ (mV)", "Interpreter", "latex");
title("Axial potential");
grid on;
box on;
saveas(fig1, fullfile(save_path, "axial_phi_" + tag + "_line.png"));

%% Dimension-aware overview
if dimension == 1
    % 3-D surface: x = cell number, y = time, z = phi.
    [Cell, Time] = meshgrid(1:Ncell, time);
    fig2 = figure("Color", "w", "Name", "Axial potential surface");
    surf(Cell, Time, phi_axial.');
    xlabel("Cell number");
    ylabel("Time (ms)");
    zlabel("$\phi_{axial}$ (mV)", "Interpreter", "latex");
    title("Axial potential over space and time");
    shading interp;
    colorbar;
    view(45, 30);
    saveas(fig2, fullfile(save_path, "axial_phi_" + tag + "_surf.png"));
else
    % 2-D / 3-D: animate the node potentials with time as the video axis.
    write_potential_video( ...
        fullfile(save_path, "phi_video_" + tag), ...
        coordinates, edges, Njunction, dimension, time, phi_axial);
end
end

function write_potential_video(base_name, coordinates, edges, Njunction, ...
    dimension, time, phi_axial)
%WRITE_POTENTIAL_VIDEO Render node potentials frame by frame to a video file.

x = coordinates(:, 1);
y = coordinates(:, 2);
if dimension == 3
    z = coordinates(:, 3);
end

% Try MPEG-4, fall back to Motion JPEG AVI where it is unavailable.
try
    writer = VideoWriter(char(base_name + ".mp4"), "MPEG-4");
catch
    writer = VideoWriter(char(base_name + ".avi"), "Motion JPEG AVI");
end
writer.FrameRate = 20;
open(writer);

% Fixed colour scale across all frames.
clim_lo = min(phi_axial(:));
clim_hi = max(phi_axial(:));
if clim_lo == clim_hi
    clim_hi = clim_lo + 1;
end

% Cap the number of rendered frames so long runs stay responsive.
Nt = numel(time);
frame_step = max(1, round(Nt / 300));

fig = figure("Color", "w", "Name", "Node potential video");
for t = 1:frame_step:Nt
    clf(fig);
    ax = axes("Parent", fig);
    hold(ax, "on");
    values = phi_axial(:, t);
    if dimension == 2
        for j = 1:Njunction
            c = edges(:, j);
            plot(ax, x(c), y(c), "-", "Color", [0.6 0.6 0.6], "LineWidth", 1.0);
        end
        scatter(ax, x, y, 160, values, "filled");
        axis(ax, "equal");
        xlabel(ax, "x");
        ylabel(ax, "y");
    else
        for j = 1:Njunction
            c = edges(:, j);
            plot3(ax, x(c), y(c), z(c), "-", ...
                "Color", [0.6 0.6 0.6], "LineWidth", 1.0);
        end
        scatter3(ax, x, y, z, 160, values, "filled");
        axis(ax, "equal");
        view(ax, 35, 25);
        xlabel(ax, "x");
        ylabel(ax, "y");
        zlabel(ax, "z");
    end
    caxis(ax, [clim_lo, clim_hi]);
    cb = colorbar(ax);
    cb.Label.String = "\phi_{axial} (mV)";
    title(ax, sprintf("t = %.1f ms", time(t)));
    hold(ax, "off");
    drawnow;
    writeVideo(writer, getframe(fig));
end
close(writer);
close(fig);
end
