function plot_graph_result(result)
%PLOT_GRAPH_RESULT Plot cell voltage, topology, and mean cleft voltage.
% Input: result returned by run_graph_simulation.
% Output: three MATLAB figures.

topology = result.topology;
coordinates = topology.coordinates;
edges = topology.junction_cells;

figure("Color", "w", "Name", "Cell membrane voltage");
plot(result.time, result.Vm_cell', "LineWidth", 1.2);
xlabel("Time (ms)");
ylabel("Cell membrane voltage (mV)");
legend(compose("Cell %d", 1:topology.Ncell), "Location", "best");
grid on;
box on;

figure("Color", "w", "Name", "Cell topology");
hold on;
if topology.dimension == 1
    x = coordinates(:, 1);
    y = zeros(topology.Ncell, 1);
    for j = 1:topology.Njunction
        cells = edges(:, j);
        plot(x(cells), y(cells), "k-", "LineWidth", 1.5);
    end
    scatter(x, y, 120, result.Vm_cell(:, end), "filled");
    text(x, y, compose("  %d", 1:topology.Ncell));
    axis padded;
elseif topology.dimension == 2
    x = coordinates(:, 1);
    y = coordinates(:, 2);
    for j = 1:topology.Njunction
        cells = edges(:, j);
        plot(x(cells), y(cells), "k-", "LineWidth", 1.5);
    end
    scatter(x, y, 140, result.Vm_cell(:, end), "filled");
    text(x, y, compose("  %d", 1:topology.Ncell));
    axis equal padded;
else
    x = coordinates(:, 1);
    y = coordinates(:, 2);
    z = coordinates(:, 3);
    for j = 1:topology.Njunction
        cells = edges(:, j);
        plot3(x(cells), y(cells), z(cells), "k-", "LineWidth", 1.5);
    end
    scatter3(x, y, z, 140, result.Vm_cell(:, end), "filled");
    text(x, y, z, compose("  %d", 1:topology.Ncell));
    axis equal;
    view(35, 25);
end
colorbar;
title("Final cell voltage and graph connectivity");
hold off;

figure("Color", "w", "Name", "Mean cleft voltage");
if topology.Njunction > 0
    plot(result.time, result.phi_cleft_mean', "LineWidth", 1.2);
    legend(compose("Junction %d", 1:topology.Njunction), ...
        "Location", "best");
else
    plot(result.time, zeros(size(result.time)), "LineWidth", 1.2);
end
xlabel("Time (ms)");
ylabel("Mean cleft potential (mV)");
grid on;
box on;
end
