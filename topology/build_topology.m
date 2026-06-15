function topology = build_topology(adjacency_matrix, cell_coordinates, cell_port_count)
%BUILD_TOPOLOGY Convert a cell adjacency matrix into graph metadata.
% Inputs: symmetric 0/1 adjacency matrix, optional plotting coordinates,
% and optional total membrane-port count for each cell.
% Output: topology with cells, junctions, graph ports, and coordinates.

if ~ismatrix(adjacency_matrix) || ...
        size(adjacency_matrix, 1) ~= size(adjacency_matrix, 2)
    error("adjacency_matrix must be a square matrix.");
end
if any(~isfinite(adjacency_matrix), "all") || ...
        any(adjacency_matrix(:) ~= 0 & adjacency_matrix(:) ~= 1)
    error("adjacency_matrix must contain only finite 0 and 1 values.");
end
if any(diag(adjacency_matrix) ~= 0)
    error("adjacency_matrix must have a zero diagonal.");
end
if ~isequal(adjacency_matrix, adjacency_matrix')
    error("adjacency_matrix must be symmetric for an undirected tissue graph.");
end

adjacency_matrix = logical(adjacency_matrix);
Ncell = size(adjacency_matrix, 1);
if Ncell == 0
    error("adjacency_matrix must contain at least one cell.");
end
component = conncomp(graph(adjacency_matrix))';
is_connected = all(component == component(1));

[cell_1, cell_2] = find(triu(adjacency_matrix, 1));
junction_cells = [cell_1'; cell_2'];
Njunction = numel(cell_1);
degree = full(sum(adjacency_matrix, 2));

if nargin < 3 || isempty(cell_port_count)
    cell_port_count = max(2, degree);
elseif isscalar(cell_port_count)
    cell_port_count = repmat(cell_port_count, Ncell, 1);
else
    cell_port_count = cell_port_count(:);
end
if numel(cell_port_count) ~= Ncell || ...
        any(cell_port_count < degree) || ...
        any(cell_port_count < 1) || ...
        any(cell_port_count ~= round(cell_port_count))
    error("cell_port_count must contain positive integers not smaller than cell degree.");
end

junction_ports = zeros(2, Njunction);
next_port = zeros(Ncell, 1);
for junction = 1:Njunction
    for side = 1:2
        cell_number = junction_cells(side, junction);
        next_port(cell_number) = next_port(cell_number) + 1;
        junction_ports(side, junction) = next_port(cell_number);
    end
end

Nports = max(cell_port_count);
boundary_mask = false(Nports, Ncell);
for cell_number = 1:Ncell
    boundary_mask(1:cell_port_count(cell_number), cell_number) = true;
end
for junction = 1:Njunction
    for side = 1:2
        boundary_mask( ...
            junction_ports(side, junction), ...
            junction_cells(side, junction)) = false;
    end
end

if nargin < 2 || isempty(cell_coordinates)
    if Ncell == 1
        cell_coordinates = 0;
    elseif is_connected && Njunction == Ncell - 1 && all(degree <= 2)
        path_order = zeros(Ncell, 1);
        current_cell = find(degree == 1, 1);
        previous_cell = 0;
        for position = 1:Ncell
            path_order(position) = current_cell;
            next_cells = find(adjacency_matrix(current_cell, :));
            next_cells(next_cells == previous_cell) = [];
            if position < Ncell
                previous_cell = current_cell;
                current_cell = next_cells(1);
            end
        end
        cell_coordinates = zeros(Ncell, 1);
        cell_coordinates(path_order) = (0:Ncell - 1)';
    else
        angle = 2 * pi * (0:Ncell - 1)' / Ncell;
        cell_coordinates = [cos(angle), sin(angle)];
    end
else
    if size(cell_coordinates, 1) ~= Ncell || ...
            size(cell_coordinates, 2) < 1 || ...
            size(cell_coordinates, 2) > 3 || ...
            any(~isfinite(cell_coordinates), "all")
        error("cell_coordinates must be a finite Ncell-by-1, 2, or 3 matrix.");
    end
end

topology.dimension = size(cell_coordinates, 2);
topology.adjacency = adjacency_matrix;
topology.coordinates = cell_coordinates;
topology.Ncell = Ncell;
topology.Njunction = Njunction;
topology.junction_cells = junction_cells;
topology.junction_ports = junction_ports;
topology.junction_faces = junction_ports;
topology.degree = degree;
topology.component = component;
topology.cell_port_count = cell_port_count;
topology.Nports = Nports;
topology.boundary_mask = boundary_mask;
end
