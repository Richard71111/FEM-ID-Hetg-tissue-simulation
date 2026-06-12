function topology = build_topology(cell_mask)
%BUILD_TOPOLOGY Convert an occupancy array into cells and junctions.
% Build topological node and edge structure.
% Input: logical/numeric 1-D, 2-D, or 3-D occupancy array.
% Example 5 cell network:
%     0, 1, 0;
%     1, 1, 1;
%     0, 1, 0
% 1 represent cell/ 0 represent empty
% 
% Important parameters
% 
% junction_cells: Mark down the connection order
% Example:
%          2     3     1     3
%          3     4     3     5
%     Cell 2 connect to cell 3, cell 3 connect to cell 4 .etc. read by col
% junction_axis: Mark down connection orientation
% Example:
%          [1 1 2 2]
%     Cell (2,3) (3,4) is first dim connection(connect through col)
% 
% Output: topology with cell coordinates and a 2-by-Njunction edge list.

if ndims(cell_mask) > 3
    error("The topology array must have at most three dimensions.");
end

cell_mask = logical(cell_mask);
array_size = size(cell_mask); % dimension array size
array_size(end + 1:3) = 1; % make mask to 3D matrix [1,2] -> [1,2,1]
cell_mask = reshape(cell_mask, array_size);

if ~any(cell_mask, "all")
    error("The topology must contain at least one occupied cell.");
end

active_dims = find(array_size > 1);
spatial_dimension = max(1, numel(active_dims));

cell_id = zeros(array_size);
cell_id(cell_mask) = 1:nnz(cell_mask);
Ncell = nnz(cell_mask);

[sub_1, sub_2, sub_3] = ind2sub(array_size, find(cell_mask));
all_coordinates = [sub_1(:), sub_2(:), sub_3(:)];
if isempty(active_dims)
    coordinates = zeros(Ncell, 1);
else
    coordinates = all_coordinates(:, active_dims);
end

junction_cells = zeros(2, 0);
junction_axis = zeros(1, 0);
% Face order is [-axis 1, +axis 1, -axis 2, +axis 2, ...].
junction_faces = zeros(2, 0);

for array_dim = active_dims
    first = {':', ':', ':'};
    second = first;
    first{array_dim} = 1:(array_size(array_dim) - 1);
    second{array_dim} = 2:array_size(array_dim);

    cell_1 = cell_id(first{:});
    cell_2 = cell_id(second{:});
    connected = cell_1 > 0 & cell_2 > 0;

    junction_cells = [junction_cells, ...
        [reshape(cell_1(connected), 1, []); ...
        reshape(cell_2(connected), 1, [])]];
    local_axis = find(active_dims == array_dim);
    junction_axis = [junction_axis, ...
        local_axis * ones(1, nnz(connected))];
    junction_faces = [junction_faces, ...
        [(2 * local_axis) * ones(1, nnz(connected)); ...
        (2 * local_axis - 1) * ones(1, nnz(connected))]];
end

Njunction = size(junction_cells, 2);
degree = accumarray(junction_cells(:), 1, [Ncell, 1]);

topology.dimension = spatial_dimension;
topology.mask = cell_mask;
topology.array_size = array_size;
topology.cell_id = cell_id;
topology.coordinates = coordinates;
topology.Ncell = Ncell;
topology.Njunction = Njunction;
topology.junction_cells = junction_cells;
topology.junction_axis = junction_axis;
topology.junction_faces = junction_faces;
topology.degree = degree;
end
