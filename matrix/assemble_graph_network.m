function network = assemble_graph_network(cfg, topology, model)
%ASSEMBLE_GRAPH_NETWORK Build tensor layouts and global sparse operators.
% Parameter Notes:
%     layout: flat global index structure.
%         .cell: 1D cell global index [1,2,...,Ncell]
%         .ID: tensorized global intracellular ID index [M, 2, Njunction]
%         .boundary: one lumped outer-disc node [Nface, Ncell], zero if connected
%         .cleft: tensorized global cleft node index [M, Njunction]
% Inputs: cfg, graph topology, and ionic model.
% Output: network containing physical tensors and hidden solver operators.

Ncell = topology.Ncell;
Njunction = topology.Njunction;

mesh_index = cfg.junction_mesh(:);
if Njunction == 0
    mesh_index = zeros(0, 1);
elseif isscalar(mesh_index)
    mesh_index = repmat(mesh_index, Njunction, 1);
end
if numel(mesh_index) ~= Njunction
    error("junction_mesh must be scalar or contain Njunction values.");
end

scale_gj = cfg.scale_gj_loc(:);
scale_chan = cfg.scale_chan_loc(:);
if Njunction == 0
    scale_gj = zeros(0, 1);
    scale_chan = zeros(0, 1);
elseif isscalar(scale_gj)
    scale_gj = repmat(scale_gj, Njunction, 1);
end
if Njunction > 0 && isscalar(scale_chan)
    scale_chan = repmat(scale_chan, Njunction, 1);
end
if numel(scale_gj) ~= Njunction || numel(scale_chan) ~= Njunction
    error("Localization scales must be scalar or contain Njunction values.");
end

if Njunction == 0
    first_mesh_index = 1;
else
    first_mesh_index = mesh_index(1);
end
first_data = load(fullfile( ...
    cfg.mesh_folder, cfg.mesh_files(first_mesh_index)), "FEM_data");
M = numel(first_data.FEM_data.partition_surface);

mesh.Gc = zeros(M, M, Njunction);              % Local cleft-cleft conductance matrices.
mesh.Gb = zeros(M, Njunction);                 % Local cleft-to-bulk conductance.
mesh.area = zeros(M, Njunction);               % Membrane patch area for each local node.
mesh.volume = zeros(M, Njunction);             % Cleft volume for each local node.
mesh.gj_weight = zeros(M, Njunction);          % Local distribution weight of total GJ conductance.
mesh.current_weight = zeros(M, model.Ncurrents, Njunction); % Local distribution weight of each ionic current.

for n_junction = 1:Njunction
    loaded = load(fullfile( ...
        cfg.mesh_folder, cfg.mesh_files(mesh_index(n_junction))), "FEM_data");
    FEM = loaded.FEM_data;
    if numel(FEM.partition_surface) ~= M
        error("All junction meshes must have the same number of patches.");
    end

    mesh.Gc(:, :, n_junction) = cfg.f_cleft * ...
        FEM.cleft_adjacency_matrix / cfg.rho_ext;
    mesh.Gb(:, n_junction) = cfg.f_bulk * ...
        FEM.bulk_adjacency_matrix(:) / cfg.rho_ext;
    mesh.area(:, n_junction) = FEM.partition_surface(:);
    mesh.volume(:, n_junction) = cfg.fVol * FEM.partition_volume(:);

    gj_weight = FEM.gj_area_norm(:) .^ scale_gj(n_junction);
    mesh.gj_weight(:, n_junction) = gj_weight / sum(gj_weight);

    area_weight = FEM.partition_surface(:);
    area_weight = area_weight / sum(area_weight);
    current_weight = area_weight * model.loc_vec;

    if isfield(model.p, "iina")
        weight = FEM.Na_area_norm(:) .^ scale_chan(n_junction);
        current_weight(:, model.p.iina) = ...
            model.loc_vec(model.p.iina) * weight / sum(weight);
    end
    if isfield(model.p, "iinak")
        weight = FEM.NKA_area_norm(:);
        current_weight(:, model.p.iinak) = ...
            model.loc_vec(model.p.iinak) * weight / sum(weight);
    end
    if isfield(model.p, "iik1")
        weight = FEM.Kir21_area_norm(:);
        current_weight(:, model.p.iik1) = ...
            model.loc_vec(model.p.iik1) * weight / sum(weight);
    end
    mesh.current_weight(:, :, n_junction) = current_weight;
end

%% Stable tensor layout exposed to the simulation
layout.cell = reshape(1:Ncell, 1, Ncell);
next_node = Ncell;
layout.ID = reshape( ...
    next_node + (1:M * 2 * Njunction), M, 2, Njunction);
next_node = next_node + M * 2 * Njunction;

Nfaces = 2 * topology.dimension;
boundary_mask = true(Nfaces, Ncell);
for n_junction = 1:Njunction
    for side = 1:2
        cell_number = topology.junction_cells(side, n_junction);
        face_number = topology.junction_faces(side, n_junction);
        boundary_mask(face_number, cell_number) = false;
    end
end

layout.boundary = zeros(Nfaces, Ncell);
Nboundary = nnz(boundary_mask);
layout.boundary(boundary_mask) = next_node + (1:Nboundary);
next_node = next_node + Nboundary;

layout.cleft = reshape( ...
    next_node + (1:M * Njunction), M, Njunction);
Nnodes = next_node + M * Njunction;

%% Resistive edges
edge_i = zeros(0, 1); % edge_i(k): start node idx of edge k
edge_j = zeros(0, 1); % edge_i(k): end node idx of edge k
edge_g = zeros(0, 1); % edge_i(k): conductance on edge k

Rmyo = cfg.rho_myo * (cfg.L / 2) / (pi * cfg.r^2);
gmyo = 1 / Rmyo;
ggap = cfg.ggap;
if isempty(ggap)
    ggap = 7.35e-4 * cfg.D;
end

% Each unconnected cell face uses the original one-node terminal-disc model.
[boundary_face, boundary_cell] = find(boundary_mask);
edge_i = [edge_i; reshape(layout.cell(boundary_cell), [], 1)];
edge_j = [edge_j; layout.boundary(boundary_mask)];
edge_g = [edge_g; repmat(gmyo, Nboundary, 1)];

for n_junction = 1:Njunction
    for side = 1:2
        cell_number = topology.junction_cells(side, n_junction);
        edge_i = [edge_i; ...
            repmat(layout.cell(cell_number), M, 1)];
        edge_j = [edge_j; layout.ID(:, side, n_junction)];
        edge_g = [edge_g; repmat(gmyo / M, M, 1)];
    end

    edge_i = [edge_i; layout.ID(:, 1, n_junction)];
    edge_j = [edge_j; layout.ID(:, 2, n_junction)];
    edge_g = [edge_g; ggap * mesh.gj_weight(:, n_junction)];

    [row, column, value] = find(triu(mesh.Gc(:, :, n_junction), 1));
    for side = 1:2
        edge_i = [edge_i; layout.ID(row, side, n_junction)];
        edge_j = [edge_j; layout.ID(column, side, n_junction)];
        edge_g = [edge_g; value / cfg.rho_ie];
    end

    edge_i = [edge_i; layout.cleft(row, n_junction)];
    edge_j = [edge_j; layout.cleft(column, n_junction)];
    edge_g = [edge_g; value];

    bulk_patch = find(mesh.Gb(:, n_junction) > 0);
    edge_i = [edge_i; layout.cleft(bulk_patch, n_junction)];
    edge_j = [edge_j; zeros(numel(bulk_patch), 1)];
    edge_g = [edge_g; mesh.Gb(bulk_patch, n_junction)];
end

Nedge = numel(edge_g);
edge_number = (1:Nedge)';
non_ground = edge_j > 0;
incidence = sparse( ...
    [edge_number; edge_number(non_ground)], ...
    [edge_i; edge_j(non_ground)], ...
    [ones(Nedge, 1); -ones(nnz(non_ground), 1)], ...
    Nedge, Nnodes);
Gmat = incidence' * spdiags(edge_g, 0, Nedge, Nedge) * incidence;

%% Membrane patches and capacitance operator
patch.axial = reshape(1:Ncell, 1, Ncell);
patch.ID = reshape( ...
    Ncell + (1:M * 2 * Njunction), M, 2, Njunction);
patch.boundary = zeros(Nfaces, Ncell);
patch.boundary(boundary_mask) = ...
    Ncell + M * 2 * Njunction + (1:Nboundary);
Npatches = Ncell + M * 2 * Njunction + Nboundary;

outside_ID = zeros(M, 2, Njunction);
cleft_local = reshape(1:M * Njunction, M, Njunction);
patch_cleft = zeros(M, 2, Njunction);
for n_junction = 1:Njunction
    outside_ID(:, 1, n_junction) = layout.cleft(:, n_junction);
    outside_ID(:, 2, n_junction) = layout.cleft(:, n_junction);
    patch_cleft(:, 1, n_junction) = cleft_local(:, n_junction);
    patch_cleft(:, 2, n_junction) = cleft_local(:, n_junction);
end

inside_node = [ ...
    layout.cell(:); ...
    layout.ID(:); ...
    layout.boundary(boundary_mask)];
outside_node = [ ...
    zeros(Ncell, 1); ...
    outside_ID(:); ...
    zeros(Nboundary, 1)];
patch_number = (1:Npatches)';
has_outside_node = outside_node > 0;
membrane_incidence = sparse( ...
    [patch_number; patch_number(has_outside_node)], ...
    [inside_node; outside_node(has_outside_node)], ...
    [ones(Npatches, 1); -ones(nnz(has_outside_node), 1)], ...
    Npatches, Nnodes);

Aax = 2 * pi * cfg.r * cfg.L;
ID_area = zeros(M, 2, Njunction);
face_area = zeros(Nfaces, Ncell);
for n_junction = 1:Njunction
    for side = 1:2
        cell_number = topology.junction_cells(side, n_junction);
        face_number = topology.junction_faces(side, n_junction);
        ID_area(:, side, n_junction) = mesh.area(:, n_junction);
        face_area(face_number, cell_number) = sum(mesh.area(:, n_junction));
    end
end

reference_face_area = sum(first_data.FEM_data.partition_surface);
boundary_area = zeros(Nfaces, Ncell);
for boundary_number = 1:Nboundary
    face_number = boundary_face(boundary_number);
    cell_number = boundary_cell(boundary_number);
    if mod(face_number, 2) == 1
        opposite_face = face_number + 1;
    else
        opposite_face = face_number - 1;
    end
    boundary_area(face_number, cell_number) = ...
        face_area(opposite_face, cell_number);
    if boundary_area(face_number, cell_number) == 0
        boundary_area(face_number, cell_number) = reference_face_area;
    end
end

patch_capacitance = cfg.Cm * [ ...
    Aax * ones(Ncell, 1); ...
    ID_area(:); ...
    boundary_area(boundary_mask)];
Cmat = membrane_incidence' * ...
    spdiags(patch_capacitance, 0, Npatches, Npatches) * ...
    membrane_incidence;

%% Ionic-current fractions and patch ownership
f_I = zeros(model.Ncurrents, Npatches);
patch_cell = zeros(Npatches, 1);
patch_cell(patch.axial) = 1:Ncell;
f_I(:, patch.axial) = repmat(1 - model.loc_vec(:), 1, Ncell);

for n_junction = 1:Njunction
    for side = 1:2
        cell_number = topology.junction_cells(side, n_junction);
        current_weight = mesh.current_weight(:, :, n_junction) / Nfaces;
        f_I(:, patch.ID(:, side, n_junction)) = current_weight';
        patch_cell(patch.ID(:, side, n_junction)) = cell_number;
    end
end

for boundary_number = 1:Nboundary
    patch_number_boundary = ...
        patch.boundary(boundary_face(boundary_number), ...
        boundary_cell(boundary_number));
    f_I(:, patch_number_boundary) = model.loc_vec(:) / Nfaces;
    patch_cell(patch_number_boundary) = boundary_cell(boundary_number);
end

patch.cleft_linear = [ ...
    zeros(Ncell, 1); ...
    patch_cleft(:); ...
    zeros(Nboundary, 1)];
patch.cell = patch_cell;

network.Nnodes = Nnodes;
network.Npatches = Npatches;
network.M = M;
network.Nfaces = Nfaces;
network.boundary_mask = boundary_mask;
network.layout = layout;
network.patch = patch;
network.mesh = mesh;
network.G = Gmat;
network.Cm = Cmat;
network.Am = membrane_incidence;
network.f_I = f_I;
network.patch_capacitance = patch_capacitance;
network.edge = struct("from", edge_i, "to", edge_j, "g", edge_g);
end
