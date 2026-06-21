function network = assemble_continuum_network(cfg, topology, model)
%ASSEMBLE_CONTINUUM_NETWORK Build layouts and physical conductance tensors.
%
% This file is intentionally the "geometry and bookkeeping" layer. It does
% not advance time and it does not evaluate ionic model dynamics. Its job is
% to translate the user-facing inputs
%   adjacency_matrix + FEM meshes + channel localization settings
% into reusable tensors/operators for the continuum-note variables:
%   phi_m, phi_pre, phi_post, phi_c, and s_c.
%
% Three sparse operators appear here:
%   G   : conductance/Laplacian operator for Ohmic graph currents.
%   Cm  : capacitance mass matrix, built from membrane patch capacitances.
%   Am  : membrane incidence/projection matrix.
%
% Am has two roles:
%   Vm = Am * phi      maps node potentials to patch membrane voltages.
%   Am' * Iion        projects patch ionic currents back into node KCL rows.
%
% The time loop keeps physical variables separated. It stacks them only at
% the sparse-solve boundary because G, Cm, and Am are most naturally sparse
% matrices over global electrical nodes.

Ncell = topology.Ncell;
Njunction = topology.Njunction;

%% Per-junction mesh and localization settings
% A topology can have many junctions. cfg.junction_mesh, cfg.scale_gj_loc,
% and cfg.scale_chan_loc may be scalars, in which case the same mesh or
% localization exponent is reused for every junction.
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

%% Load the FEM meshes once
% Each FEM_data file provides patch-level geometry and measured/constructed
% adjacency data on one ID cleft:
%   partition_surface       -> membrane area of each ID patch
%   partition_volume        -> local cleft volume used in ds_c/dt
%   cleft_adjacency_matrix  -> lateral cleft conductance graph
%   bulk_adjacency_matrix   -> cleft-to-bulk conductance per patch
%   *_area_norm             -> channel/GJ localization weights
if Njunction == 0
    ref = load(fullfile(cfg.mesh_folder, cfg.mesh_files(1)), "FEM_data");
    loaded_meshes = {ref.FEM_data};
    file_map = zeros(0, 1);
    M = numel(ref.FEM_data.partition_surface);
    reference_face_area = sum(ref.FEM_data.partition_surface);
else
    files_for_junction = cfg.mesh_files(mesh_index);
    [unique_files, ~, file_map] = unique(files_for_junction(:), "stable");
    loaded_meshes = cell(numel(unique_files), 1);
    for u = 1:numel(unique_files)
        d = load(fullfile(cfg.mesh_folder, unique_files(u)), "FEM_data");
        loaded_meshes{u} = d.FEM_data;
    end
    M = numel(loaded_meshes{file_map(1)}.partition_surface);
    reference_face_area = sum(loaded_meshes{file_map(1)}.partition_surface);
end

%% FEM data -> physical tensors on each cleft surface
% mesh.Gc is the discrete version of the cleft surface operator:
%   div_Omega(sigma_c grad_Omega eta_c).
% It is kept as pairwise conductances because the code evaluates fluxes by
% patch-to-patch differences.
%
% mesh.Gb is the local bulk exchange conductance bar_g_b.
% mesh.current_weight tells the ionic model how much of each current density
% is assigned to each ID patch. This is where "most Na current lives in the
% cleft/ID" enters numerically through Na_area_norm and cfg.locINa.
mesh.Gc = zeros(M, M, Njunction);
mesh.Gb = zeros(M, Njunction);
mesh.area = zeros(M, Njunction);
mesh.volume = zeros(M, Njunction);
mesh.gj_weight = zeros(M, Njunction);
mesh.current_weight = zeros(M, model.Ncurrents, Njunction);

for n_junction = 1:Njunction
    FEM = loaded_meshes{file_map(n_junction)};
    if numel(FEM.partition_surface) ~= M
        error("All junction meshes must have the same number of patches.");
    end

    % Convert FEM adjacency/resistance data into electrical conductances.
    % rho_ext scales extracellular cleft and bulk paths.
    mesh.Gc(:, :, n_junction) = cfg.f_cleft * ...
        FEM.cleft_adjacency_matrix / cfg.rho_ext;
    mesh.Gb(:, n_junction) = cfg.f_bulk * ...
        FEM.bulk_adjacency_matrix(:) / cfg.rho_ext;

    % Area affects membrane capacitance and ID current distribution.
    % Volume affects the cleft concentration equation ds_c/dt.
    mesh.area(:, n_junction) = FEM.partition_surface(:);
    mesh.volume(:, n_junction) = cfg.fVol * FEM.partition_volume(:);

    % Gap junction conductance is total conductance times a local weight.
    % scale_gj_loc controls how strongly GJ conductance follows gj_area_norm.
    gj_weight = FEM.gj_area_norm(:) .^ scale_gj(n_junction);
    mesh.gj_weight(:, n_junction) = gj_weight / sum(gj_weight);

    % Default: distribute localized current by membrane area.
    % Selected currents override this with experimentally motivated channel
    % maps: Na_area_norm, NKA_area_norm, and Kir21_area_norm.
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

%% Global electrical-node layout
% The continuum code thinks in separated fields:
%   phi_m      Ncell-by-1
%   phi_pre    M-by-Njunction
%   phi_post   M-by-Njunction
%   phi_c      M-by-Njunction
%   phi_boundary for uncoupled terminal faces
%
% Sparse linear algebra still needs one global ordering. layout records that
% ordering. These indices are bookkeeping only; they are not the conceptual
% state representation used by the run loop.
layout.cell = reshape(1:Ncell, 1, Ncell);
next_node = Ncell;
layout.ID = reshape( ...
    next_node + (1:M * 2 * Njunction), M, 2, Njunction);
next_node = next_node + M * 2 * Njunction;

Nports = topology.Nports;
boundary_mask = topology.boundary_mask;
layout.boundary = zeros(Nports, Ncell);
Nboundary = nnz(boundary_mask);
layout.boundary(boundary_mask) = next_node + (1:Nboundary);
next_node = next_node + Nboundary;

layout.cleft = reshape( ...
    next_node + (1:M * Njunction), M, Njunction);
Nnodes = next_node + M * Njunction;

%% Conductive graph edges -> G matrix
% Every edge represents an Ohmic current g*(phi_i - phi_j). The incidence
% matrix B gives edge voltage drops B*phi, and
%   G = B' * diag(g) * B
% is the global graph Laplacian/conductance operator. This is the discrete
% counterpart of the div(sigma grad phi) terms in the note.
Rmyo = cfg.rho_myo * (cfg.L / 2) / (pi * cfg.r^2);
gmyo = 1 / Rmyo;
ggap = cfg.ggap;
if isempty(ggap)
    ggap = 7.35e-4 * cfg.D;
end

seg_i = cell(1, 1 + 7 * Njunction);
seg_j = cell(1, 1 + 7 * Njunction);
seg_g = cell(1, 1 + 7 * Njunction);
k = 0;

% Boundary ports are terminal ID faces with no neighboring cell. They keep
% the same one-node terminal-disc model used by the older 1-D code.
[boundary_face, boundary_cell] = find(boundary_mask);
k = k + 1;
seg_i{k} = reshape(layout.cell(boundary_cell), [], 1);
seg_j{k} = layout.boundary(boundary_mask);
seg_g{k} = repmat(gmyo, Nboundary, 1);

for n_junction = 1:Njunction
    % Axial/myoplasm-to-ID coupling:
    %   gmyo/M * (phi_m(cell) - phi_pre/post(patch)).
    % This corresponds to the g_myo(phi_m - <phi_pre/post>) term after
    % summing over all patches on a junction side.
    for side = 1:2
        cell_number = topology.junction_cells(side, n_junction);
        k = k + 1;
        seg_i{k} = repmat(layout.cell(cell_number), M, 1);
        seg_j{k} = layout.ID(:, side, n_junction);
        seg_g{k} = repmat(gmyo / M, M, 1);
    end

    % Patchwise gap-junction coupling phi_pre <-> phi_post.
    % The total GJ conductance is distributed over patches by gj_weight.
    k = k + 1;
    seg_i{k} = layout.ID(:, 1, n_junction);
    seg_j{k} = layout.ID(:, 2, n_junction);
    seg_g{k} = ggap * mesh.gj_weight(:, n_junction);

    % Lateral intracellular ID conduction on each side of the disc.
    % The same FEM cleft adjacency is reused as a patch graph, scaled by
    % rho_ie to allow ID membrane-side resistance to differ from cleft
    % extracellular resistance.
    [row, column, value] = find(triu(mesh.Gc(:, :, n_junction), 1));
    for side = 1:2
        k = k + 1;
        seg_i{k} = layout.ID(row, side, n_junction);
        seg_j{k} = layout.ID(column, side, n_junction);
        seg_g{k} = value / cfg.rho_ie;
    end

    % Lateral cleft conduction between cleft patches.
    k = k + 1;
    seg_i{k} = layout.cleft(row, n_junction);
    seg_j{k} = layout.cleft(column, n_junction);
    seg_g{k} = value;

    % Cleft-to-bulk exchange. A zero endpoint means electrical ground/bulk
    % reference, so only the cleft node appears in the finite node vector.
    bulk_patch = find(mesh.Gb(:, n_junction) > 0);
    k = k + 1;
    seg_i{k} = layout.cleft(bulk_patch, n_junction);
    seg_j{k} = zeros(numel(bulk_patch), 1);
    seg_g{k} = mesh.Gb(bulk_patch, n_junction);
end

edge_i = vertcat(seg_i{1:k});
edge_j = vertcat(seg_j{1:k});
edge_g = vertcat(seg_g{1:k});

Nedge = numel(edge_g);
edge_number = (1:Nedge)';
non_ground = edge_j > 0;
incidence = sparse( ...
    [edge_number; edge_number(non_ground)], ...
    [edge_i; edge_j(non_ground)], ...
    [ones(Nedge, 1); -ones(nnz(non_ground), 1)], ...
    Nedge, Nnodes);
Gmat = incidence' * spdiags(edge_g, 0, Nedge, Nedge) * incidence;

%% Membrane patch layout
% Electrical nodes and membrane patches are different objects:
%   node  = unknown electrical potential
%   patch = membrane surface area with capacitance and ionic currents
%
% patch.* gives the patch ordering expected by the ionic model. It mirrors
% the node layout where possible, but boundary patches and cleft mapping are
% patch-level concepts.
patch.axial = reshape(1:Ncell, 1, Ncell);
patch.ID = reshape( ...
    Ncell + (1:M * 2 * Njunction), M, 2, Njunction);
patch.boundary = zeros(Nports, Ncell);
patch.boundary(boundary_mask) = ...
    Ncell + M * 2 * Njunction + (1:Nboundary);
Npatches = Ncell + M * 2 * Njunction + Nboundary;

% For ID membrane patches, the outside node is the local cleft potential.
% For axial/boundary membrane patches, the outside is the bulk/reference
% potential and is represented by zero in outside_node.
outside_ID = zeros(M, 2, Njunction);
cleft_local = reshape(1:M * Njunction, M, Njunction);
patch_cleft = zeros(M, 2, Njunction);
for n_junction = 1:Njunction
    outside_ID(:, 1, n_junction) = layout.cleft(:, n_junction);
    outside_ID(:, 2, n_junction) = layout.cleft(:, n_junction);
    patch_cleft(:, 1, n_junction) = cleft_local(:, n_junction);
    patch_cleft(:, 2, n_junction) = cleft_local(:, n_junction);
end

%% Membrane incidence Am and capacitance matrix Cm
% Each membrane patch contributes a voltage difference:
%   V_patch = phi_inside - phi_outside.
% membrane_incidence is the matrix Am that encodes these signs.
%
% This is why Am appears twice in the time loop:
%   V_membrane = Am * phi
%       gives the voltage passed to the ionic model.
%   -Am' * Iion
%       injects equal-and-opposite ionic currents into the inside/outside
%       electrical nodes for the KCL equations.
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

% Patch capacitance is capacitance density times physical membrane area.
% Axial patches use the cylindrical side area. ID patches use the FEM patch
% areas. Unconnected boundary ports reuse a representative ID face area.
Aax = 2 * pi * cfg.r * cfg.L;
ID_area = zeros(M, 2, Njunction);
port_area = zeros(Nports, Ncell);
for n_junction = 1:Njunction
    for side = 1:2
        cell_number = topology.junction_cells(side, n_junction);
        port_number = topology.junction_ports(side, n_junction);
        ID_area(:, side, n_junction) = mesh.area(:, n_junction);
        port_area(port_number, cell_number) = sum(mesh.area(:, n_junction));
    end
end

boundary_area = zeros(Nports, Ncell);
for cell_number = 1:Ncell
    connected_area = port_area(port_area(:, cell_number) > 0, cell_number);
    if isempty(connected_area)
        local_boundary_area = reference_face_area;
    else
        local_boundary_area = mean(connected_area);
    end
    boundary_area(boundary_mask(:, cell_number), cell_number) = ...
        local_boundary_area;
end

patch_capacitance = cfg.Cm * [ ...
    Aax * ones(Ncell, 1); ...
    ID_area(:); ...
    boundary_area(boundary_mask)];

% Cm is the mass matrix in the semi-discrete electrical equation:
%   Cm * dphi/dt = -G*phi - Am'*Iion + source.
% It is generally singular because some rows are algebraic constraints
% rather than independent ODE states.
Cmat = membrane_incidence' * ...
    spdiags(patch_capacitance, 0, Npatches, Npatches) * ...
    membrane_incidence;

%% Ionic-current fractions and patch ownership
% f_I(current, patch) tells the ionic model what fraction of a cell's total
% current density belongs on each membrane patch.
%
% Axial patches receive the non-ID fraction 1 - loc_vec.
% ID patches receive mesh.current_weight, split over however many membrane
% ports that cell has. Boundary terminal discs receive the remaining ID
% fraction for ports that are not connected to another cell.
f_I = zeros(model.Ncurrents, Npatches);
patch_cell = zeros(Npatches, 1);
patch_cell(patch.axial) = 1:Ncell;
f_I(:, patch.axial) = repmat(1 - model.loc_vec(:), 1, Ncell);

for n_junction = 1:Njunction
    for side = 1:2
        cell_number = topology.junction_cells(side, n_junction);
        current_weight = mesh.current_weight(:, :, n_junction) / ...
            topology.cell_port_count(cell_number);
        f_I(:, patch.ID(:, side, n_junction)) = current_weight';
        patch_cell(patch.ID(:, side, n_junction)) = cell_number;
    end
end

for boundary_number = 1:Nboundary
    patch_number_boundary = ...
        patch.boundary(boundary_face(boundary_number), ...
        boundary_cell(boundary_number));
    cell_number = boundary_cell(boundary_number);
    f_I(:, patch_number_boundary) = model.loc_vec(:) / ...
        topology.cell_port_count(cell_number);
    patch_cell(patch_number_boundary) = cell_number;
end

patch.cleft_linear = [ ...
    zeros(Ncell, 1); ...
    patch_cleft(:); ...
    zeros(Nboundary, 1)];
patch.cell = patch_cell;

network.Nnodes = Nnodes;
network.Npatches = Npatches;
network.M = M;
network.Nports = Nports;
network.Nfaces = Nports;
network.boundary_mask = boundary_mask;
network.layout = layout;
network.patch = patch;
network.mesh = mesh;
network.G = Gmat;
network.Cm = Cmat;
network.gmyo = gmyo;
network.ggap = ggap;
network.Am = membrane_incidence;
network.f_I = f_I;
network.patch_capacitance = patch_capacitance;
network.edge = struct("from", edge_i, "to", edge_j, "g", edge_g);
end
