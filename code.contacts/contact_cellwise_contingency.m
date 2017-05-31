function [cell1_ids, cell2_ids, contingency_mat, contingency_mat_norm, cell1_types, cell2_types, ...
	distance_map, cells1_surface, cell1_n_edges] = ...
	contact_cellwise_contingency(cell_info, varargin)
	% varargin: types1, types2


nvarargin = length(varargin);
optargs = {1, ''};	 % 1 = class 1 = all GCs
optargs(1:nvarargin) = varargin;
[types1, types2] = optargs{:};

if isempty(types2)
	types2 = types1;
end


[~, surface_area, contingency_mat, cell_ids] = load_contact_vars();
surface_area = sortrows(surface_area.');	  % assuming the fact that build_cell_contingency() returns the cell ids sorted, hence the latter assert()
cell_ids = cell_ids(:);
surface_area = surface_area(ismember(surface_area(:,1), int64(cell_ids)), :);
assert(isequal(surface_area(:,1), cell_ids(:)))
% class(surface_area)
% class(contingency_mat)
surface_area = double(surface_area);

% get the union of the two groups in preparation for returning group assignments from a common pool of groups (types)
cells1 = get_cell_info(cell_info, types1);
cells2 = get_cell_info(cell_info, types2);
cell_info = get_cell_info(cell_info, unique([cells1.cell_id cells2.cell_id]));
cell_info_table = get_cell_info_table(cell_info);
cell_info_table.type = categorical(cell_info_table.type);

% TODO: this should be its own function
if ~isfield(cell_info, 'soma_coords_warped_mip2_zscaled')
	if ~isfield(cell_info, 'soma_coord_omni')
		cell_info_table.soma_coord = nan(size(cell_info_table,1), 3);
		cell_info = cell_info_get_soma_coord_omni(table2struct(cell_info_table));

		cell_info_table = struct2table(cell_info);
	elseif ~strcmp(class(cell_info_table.soma_coord), 'double')
		cell_info_table.soma_coord(cellfun(@isempty,cell_info_table.soma_coord)) = {nan(1,3)};
		% make the class double
		cell_info_table.soma_coord = cell2mat(cell_info_table.soma_coord);
	end

	soma_coords = cell_info_table.soma_coord;

	soma_coords = get_m2_warped_omni_coords(soma_coords);
	soma_coords = soma_coords(:, [2 3 1]);
	soma_coords(:, 2) = soma_coords(:, 2) * 23/16.5;

	cell_info_table.soma_coords_warped_mip2_zscaled = soma_coords;
elseif ~strcmp(class(cell_info_table.soma_coords_warped_mip2_zscaled), 'double')
	disp('soma_coords_warped_mip2_zscaled already computed but not all valid')
	cell_info_table.soma_coords_warped_mip2_zscaled(cellfun(@isempty,cell_info_table.soma_coords_warped_mip2_zscaled)) = {nan(1,3)};
	%class(cell_info_table.soma_coords_warped_mip2_zscaled)
	% make the class double
	cell_info_table.soma_coords_warped_mip2_zscaled = cell2mat(cell_info_table.soma_coords_warped_mip2_zscaled);
	%class(cell_info_table.soma_coords_warped_mip2_zscaled)
end
%cell_info = table2struct(cell_info_table);
cell_info = cell_info_table;


cells1 = get_cell_info(cell_info, [cells1.cell_id]);
idx1 = ismember(cell_ids, cells1.cell_id);
cell1_ids = cell_ids(idx1);
n1 = length(cell1_ids);
cells1 = get_cell_info(cell_info, cell1_ids);
cell1_types = cells1.type;


cells2 = get_cell_info(cell_info, [cells2.cell_id]);
idx2 = ismember(cell_ids, cells2.cell_id);
cell2_ids = cell_ids(idx2);
n2 = length(cell2_ids);
cells2 = get_cell_info(cell_info, cell2_ids);
cell2_types = cells2.type;


contingency_mat = contingency_mat(idx1, idx2);

cells1_surface = surface_area(idx1, 2);
%cells2_surface = surface_area(idx2, 2);

contingency_mat_norm = contingency_mat ./ repmat(cells1_surface(:), 1, n2);


contact_heat_map(contingency_mat, cell1_types, cell2_types);





cells1_soma = reshape(cells1.soma_coords_warped_mip2_zscaled, [], 1, 3);
cells2_soma = reshape(cells2.soma_coords_warped_mip2_zscaled, 1, [], 3);

%size(cells1_soma)
%size(cells2_soma)
distance_map = repmat(cells1_soma, 1, n2) - repmat(cells2_soma, n1, 1);
distance_map(:,:,3) = [];  % looking at 2D projection only, drop 3rd dimension
distance_map = sqrt(distance_map(:,:,1).^2 + distance_map(:,:,2).^2);

figure; 
mask = contingency_mat > 0;
nbins = 25;
nbins = 100;
nrows = 4;

contacthist = contingency_mat;
contacthist = contingency_mat_norm;
subplot(nrows,1,1:nrows-1)
histogram2(contacthist(mask), distance_map(mask), nbins);
xlabel('contact')
ylabel('soma 2D distance')
xlims = xlim();
try
title(varargin)
catch
end

subplot(nrows,1,nrows)
histogram(contacthist(mask), nbins);
xlabel('contact')
xlim(xlims)

% remove cells with 0 contact
%{
idx = any(contingency_mat);

cell_ids = cell_ids(idx);
contingency_mat = contingency_mat(idx, idx);
size(contingency_mat)

contingency_mat_raw = contingency_mat;
%}



