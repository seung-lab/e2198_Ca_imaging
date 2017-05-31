function [all, surfacearea, contingency_mat, cell_ids, surface_at_grid_keys, surface_at_grid_vals, ...
	cell_hull] = load_contact_vars()
% warning: surfacearea may contain more cells than contingency_mat

basepath = 'contacts/raw3d_445-1167';

surfacearea_file_path = fullfile(basepath, 'surface_area.mat');
load(surfacearea_file_path);  % var: 'surfacearea'

%tic
contingency_file_path = fullfile(basepath, 'surface_contingency.mat');
load(contingency_file_path);  % var: confusion_mat
[contingency_mat, cell_ids] = build_cell_contingency(confusion_mat);
%toc

%{
ids = surfacearea(1,:);
ids(find(~ismember(ids , int64(cell_ids))))
%}

%tic
surface_grid_file_path = fullfile(basepath, 'surface_grid.mat');
load(surface_grid_file_path);  % var: surface_at_grid_keys, surface_at_grid_vals
%toc

if nargout>6 || nargout==1 % exist('cell_hull', 'var')
	convex_hulls_file_path = '~/dev/e2198-gc-analysis/cell_features/cell_hull.mat';
	load(convex_hulls_file_path);  % var: cell_hull
	all.cell_hull = cell_hull;
end

	all.surfacearea = surfacearea;
	all.contingency_mat = contingency_mat;
	all.cell_ids = cell_ids;
	all.surface_at_grid_keys = surface_at_grid_keys;
	all.surface_at_grid_vals = surface_at_grid_vals;
