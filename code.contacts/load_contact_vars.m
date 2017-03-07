function [all, surfacearea, contingency_mat, cell_ids, surface_at_grid_keys, surface_at_grid_vals] = load_contact_vars()

basepath = 'contacts/raw3d_445-1167';

surfacearea_file_path = fullfile(basepath, 'surface_area.mat');
load(surfacearea_file_path);  % var: 'surfacearea'

%tic
contingency_file_path = fullfile(basepath, 'surface_contingency.mat');
load(contingency_file_path);  % var: confusion_mat
[contingency_mat, cell_ids] = build_cell_contingency(confusion_mat);
%toc

%tic
surface_grid_file_path = fullfile(basepath, 'surface_grid.mat');
load(surface_grid_file_path);  % var: surface_at_grid_keys, surface_at_grid_vals
%toc

	all.surfacearea = surfacearea;
	all.contingency_mat = contingency_mat;
	all.cell_ids = cell_ids;
	all.surface_at_grid_keys = surface_at_grid_keys;
	all.surface_at_grid_vals = surface_at_grid_vals;
