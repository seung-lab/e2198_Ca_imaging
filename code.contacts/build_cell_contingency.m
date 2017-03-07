function [contingency_mat, cell_ids] = build_cell_contingency(confusion_mat)
	% confusion_mat: list of instances

%{
basepath = 'contacts/raw3d_445-1167';

surfacearea_file_path = fullfile(basepath, 'surface_area.mat');
load(surfacearea_file_path);  % var: 'surfacearea'

contingency_file_path = fullfile(basepath, 'surface_contingency.mat');
load(contingency_file_path);  % var: confusion_mat

surface_grid_file_path = fullfile(basepath, 'surface_grid.mat');
load(surface_grid_file_path);  % var: surface_at_grid_keys, surface_at_grid_vals
%}

[cell1ids, ~, ind1] = unique(confusion_mat(1,:));
[cell2ids, ~, ind2] = unique(confusion_mat(2,:));
assert(isequal(cell1ids, cell2ids))
cell_ids = cell1ids;
n = length(cell_ids);

% build cross tabulation (giant matrix)
contingency_mat = zeros(n);
%contingency_mat = sparse(n,n);	% makes the next step slow, 45s
contingency_mat(sub2ind([n n], ind1, ind2)) = confusion_mat(3,:);
%contingency_mat = sparse(contingency_mat); % 2.5s

%crosstab


