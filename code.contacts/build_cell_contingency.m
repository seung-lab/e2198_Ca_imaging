function [contingency_mat, cell1ids, cell2ids] = build_cell_contingency(confusion_mat, make_square_matrix)
	% confusion_mat: list of instances in the format of [id1 id2 value(s)]

% backward compatibility with the older call
if size(confusion_mat, 1) == 3 && size(confusion_mat, 2) ~= 3
	confusion_mat = confusion_mat.';
end

%{
basepath = 'contacts/raw3d_445-1167';

surfacearea_file_path = fullfile(basepath, 'surface_area.mat');
load(surfacearea_file_path);  % var: 'surfacearea'

contingency_file_path = fullfile(basepath, 'surface_contingency.mat');
load(contingency_file_path);  % var: confusion_mat

surface_grid_file_path = fullfile(basepath, 'surface_grid.mat');
load(surface_grid_file_path);  % var: surface_at_grid_keys, surface_at_grid_vals
%}

if ~exist('make_square_matrix', 'var')
	make_square_matrix = false;
end

if make_square_matrix
	[ids, ~, ind12] = unique(confusion_mat(:,1:2));
	cell1ids = ids;
	cell2ids = ids;
	ind12 = reshape(ind12, [], 2);
	ind1 = ind12(:,1);
	ind2 = ind12(:,2);
else
	[cell1ids, ~, ind1] = unique(confusion_mat(:,1));
	[cell2ids, ~, ind2] = unique(confusion_mat(:,2));
end

if nargout < 3
	assert(isequal(cell1ids, cell2ids))
	cell_ids = cell1ids;
end

n1 = length(cell1ids);
n2 = length(cell2ids);
d = size(confusion_mat,2)-2;

% build cross tabulation (giant matrix)
contingency_mat = zeros(n1, n2, d);
%contingency_mat = sparse(n,n);	% makes the next step slow, 45s
for k = 1:d
	kk = repmat(k, size(ind1));
	contingency_mat(sub2ind([n1 n2 d], ind1, ind2, kk)) = confusion_mat(:, 2+k);
end
%contingency_mat = sparse(contingency_mat); % 2.5s

%crosstab


