function [aggregate, individual, aggregate2, individual2] = contact_summary(cell_info, group1, group2, varargin) %direction, contact_vars)
% direction: 1/0 = normal, -1 = reverse, 2 = both directions (bi directional).

nvarargin = length(varargin);
optargs = {1, []};
optargs(1:nvarargin) = varargin;
[direction, contact_vars] = optargs{:};

if direction < 0 % inverse direction
	[aggregate, individual] = contact_summary(cell_info, group2, group1, 1, contact_vars);
	return;
elseif direction > 1 && nargout > 2 
	[aggregate2, individual2] = contact_summary(cell_info, group2, group1, 1, contact_vars);
else
	aggregate2 = []; individual2 = [];
end

if isstruct(contact_vars)
	s = contact_vars;
	surfacearea = s.surfacearea;
	contingency_mat = s.contingency_mat;
	cell_ids = s.cell_ids;
	surface_at_grid_keys = s.surface_at_grid_keys;
	surface_at_grid_vals = s.surface_at_grid_vals;
else
	[~, surfacearea, contingency_mat, cell_ids, surface_at_grid_keys, surface_at_grid_vals] = ...
		load_contact_vars();
end

cells1 = get_cell_info(cell_info, group1); cells1 = [cells1.cell_id];
cells2 = get_cell_info(cell_info, group2); cells2 = [cells2.cell_id];

subset = intersect(cells1, cells2);
if length(subset) == length(cells1)
	issubset = true;	% need to exclude each cell when computing the mask for group2
	if length(cells1) == 1
		% one cell to itself, makes no sense, return empty
		display(['ignored: one cell to itself: ' mat2str(cells1)])  % R2016b: string(group1)
		aggregate = []; individual = [];
		return;
	end
elseif length(subset)>0
	error(['Two cell groups have common members\n' mat2str(cells1) '\n' mat2str(cells2)])  % R2016b: string(group1)
else
	issubset = false;
end

[~,~,contacts] = find_group_contingency(contingency_mat, cell_ids, cells1, cells2);
%individual = contingency_mat(ismember2(cells1, contingency_labels), idx2);

%total = sum(individual(:));

%surface_a = sum(look_up(surfacearea.', cells1));

idx1 = look_up(surface_at_grid_keys, cells1);
surface1 = cat(3, surface_at_grid_vals{idx1});
surface1 = double(surface1);
%surface1 = sum(double(surface1), 3);

idx2 = look_up(surface_at_grid_keys, cells2);
surface2 = cat(3, surface_at_grid_vals{idx2});
surface2aggregated = sum(double(surface2), 3);	% summed over group2

%collocate = any(surface1, 3) & any(surface2, 3);	%mask

surface2total = sum(surface2aggregated(:));
maskedSurface2aggregate = surface2aggregated(any(surface1, 3));
maskedSurface2aggregate = sum(maskedSurface2aggregate(:));	% colocalized surface2aggregate

maskedSurface2 = [];
maskedSurface1 = [];
for ii = 1:length(idx1)
	if issubset % need to recompute group2 by excluding the group1 cell
		newcells2 = setdiff(cells2,cells1(ii));
		idx2 = look_up(surface_at_grid_keys, newcells2);
		surface2 = cat(3, surface_at_grid_vals{idx2});
		surface2aggregated = sum(double(surface2), 3);	% summed over group2
	end
	masked = surface2aggregated(logical(surface1(:,:,ii)));
	% size(masked)  % apparently using a logical mask automatically makes it a vector..which makes sense
	maskedSurface2(end+1) = sum(masked(:));

	masked = surface1(:,:,ii);
	masked = masked(surface2aggregated>0);
	maskedSurface1(end+1) = sum(masked(:));
end
% matrix form operations:
% tmp = repmat(surface2aggregated, 1, 1, length(idx1));
% maskedSurface2 = sum(sum(tmp(surface1), 1), 2)
maskedSurface1 = maskedSurface1(:);
maskedSurface2 = maskedSurface2(:);

stats = array2table(cells1(:), 'VariableNames', {'id'});
stats.contact = contacts;
stats.maskedSurface1 = maskedSurface1;
stats.maskedSurface2 = maskedSurface2;

surf1 = sum(sum(surface1, 1), 2);
surf1 = surf1(:);
assert(isequal(surf1, look_up(surfacearea.', cells1)))
stats.surface1 = surf1;
stats = compute_percentage(stats, surface2total);


aggregate = stats([], :);
aggregate(1, :) = {0};

aggregate.contact = sum(contacts);
aggregate.maskedSurface1 = sum(maskedSurface1);
if issubset
	aggregate.maskedSurface2 = sum(maskedSurface2);  % in this case this one doesn't actually make sense and probably shouldn't be used
	% same is true for the aggregate.contact value
	if length(subset) == length(cells2)
		% two sets are equal, make it make sense
		aggregate.contact = aggregate.contact/2;  % we've double counted the numbers
		aggregate.maskedSurface2 = aggregate.maskedSurface1;
	elseif length(subset)>1	 % we are just look at a single cell wrt its own type
		warning('bogus aggregate.maskedSurface2 and aggregate.contact values')
	end
else
	aggregate.maskedSurface2 = maskedSurface2aggregate;
end
aggregate.surface1 = sum(surf1);

aggregate = compute_percentage(aggregate, surface2total);
aggregate.surface2 = surface2total;


%stats.maskedSurface2 = % sumOverGroup2
%cat(3, surface1, surface2))
individual = stats;


function out = compute_percentage(out, surface2total)
	out.surf1percent = out.maskedSurface1 ./ out.surface1;
	out.surf2percent = out.maskedSurface2 ./ surface2total;
	out.contact_surf1 = out.contact ./ out.surface1;
	out.contact_surf2 = out.contact ./ surface2total;
	out.contact_maskedSurf1 = out.contact ./ out.maskedSurface1;
	out.contact_maskedSurf2 = out.contact ./ out.maskedSurface2;
