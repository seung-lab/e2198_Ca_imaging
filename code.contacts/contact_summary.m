function [aggregate, individual, aggregate2, individual2] = contact_summary(cell_info, group1, group2, varargin) %direction, contact_vars)
% direction: 1/0 = normal, -1 = reverse, 2 = both directions (bi directional).

% Note I'm using the opposite convention as plot_all_contacts() and all_contacts() as to which is the 1 and the 2 in the pair,
% in the sense that 1 is the breakable group and 2 is the constant (aggregate) group here, where as in 
% plot_all_contacts() group 2 is the breakable group.
% TODO: probably fix plot_all_contacts() and all_contacts() to avoid further confusion, which will also make 
% all_contacts() consistent with gc_bc_contacts(). Looks like the reason I swapped the order when doing all_contacts()
% was only to have the default order for the returned cell stats table a convienent one.

nvarargin = length(varargin);
optargs = {1, [], []};
optargs(1:nvarargin) = varargin;
[direction, contact_vars, contactvoxels] = optargs{:};

if direction < 0 % inverse direction
	[aggregate, individual] = contact_summary(cell_info, group2, group1, 1, contact_vars);
	aggregate2 = []; individual2 = [];
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
	cell_hull = s.cell_hull;
else
	[~, surfacearea, contingency_mat, cell_ids, surface_at_grid_keys, surface_at_grid_vals, cell_hull] = ...
		load_contact_vars();
end

cells1 = get_cell_info(cell_info, group1); cells1 = [cells1.cell_id];
cells2 = get_cell_info(cell_info, group2); cells2 = [cells2.cell_id];

% debugging only
if 0 && isempty(contactvoxels)
contactvoxels = load_raw_contacts(cells1);
end

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
gridcounts1 = logical(surface1);
gridcounts1_sum = sum(gridcounts1, 3);
gridcounts1_any = any(gridcounts1, 3);
%gridcounts1_individualized

idx2 = look_up(surface_at_grid_keys, cells2);
surface2 = cat(3, surface_at_grid_vals{idx2});
surface2aggregated = sum(double(surface2), 3);	% summed over group2

%collocate = any(surface1, 3) & any(surface2, 3);	%mask

gridcounts2 = logical(surface2);
gridcounts2_sum = sum(gridcounts2, 3);  %gridcounts2aggregated
gridcounts2_any = any(gridcounts2, 3);
%gridcounts2_individualized
grids2total_sum = sum(gridcounts2_sum(:));
grids2total_any = sum(gridcounts2_any(:));


surface2total = sum(surface2aggregated(:));
maskedSurface2aggregate = surface2aggregated(any(surface1, 3));
maskedSurface2aggregate = sum(maskedSurface2aggregate(:));	% colocalized surface2aggregate

maskedSurface2 = [];
maskedSurface1 = [];
maskedGrids = [];
maskedGrids2sum = [];
hull_overlaps = [];	% for issubset==true only
hull_area = [];
hull_union_xx = [];
hull_union_yy = [];
hull_overlaps_union_xx = [];
hull_overlaps_union_yy = [];
grids1wContacts = [];
surf1inGrids1wContacts = [];
contactmaskAggregate = false(size(surface2aggregated));
for ii = 1:length(idx1)
	if issubset % need to recompute group2 by excluding the group1 cell
		newcells2 = setdiff(cells2,cells1(ii));
		idx2 = look_up(surface_at_grid_keys, newcells2);
		surface2 = cat(3, surface_at_grid_vals{idx2});
		surface2aggregated = sum(double(surface2), 3);	% summed over group2

		xx = []; yy=[];
		for jj = 1:length(newcells2)
			xxyy = flip(cell_hull{newcells2(jj)});
			[xx, yy] = polybool('union', xx, yy, xxyy(:,1), xxyy(:,2));
		end

		xxyy = flip(cell_hull{cells1(ii)});
		hull_area(end+1) = polyarea(xxyy(:,1), xxyy(:,2));
		[hull_union_xx, hull_union_yy] = polybool('union', xxyy(:,1), xxyy(:,2), hull_union_xx, hull_union_yy);

		[xx, yy] = polybool('intersection', xx, yy, xxyy(:,1), xxyy(:,2));
		hull_overlaps(end+1) = nanpolyarea(xx, yy);
		if isnan(hull_overlaps(end))
			xx
			yy
		end

		[hull_overlaps_union_xx, hull_overlaps_union_yy] = polybool('union', xx, yy, hull_overlaps_union_xx, hull_overlaps_union_yy);
	end
	masked = surface2aggregated(logical(surface1(:,:,ii)));
	% size(masked)  % apparently using a logical mask automatically makes it a vector..which makes sense
	maskedSurface2(end+1) = sum(masked(:));

	masked = gridcounts2_sum(logical(surface1(:,:,ii)));
	maskedGrids2sum(end+1) = sum(masked(:));

	masked = gridcounts1(:,:,ii);
	masked = masked(surface2aggregated>0);
	maskedGrids(end+1) = sum(masked(:));

	masked = surface1(:,:,ii);
	masked = masked(surface2aggregated>0);
	maskedSurface1(end+1) = sum(masked(:));

	% compute grids1wContacts and surf1inGrids1wContacts (2d grid)
	if ~isempty(contactvoxels)
		contactmask = false(size(contactmaskAggregate));
		vox = contactvoxels{cells1(ii)};
		vox = vox.';
		vox = vox(ismember(vox(:,1), cells2), 2:end);
		if 0  % debug. Arrrh we again ran into the issue of Int arithmatics with magical rounding in matlab
			class(vox)	%int32
			ceil((vox(:,1)-1)./15)
			ceil((vox(:,2)-1)./11)
		end
		vox = double(vox);

		%contactmask(sub2ind(size(contactmask), ceil((vox(:,1)-1)./15), ceil((vox(:,2)-1)./11))) = true;	% grid size was (15,11) with a 1 vox offset
		gridx = reshape(ceil((vox(:,1)-1)./15), 1, 1, []);	% grid size was (15,11) with a 1 vox offset
		gridy = reshape(ceil((vox(:,2)-1)./11), 1, 1, []);
		% expand, roughly equivalent to having larger grids. with parameter r (r=0: no expansion)
		n = length(gridx); sz = size(contactmask); r=0; d=2*r+1;
		gridx = repmat([-r:r], d, 1, n) + repmat(gridx, d, d, 1);
		gridy = repmat([-r:r].', 1, d, n) + repmat(gridy, d, d, 1);
		gridx(gridx<1) = 1; gridx(gridx>sz(1)) = sz(1);
		gridy(gridy<1) = 1; gridy(gridx>sz(2)) = sz(2);
		contactmask(sub2ind(size(contactmask), gridx(:), gridy(:))) = true;
		grids1wContacts(end+1) = sum(contactmask(:));	% to be used... need to individulize to compute the aggregate
		contactmaskAggregate = contactmaskAggregate | contactmask;
		
		masked = surface1(:,:,ii);
		if 0	% debug
			figure;imshow(masked)
		end
		masked = masked(contactmask);
		surf1inGrids1wContacts(end+1) = sum(masked(:));
	end
end

if 0	% debug
figure;imshow(contactmask)
figure;imshow(surface2aggregated>0)
end
% matrix form operations:
% tmp = repmat(surface2aggregated, 1, 1, length(idx1));
% maskedSurface2 = sum(sum(tmp(surface1), 1), 2)
maskedSurface1 = maskedSurface1(:);
maskedSurface2 = maskedSurface2(:);
maskedGrids = maskedGrids(:);
maskedGrids2sum = maskedGrids2sum(:);


hull_overlaps_union = nanpolyarea(hull_overlaps_union_xx, hull_overlaps_union_yy);
hull_union = nanpolyarea(hull_union_xx, hull_union_yy);
		if isnan(hull_union) || hull_union<0
			eeee
		end

stats = array2table(cells1(:), 'VariableNames', {'id'});
stats.contact = contacts;
stats.maskedSurface1 = maskedSurface1;
stats.maskedSurface2 = maskedSurface2;

surf1 = sum(sum(surface1, 1), 2);
surf1 = surf1(:);
assert(isequal(surf1, look_up(surfacearea.', cells1)))
stats.surface1 = surf1;
grids1 = sum(sum(gridcounts1, 1), 2);
grids1 = grids1(:);
stats.grids1 = grids1;
stats.maskedGrids = maskedGrids;
%size(maskedGrids2sum)
%size(stats)
stats.maskedGrids2sum = maskedGrids2sum;
if 0 && ~isempty(hull_overlaps)
stats.hull_overlaps = hull_overlaps(:);
end
stats.surf1inGrids1wContacts = surf1inGrids1wContacts(:);
stats.grids1wContacts = grids1wContacts(:);

stats = compute_percentage(stats, surface2total, grids2total_sum, grids2total_any, false);


aggregate = stats([], :);
aggregate(1, :) = {0};

aggregate.contact = sum(contacts);
aggregate.maskedSurface1 = sum(maskedSurface1);
aggregate.maskedGrids1sum = sum(maskedGrids);
aggregate.grids1 = sum(gridcounts1_any(:)); %aggregate.grids1_any
aggregate.grids1_sum = sum(grids1);
maskedGrids_mutualmin = min(gridcounts1_sum, gridcounts2_sum);
aggregate.surf1inGrids1wContacts = sum(surf1inGrids1wContacts(:));
aggregate.grids1wContacts = sum(contactmaskAggregate(:));
aggregate.grids1wContacts_sum = sum(grids1wContacts(:));
if issubset
	aggregate.maskedSurface2 = sum(maskedSurface2);  % in this case this one doesn't actually make sense and probably shouldn't be used
	aggregate.maskedGrids2sum = sum(maskedGrids2sum);  % ditto
	% same is true for the aggregate.contact value
	if length(subset) == length(cells2)
		% two sets are equal, make it make sense
		aggregate.contact = aggregate.contact; %/2;  % we've double counted the numbers, or did we?
		aggregate.maskedSurface2 = aggregate.maskedSurface1;
		aggregate.maskedGrids2sum = aggregate.maskedGrids1sum;
		aggregate.maskedGrids_mutualmin = aggregate.maskedGrids1sum;
	elseif length(subset)>1	 % we are just looking at individual cells wrt the type they belong to
		warning('bogus aggregate.maskedSurface2, aggregate.maskedGrids2sum, and aggregate.contact values')
	end
	tmp = maskedGrids_mutualmin > 1;	% at any location, a group1 cell itself plus one additional cell
	aggregate.maskedGrids = sum(tmp(:));
	aggregate.maskedGrids = sum(tmp(:));
%aggregate.grids1 
%hull_overlaps_union 
%hull_union
	%aggregate.grids1 = hull_overlaps_union * 16.5/23 / 15/11;
	%aggregate.grids1_sum = sum(hull_overlaps(:)) * 16.5/23 / 15/11;
	%Estimates on trying to only count grids in the convhull overlap regions
	aggregate.grids1 = aggregate.grids1 * hull_overlaps_union / hull_union;
	aggregate.grids1_sum = aggregate.grids1_sum * sum(hull_overlaps) / sum(hull_area);
else
	aggregate.maskedGrids_mutualmin = sum(maskedGrids_mutualmin(:));
	aggregate.maskedSurface2 = maskedSurface2aggregate;
	tmp = logical(maskedGrids_mutualmin);
	aggregate.maskedGrids = sum(tmp(:));
	if 1
	maskedGrids1_any = gridcounts1_any(surface2aggregated>0);
	maskedGrids1_any = sum(maskedGrids1_any(:));
	assert(aggregate.maskedGrids == maskedGrids1_any);
	end
	tmp = gridcounts2_sum(gridcounts1_any); %any(surface1, 3))
	aggregate.maskedGrids2sum = sum(tmp(:));
end
aggregate.surface1 = sum(surf1);

aggregate = compute_percentage(aggregate, surface2total, grids2total_sum, grids2total_any, true);
aggregate.surface2 = surface2total;
aggregate.grids2_any = grids2total_any;
aggregate.grids2_sum = grids2total_sum;


%stats.maskedSurface2 = % sumOverGroup2
%cat(3, surface1, surface2))
individual = stats;


function out = compute_percentage(out, surface2total, grids2total_sum, grids2total_any, isAggregate)
	out.surf1percent = out.maskedSurface1 ./ out.surface1;
	out.surf2percent = out.maskedSurface2 ./ surface2total;
	out.surf12percent = sqrt(out.surf1percent .* out.surf2percent);
	out.contact_surf1 = out.contact ./ out.surface1;
	out.contact_surf2 = out.contact ./ surface2total;
	out.contact_surf12 = out.contact ./ out.surface1 / surface2total;
	out.contact_maskedSurf1 = out.contact ./ out.maskedSurface1;
	out.contact_maskedSurf2 = out.contact ./ out.maskedSurface2;
	out.contact_maskedSurf12 = sqrt(out.contact_maskedSurf1 .* out.contact_maskedSurf2);
	out.grid1percent = out.maskedGrids ./ out.grids1;
	out.grid2percent = out.maskedGrids ./ grids2total_any; %out.grids2;
	out.grid12percent = sqrt(out.grid1percent .* out.grid2percent);
	out.grid2percent_lb = out.maskedGrids ./ grids2total_sum; %out.grids2_sum;
	out.grid2percent_ub = out.maskedGrids2sum ./ grids2total_sum; %out.grids2_sum;
	
	out.fasciculation = out.contact ./ out.surf1inGrids1wContacts; % surf2inGrids1wContacts
	out.grids1wContacts_maskedGrids = out.grids1wContacts ./ out.maskedGrids;
	% We don't have surf2inGrids2wContacts
	%out.surfInGridsWContacts_surfInMaskedGrids = sqrt(out.surf1inGrids1wContacts .* surf2inGrids2wContacts ./ out.maskedSurface1 ./ out.maskedSurface2);

	if isAggregate
		out.grid1percent_ub = out.maskedGrids1sum ./ out.grids1_sum;
		out.grid1percenta = out.maskedGrids_mutualmin ./ out.grids1_sum;
		out.grid2percenta = out.maskedGrids_mutualmin ./ grids2total_sum; %out.grids2_sum;

		out.grids1wContacts_maskedGrids1_sum = out.grids1wContacts_sum ./ out.maskedGrids1sum;
	end

function a = nanpolyarea2(x, y)	% hopefully there's no self intersections
	% confirmed this produces "numerically identical" (i.e. identical other than numerical error) results as this:
	% https://www.mathworks.com/matlabcentral/fileexchange/49179-nanpolyarea-x-y-/content/nanpolyarea.m
	% , at least for our data used here.
	[x, y] = polysplit(x, y);
	a = 0;
	for ii = 1:length(x)
		if ispolycw(x{ii}, y{ii})
			a = a + polyarea(x{ii}, y{ii});
		else
			a = a - polyarea(x{ii}, y{ii});
		end
	end
