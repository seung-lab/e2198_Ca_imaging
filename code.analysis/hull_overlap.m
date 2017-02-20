function ratio = hull_overlap(cell_info, cellquery, against_type)

if ~exist('against_type', 'var')
	against_type = [];
end

hull_file = '~/dev/e2198-gc-analysis/cell_features/cell_hull.mat';  % cell_hull{cell_id}
load(hull_file);

cells = get_cell_info_table(cell_info, cellquery);

if isempty(against_type) || isequal(cellquery, against_type)  % intra type cell overlaps, similar to coverage factor
	x = [];
	y = [];
	area_sum = 0;
	for cid = cells.cell_id(:).'
		hull = cell_hull{cid};
		hull = hull(end:-1:1, :); % convert to clock-wise convention for polybool()
		[x, y] = polybool('union', x, y, hull(:,1),hull(:,2));
		area_sum = area_sum + polyarea(hull(:,1),hull(:,2));
	end
	ratio = area_sum / polyarea(x, y) - 1;
else
	body
end
