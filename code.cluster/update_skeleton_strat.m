function cell_info = update_skeleton_strat(cell_info)

load('skel_strat_real.20160620.mat', 'skel_strat_real')
skel_strat = skel_strat_real;
load('skel_strat.20160620.mat', 'skel_strat')
%ready = {'2i', '3i', '6t', '28', '58', '82n', '83', '91n'};

%[ready, ~, ~] = list_types(cell_info);
[~, ~, ~, ready] = list_types(cell_info);

to_update = get_cell_info(cell_info, ready);

binsize = to_update(1).strat_unrml(:,1);
binsize = abs(binsize(2) - binsize(1));

for id = [to_update.cell_id]
	idx = find(vertcat(cell_info.cell_id)==id);
	strat = skel_strat{id};
	%strat(:,1) = 100 - strat(:,1);
	if isempty(strat)
		warning(sprintf('not found: %d', id))
	end
	cell_info(idx).strat_unrml = strat;
	%strat(:, 2) = strat(:, 2) / sum(strat(:, 2)) / abs(strat(2,1) - strat(1,1));
	strat(:, 2) = strat(:, 2) / sum(strat(:, 2)) / binsize;
	cell_info(idx).strat_nrml = strat;
end

