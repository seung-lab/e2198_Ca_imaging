function plot_strat(cell_info,query, average, use_normalized_strat)


if ~exist('use_normalized_strat', 'var')
    use_normalized_strat = true;
end

if use_normalized_strat
    stratname = 'strat_nrml';
else
    stratname = 'strat_unrml';
end


x=cell_info(1).(stratname)(:,1);

ngroups = numel(query);
colors=distinguishable_colors(numel(query));

if ~exist('average', 'var')
	average = false;
end

figure; hold on;
legends = {};
entries = [];
ncells = 0;
for k = 1:ngroups
	group = query{k};
    cells = get_cell_info(cell_info, group);
    strats = cat(3, cells.(stratname));
    strats = squeeze(strats(:, 2, :));
    if average
    	strats = mean(strats, 2);
    end
    if ngroups > 1 || average
        h = plot(x, strats, 'Color', colors(mod(k-1,size(colors,1))+1,:), 'LineWidth', 2);
        %ncells = ncells + size(strats,2);
        %entries(end+1) = ncells;
        entries(end+1) = h(end);
        legends{end+1} = evalc('disp(group)');
    else
        h = plot(x, strats);
        entries = h;
        legends = cellstr(num2str(vertcat(cells.cell_id)));
    end
end
xlim([-20, 120])
legend(entries, legends);

