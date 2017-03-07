function plot_all_contacts(stat_table, types1, types2, value_var, value_scale, cell_info)

% x: cell1_type, or value for individual cell2 cells
% y: cell2 or cell2_type

% filtering_types, selected_types, stat_type, grouping_types
%

allcelltypes = categories(stat_table.cell2_type);

if ~exist('value_var', 'var') || isempty(value_var)
	value_var = 'normalized_sum';
	normalized = 1;
elseif ~value_var
	value_var = 'sum_count';
	normalized = 0;
else
	normalized = 1;
end

if ~exist('value_scale', 'var') || isempty(value_scale)
	value_scale = 1;
end

yvalue = 'cell2_type';

individual_cells = false;
if exist('types2', 'var') && ~isempty(types2)
	types2 = types2(:).';

	if any(cellfun(@isnumeric, types2)) && exist('cell_info', 'var') % has individual cell id, use cell id.
		individual_cells = true;
		cells = get_cell_info(cell_info, types2);
		cells = [cells.cell_id];
		cells = cells(:).';
		idx = repmat(stat_table.cell2, size(cells)) == repmat(cells, size(stat_table.cell2));
	else % type names only, display as types
	idx = repmat(stat_table.cell2_type, size(types2)) == repmat(types2, size(stat_table.cell2_type));
	end

	idx = any(idx, 2);
	stat_table = stat_table(idx,:);
	size(stat_table)

	% this section only matters if showing types rather than individual_cells.
	only_show_selected_types = 1;
	if only_show_selected_types
		% remap the ordinal numbers in the categories
		stat_table.cell2_type = categorical(cellstr(stat_table.cell2_type));
		allcelltypes = categories(stat_table.cell2_type);
		if length(allcelltypes) ~= length(types2)
			warning('not all types are found')
		end
	end

	if individual_cells || length(types2)==1
		yvalue = 'cell2';
		only_show_selected_types = 1;
		stat_table.(yvalue) = categorical(cellstr(num2str(stat_table.(yvalue))));
	end

	diagline = 0;
else
	diagline = 1;
end
counts = stat_table;

lim = [min(double(counts.(yvalue)))-0.9, max(double(counts.(yvalue)))+0.9];


figure

if ~exist('types1', 'var') || isempty(types1)	% contacts by cell2_type
	%normalized = 1;
	if diagline
		%plot(lim, lim, 'Color', 0.88*[1 1 1]);%+normalized*0.1);
		plot(lim(1), lim(1))
		hold on
		plot3(lim, lim, [-1e10 -1e10], 'Color', 0.88*[1 1 1]);%+normalized*0.1);
	end
	% visual landmarks for easier comparison of two figs
	if 0
		ll = plot(repmat(4.5:5:65, 2, 1), lim, ':');%, 'Color', 0.88*[1 1 1]+normalized*0.1);
		for l = ll(:).'
			l.Color(4) = 0.3;
		end
		hold on
	end
	switch yvalue
	case {'cell2_type'}
		%nothing
	case {'cell2'}
		% highlight the type of itself
		ind = find(stat_table.cell1_type == stat_table.cell2_type(1), 1);
		ll = line(repmat(stat_table.cell1_type(ind), 1, 2), [0 100], 'Color', 'r', 'LineWidth', 0.1);
		ll.Color(4) = 0.1;
		hold on
	end


	if normalized
		counts.(value_var)(counts.(value_var)==0) = 1e-10;
	%scatter(counts.cell1_type, counts.(yvalue), counts.(value_var)*200*value_scale, counts.(value_var), 'filled')
	scatter3(counts.cell1_type, counts.(yvalue), -double(counts.(value_var)), counts.(value_var)*200*value_scale, counts.(value_var), 'filled')
	else
	%scatter(counts.cell1_type, counts.(yvalue), counts.sum_count/1000, counts.sum_count, 'filled')
	%scatter3(counts.cell1_type, counts.(yvalue), counts.sum_count, counts.sum_count/1000, counts.sum_count, 'filled')
	scatter3(counts.cell1_type, counts.(yvalue), counts.sum_count, counts.sum_count/100, counts.sum_count, 'filled')
	%scatter3(counts.cell1_type, counts.(yvalue), counts.sum_count, (counts.sum_count.^0.8)/10, counts.sum_count, 'filled')
	%surf(counts.cell1_type, counts.(yvalue), counts.sum_count)
	colorbar()
	end
	ax = gca;
	ax.XTick = 1:length(categories(counts.cell1_type));
	ax.XTickLabels = categories(counts.cell1_type);
	ax.XTickLabelRotation = 90;

else  	% contacts by cell, for the specific bc type

	%types1 = 'bc5i';
	cellcontacts = counts(counts.cell1_type == types1, :);
	scatter(cellcontacts.(value_var), cellcontacts.cell2_type, 30)
	title(types1);
	xlabel(types1);

end


ax = gca;

title(value_var, 'Interpreter', 'none')

switch yvalue
case {'cell2_type'}
	figure_size_x2([2 2])
	ax.YTick = 1:length(allcelltypes);
	ax.YTickLabels = allcelltypes;
case {'cell2'}
	figure_size_x2([2 1])
	ax.YTickLabels = categories(counts.(yvalue));
	ax.YTick = 1:length(ax.YTickLabels);
	%ylabel(allcelltypes)
	title(allcelltypes)
otherwise
	jalsdfklasd
end

%ylim([0, 1+length(allcelltypes)])
ylim(lim)
ax.GridAlpha = 0.05;
grid on
