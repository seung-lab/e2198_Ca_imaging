function filtered_stats = summarize_intratype_contacts(stat_table, cell_info)

% cell_info should be the cell_info_new returned by build_contingency_stats

%{
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

	idx = repmat(stat_table.cell2_type, size(types2)) == repmat(types2, size(stat_table.cell2_type));
	idx = any(idx, 2);
	stat_table = stat_table(idx,:);
	size(stat_table)
%}

filtered_stats = stat_table(stat_table.cell1_type==stat_table.cell2_type,:);

%writetable(summarize_intratype_contacts(type_stats_w_gridcounts2), 'contacts/intratype.csv')
return
outpath = 'contacts/intratype/';
%for ctype = categories(stat_table.cell2_type).'	% this should be cell array of char strings as of R2016a
%	ctype = ctype{1};
for ctype = filtered_stats.cell2_type.'
	ctype = char(ctype);
	display(ctype)
	contact_investigate(cell_info, ctype,  'self');
	fig = gcf();
	saveas(fig, [outpath ctype '.png'])
	close(fig);
end
