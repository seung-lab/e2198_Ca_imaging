function plot_bc_contacts(stat_table, bctype, celltypes)

allcelltypes = categories(stat_table.celltype);


if exist('celltypes', 'var') && ~isempty(celltypes)
	celltypes = celltypes(:).';
	idx = repmat(stat_table.celltype, size(celltypes)) == repmat(celltypes, size(stat_table.celltype));
	idx = any(idx, 2);
	stat_table = stat_table(idx,:);
	size(stat_table)

	only_show_selected_types = 1;
	if only_show_selected_types
		% remap the ordinal numbers in the categories
		stat_table.celltype = categorical(cellstr(stat_table.celltype));
		allcelltypes = categories(stat_table.celltype);
		if length(allcelltypes) ~= length(celltypes)
			warning('not all types are found')
		end
	end
end

counts = stat_table;

figure

if ~exist('bctype', 'var') || isempty(bctype)	% contacts by celltype
	scatter(counts.bctype, counts.celltype, counts.normalized_sum*100, counts.normalized_sum, 'filled')
	ax = gca;
	ax.XTick = 1:length(categories(counts.bctype));
	ax.XTickLabels = categories(counts.bctype);
	ax.XTickLabelRotation = 30;

else  	% contacts by cell, for the specific bc type

	%bctype = 'bc5i';
	cellcontacts = counts(counts.bctype == bctype, :);
	scatter(cellcontacts.normalized_sum, cellcontacts.celltype, 30)
	title(bctype);
	xlabel(bctype);

end

ax = gca;


ax.YTick = 1:length(allcelltypes);
ax.YTickLabels = allcelltypes;
%ylim([0, 1+length(allcelltypes)])
ylim([min(double(counts.celltype))-0.9, max(double(counts.celltype))+0.9])
ax.GridAlpha = 0.05;
grid on
