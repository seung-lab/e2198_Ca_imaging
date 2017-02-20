function [contact_table, stats, cell_stats, cell_contingency] = gc_bc_contacts(cell_info)

contact_file_path = 'gc_bc_contacts.mat';
contact_file_path = 'gc_bc_contacts.20161028.mat';

load(contact_file_path);    % confusion_mat

tic
bc_ids = unique(confusion_mat(1,:));
cell_ids = unique(confusion_mat(2,:));

cells = get_cell_info(cell_info, bc_ids);
cells = struct2table(cells);
cells = cells(:,{'cell_id', 'type'});
cells.type(cellfun(@isempty,cells.type)) = {''};
cells.type = categorical(cells.type);
bc = cells;
bc.Properties.VariableNames{'cell_id'} = 'bc';
bc.Properties.VariableNames{'type'} = 'bctype';

cells = get_cell_info(cell_info, cell_ids);
cells = struct2table(cells);
cells = cells(:,{'cell_id', 'type'});
cells.type(cellfun(@isempty,cells.type)) = {''};
cells.type = categorical(cells.type);
cells.Properties.VariableNames{'cell_id'} = 'cell';
cells.Properties.VariableNames{'type'} = 'celltype';

contact_table = array2table(confusion_mat.', 'VariableNames', {'bc', 'cell', 'count'});
contact_table = join(contact_table, bc);
contact_table = innerjoin(contact_table, cells); %some cells are not in our cell list, 89xxx
toc

%{
%contact_table = cell2table(num2cell(confusion_mat.'));
contact_table = num2cell(confusion_mat.');
contact_table(:,4:5) = {''};
contact_table(1:3,4:5)
size(contact_table)
%contact_table(:,4:5) = repmat({''}, size(contact_table, 1), 2);
%contact_table(1:3,4:5)
%contact_table{1,4}={'e'};
%contact_table(1:3,4:5)

tic
for jj = 1:2
	for ii = 1:size(contact_table, 1)
		id = contact_table{ii, jj};
		c = get_cell_info(cell_info, id);
		if isempty(c)
			if id < 89000 || id > 90000
				warning(sprintf('non exist: %d', id))
			end
			%contact_table{ii, 3+jj} = '';
			continue
		end
		ctype = c.type;
		if isempty(ctype)
			ctype = '';
		end
		contact_table{ii, 3+jj} = ctype; %cell array
		%contact_table{ii, 3+jj} = {ctype}; %table
	end
end
toc
contact_table = cell2table(contact_table, 'VariableNames', {'bc', 'cell', 'count', 'bctype', 'celltype'});
contact_table.bctype = categorical(contact_table.bctype);
contact_table.celltype = categorical(contact_table.celltype);
%}

% by type
stats = grpstats(contact_table, {'celltype','bctype'}, {'sum', 'mean', 'min', 'max'}, 'DataVars',{'count'});

sums = grpstats(stats, {'celltype'}, {'sum'}, 'DataVars',{'sum_count'});
normalized_sum = rowfun(@(sum_count, celltype) sum_count / sums.sum_sum_count(char(celltype)), stats(:, {'sum_count', 'celltype'}));
stats.normalized_sum = normalized_sum{:,:};

% by cell
cell_stats = grpstats(contact_table, {'cell','bctype'}, {'sum'}, 'DataVars',{'count'});
sums_cell = grpstats(cell_stats, {'cell'}, {'sum'}, 'DataVars',{'sum_count'});
normalized_sum = rowfun(@(sum_count, cell) sum_count / sums_cell.sum_sum_count(num2str(cell)), cell_stats(:, {'sum_count', 'cell'}));
cell_stats.normalized_sum = normalized_sum{:,:};
cell_stats = join(cell_stats, cells);


allbctypes = categories(cell_stats.bctype);
allbctypes = cellfun(@(s) {s(s~='/')}, allbctypes);  % bc8/9 => bc89
cell_contingency = sums_cell(:,[]);
cell_contingency{:, allbctypes} = zeros(size(sums_cell,1),length(allbctypes));
tic
for row = 1:size(cell_stats,1)
	row = cell_stats(row, :);
	cell_contingency{num2str(row.cell), int64(row.bctype)} = row.normalized_sum; %row.sum_count;
end
toc
%cell_contingency.total = sum(cell_contingency{:,:}, 2);
%cell_contingency{:, 1:end-1} = 
cell_contingency = join(cell_contingency, sums_cell(:, {'sum_sum_count', 'cell'}), 'Keys','RowNames');
cell_contingency.Properties.VariableNames{'sum_sum_count'} = 'total';
%cell_contingency = join(cell_contingency, cells, 'LeftKeys', 'RowNames', 'RightKeys', 'cell');
cell_contingency = join(cell_contingency, cells);

plot_bc_contacts(stats);
%{
figure
scatter(stats.bctype, stats.celltype, stats.normalized_sum*100, stats.normalized_sum, 'filled')
ax = gca;
ax.XTick = 1:length(categories(stats.bctype));
ax.XTickLabels = categories(stats.bctype);
ax.YTick = 1:length(categories(stats.celltype));
ax.YTickLabels = categories(stats.celltype);
ylim([0, 1+length(categories(stats.celltype))])
ax.XTickLabelRotation = 30;
ax.GridAlpha = 0.05;
grid on
%}


%{
hold on
%for cate = categories(stats.celltype)
for cate = unique(stats.celltype).'
	rowss = stats(stats.celltype == cate, :);
	plot(rowss.bctype, rowss.normalized_sum);
end
%}
