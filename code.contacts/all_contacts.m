function [contact_table, stats, cell_stats] = all_contacts(cell_info)

contact_file_path = 'all_contacts.mat';
contact_file_path = 'all_contacts.600-1167.20161219.mat';
surfacearea_file_path = 'contacts/raw3d_445-1167/surface_area.mat';  %mismatch

load(contact_file_path);    % confusion_mat

tic

confusion_mat = [confusion_mat confusion_mat([2 1 3], :)];  % append a copy with swapped cell ids

cell1_ids = unique(confusion_mat(1,:));
cell2_ids = unique(confusion_mat(2,:));

cells = get_cell_info(cell_info, cell1_ids);
cells = struct2table(cells);
cells = cells(:,{'cell_id', 'type'});
cells.type(cellfun(@isempty,cells.type)) = {''};
cells.type = categorical(cells.type);
bc = cells;
bc.Properties.VariableNames{'cell_id'} = 'cell1';
bc.Properties.VariableNames{'type'} = 'cell1_type';

cells = get_cell_info(cell_info, cell2_ids, [], false);
cells = struct2table(cells);
cells = cells(:,{'cell_id', 'type'});
cells.type(cellfun(@isempty,cells.type)) = {''};
cells.type = categorical(cells.type);
cells.Properties.VariableNames{'cell_id'} = 'cell2';
cells.Properties.VariableNames{'type'} = 'cell2_type';
length(cells.cell2)
length(unique(cells.cell2))

contact_table = array2table(confusion_mat.', 'VariableNames', {'cell1', 'cell2', 'count'});
contact_table = innerjoin(contact_table, bc);
contact_table = innerjoin(contact_table, cells); % innerjoin: some cells are not in our cell list, 89xxx
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

tic
% remove SACs which are the major contacting sources
if 0	% arh apparently I didn't remove other non-SAC ACs...
	contacts = table_sql_where_in(contact_table, 'cell1_type', {'ON SAC', 'OFF SAC'}, 'not');
else
	contacts = contact_table;
end
toc

tic
% by type
stats = grpstats(contacts, {'cell1_type','cell2_type'}, {'sum', 'mean', 'min', 'max'}, 'DataVars',{'count'});

sums = grpstats(stats, {'cell2_type'}, {'sum'}, 'DataVars',{'sum_count'});
normalized_sum = rowfun(@(sum_count, celltype) sum_count / sums.sum_sum_count(char(celltype)), stats(:, {'sum_count', 'cell2_type'}));
stats.normalized_sum = normalized_sum{:,:};
toc

tic
% by cell
cell_stats = grpstats(contacts, {'cell2','cell1_type'}, {'sum'}, 'DataVars',{'count'});
sums_cell = grpstats(cell_stats, {'cell2'}, {'sum'}, 'DataVars',{'sum_count'});
normalized_sum = rowfun(@(sum_count, cell) sum_count / sums_cell.sum_sum_count(num2str(cell)), cell_stats(:, {'sum_count', 'cell2'}));
cell_stats.normalized_sum = normalized_sum{:,:};
cell_stats = join(cell_stats, cells);
% add sum_sum_count too for investigation purpose
cell_stats = join(cell_stats, sums_cell(:,[1 3]));  % sums_cell: cell2    GroupCount    sum_sum_count
toc

%{
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
%}

plot_all_contacts(stats);
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
