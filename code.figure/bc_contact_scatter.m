function bc_contact_scatter(cell_contingency, bc_types, celltypes, print_id, raw_vox_count)

filtered = table_sql_where_in(cell_contingency, 'celltype', celltypes);

allcols = filtered.Properties.VariableNames(:).';
D = length(bc_types);
for d = 1:D
    if ~iscell(bc_types{d})
        bc_types{d} = {bc_types{d}};
    end
    tmp = cellfun(@(t) {strncmp(allcols, t, length(t))}, bc_types{d});
    tmp = cell2mat(tmp);
    col = any(tmp, 1);
    xyz(:,d) = sum(filtered{:,col}, 2);
end
if exist('raw_vox_count', 'var') && ~isempty(print_id) && print_id
    xyz(:,:) = xyz(:,:) .* repmat(filtered.total, 1, size(xyz,2));
end

size(xyz, 1)
   % scatter(, cellcontacts.celltype, 30)

allcelltypes = categories(filtered.celltype);
markerAlpha = filtered.total / 1e4;
markerAlpha(markerAlpha>1) = 1;
markerAlpha = markerAlpha / 1.25 + 0.2;

figure;
if exist('print_id', 'var') && ~isempty(print_id) && print_id
    if print_id==1
        labels = strcat(cellstr(filtered.celltype), {' '}, filtered.Properties.RowNames);  % or num2str(filterd.cell)
    elseif print_id==2
        labels = strcat(filtered.Properties.RowNames);
    else 
        labels = strcat(cellstr(num2str(filtered.total)));
    end
else
    labels = char(filtered.celltype);
end
if D<3
    %s = scatter(xyz(:,1),xyz(:,2), [], filtered.celltype, 'filled');
    for ii = 1:size(xyz,1)
        scatter(xyz(ii,1),xyz(ii,2), [], filtered.celltype(ii), 'filled', 'MarkerFaceAlpha', markerAlpha(ii), 'MarkerEdgeColor', 0.9*[1 1 1]);
        hold on
    end
    text(xyz(:,1),xyz(:,2), labels);
else
    %s = scatter3(xyz(:,1),xyz(:,2),xyz(:,3), [], filtered.celltype, 'filled');
    for ii = 1:size(xyz,1)
        scatter3(xyz(ii,1),xyz(ii,2),xyz(ii,3), [], filtered.celltype(ii), 'filled', 'MarkerFaceAlpha', markerAlpha(ii), 'MarkerEdgeColor', 0.9*[1 1 1]);
        hold on
    end
    text(xyz(:,1),xyz(:,2),xyz(:,3), labels);
    grid on
end

xlabel(strjoin(bc_types{1}, ' '));
ylabel(strjoin(bc_types{2}, ' '));
maxi = max(xyz,[],1);
thres = 0.3;
if maxi(1) < thres
    xlim([0 thres]);
end
if maxi(2) < thres
    ylim([0 thres]);
end

if D>2
    zlabel(strjoin(bc_types{3}, ' '));
    if maxi(3) < thres
        zlim([0 thres]);
    end
end
%{
function bc_contact_scatter(contact_table, bctype, celltypes)

% aaaaaaaaaaaaaaaarhhhhhh I want R.
allcelltypes = categories(contact_table.celltype);

D = length(bctype);

        for d = 1:D
            filtered = table_sql_where_in(contact_table, bctype{d})
            sums = filtered;
            sums.Properties.VariableNames{'normalized_sum'} = bctype{d};

            xyz(:,d) = 
            scatter(, cellcontacts.celltype, 30)
            cell_info_get_strat_property(cells, stat_types{d}, false, strat_bc{d});
        end

sums = grpstats(stats, {'cell'}, {'sum'}, 'DataVars',{'normalized_sum'});
sums.Properties.VariableNames{'sum_normalized_sum'} = 'normalized_sum';

cell_stats = join(cell_stats, cells);
cellcontacts.normalized_sum




    if D<3
        scatter(xyz(:,1),xyz(:,2), 'filled','MarkerFaceColor',colors(k,:),'MarkerEdgeColor',colors(k,:));
        for j = 1:size(xyz, 1)
            text(xyz(:,1),xyz(:,2), type_names{k});
        end
    else
        scatter3(xyz(:,1),xyz(:,2),xyz(:,3),'filled','MarkerFaceColor',colors(k,:),'MarkerEdgeColor',colors(k,:));
        for j = 1:size(xyz, 1)
            text(xyz(:,1),xyz(:,2),xyz(:,3), type_names{k});
        end
        grid on
    end
%}
