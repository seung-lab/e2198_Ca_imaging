function fig_2(cell_info, cell_dict_j, roi_sums_means_flatten)

leaves = ...
{
    '1ni'   3   2   1   1   2   2   -1  -1  -1
    '1no'   3   2   1   1   2   1   -1  -1  -1
    '1ws'   3   2   1   1   1   -1  -1  -1  -1
    '1wt'   3   2   1   2   1   1   -1  -1  -1
    '25'    3   2   1   2   2   1   1   -1  -1
    '27'    3   2   1   2   2   2   1   1   -1
    '28'    3   2   1   2   2   2   1   2   -1
    '2an'   3   2   1   2   2   1   2   -1  -1
    '2aw'   3   2   1   2   2   2   2   2   1
    '2i'    3   2   1   2   2   2   2   1   -1
    '2o'    3   2   1   2   1   2   -1  -1  -1
    '37c'   1   1   -1  -1  -1  -1  -1  -1  -1
    '37d'   1   2   -1  -1  -1  -1  -1  -1  -1
    '37r'   1   3   -1  -1  -1  -1  -1  -1  -1
    '37v'   1   4   -1  -1  -1  -1  -1  -1  -1
    '3i'    3   2   1   2   2   2   1   -1  -1
    '3o'    3   2   1   2   2   2   2   2   2
    '3x'    3   2   1   2   2   2   2   2   3
    '4i'    3   1   4   1   1   -1  -1  -1  -1
    '4on'   3   1   4   1   2   -1  -1  -1  -1
    '4ow'   3   1   4   2   -1  -1  -1  -1  -1
    '51'    3   1   2   2   1   -1  -1  -1  -1
    '58'    3   2   2   1   1   2   2   3   -1
    '5si'   3   1   2   1   1   -1  -1  -1  -1
    '5so'   3   1   2   1   2   -1  -1  -1  -1
    '5ti'   3   1   2   2   2   -1  -1  -1  -1
    '5to'   3   1   3   -1  -1  -1  -1  -1  -1
    '63'    3   1   2   2   3   -1  -1  -1  -1
    '6sn'   3   1   1   2   1   -1  -1  -1  -1
    '6sw'   3   1   1   2   2   -1  -1  -1  -1
    '6t'    3   1   1   1   -1  -1  -1  -1  -1
    '72n'   3   2   2   1   1   1   1   -1  -1
    '72w'   3   2   2   1   1   1   2   -1  -1
    '7id'   2   1   11  -1  -1  -1  -1  -1  -1
    '7ic'   2   1   12  -1  -1  -1  -1  -1  -1
    '7iv'   2   1   13  -1  -1  -1  -1  -1  -1
    '7o'    2   2   -1  -1  -1  -1  -1  -1  -1
    '81i'   3   2   2   1   1   2   1   1   1
    '81o'   3   2   2   1   1   2   2   2   1
    '82n'   3   2   2   1   1   2   1   2   1
    '82wi'  3   2   2   1   1   2   1   2   2
    '82wo'  3   2   2   1   1   2   1   1   2
    '83'    3   2   2   1   1   2   2   2   2
    '8w'    3   2   2   2   -1  -1  -1  -1  -1
    '91n'   3   2   2   1   1   2   2   1   1
    '91w'   3   2   2   1   1   2   2   1   2
    '9m'    3   2   2   1   2   1   -1  -1  -1
    '9n'    3   2   2   1   2   2   -1  -1  -1
};

nullval = -1;
leaves2 = cell2mat(leaves(:, 2:end));
% push division nodes as far right as possible
[nrow ncol] = size(leaves2);
for k = 2:ncol
    for l = 1:nrow
        % shift remaining columns of all siblings (rows with the same parent) towards the right side
        diffparent = leaves2(:, 1:k-1) ~= repmat(leaves2(l, 1:k-1), nrow, 1);
        diffparent = sum(diffparent, 2);
        shift = sum(leaves2(~diffparent, :) ~= nullval, 1);
        shift = find(shift(end:-1:1), 1) - 1;
        if isempty(shift)   % 91n hack for "no group", no shift
            continue;
        end
        leaves2(~diffparent, k:end) = circshift(leaves2(~diffparent, k:end), shift, 2);
    end
end
%leaves2
%return;

base = 10;  % To be safe(depending on the metric/method), need to be bigger than value differences at any level
multiplier = base.^[size(leaves,2)-1:-1:1];
leaves2 = leaves2 .* repmat(multiplier, size(leaves,1), 1);
tree = linkage(leaves2);
%tree = linkage(leaves2, 'single', 'chebychev');    %not really necessary here

figure;
ntypes = nrow;
nfigrows = ntypes;
nfigcols = 3;
figcol = 1;
subplot(nfigrows, nfigcols, nfigcols * [0:nfigrows-1] + figcol);

tree(:, 3) = log(tree(:, 3));
reorder = sum(leaves2, 2);
[~, reorder] = sort(reorder);
reorderedTypes = leaves(reorder, 1);
dendrogram(tree, Inf, 'Orientation', 'left', 'Labels', leaves(:,1), 'Reorder', reorder(end:-1:1));


%{
    '1wt'   2   1   1   1   1   1   1  -1  -1
    '25'    2   1   2   -1  -1  -1  -1  -1  -1
    %'27'    2   2   2   -01 -1  -1  -1  -1  -1
    '27'    2   1+sqrt(.5)   2   -01 -1  -1  -1  -1  -1
%}

%{
strat_groups = {{'1ni' '1no' '1ws'}, {'1wt' '2an' '2aw' '2i' '2o'  '3i' '3o',  '25'}, ...
    {'27' '28' '72w' '81' '82w'}, {'4i' '4on' '4ow' '5ti' '5to' '51'}, {'63', '37', '72n', '83'}, ...
    {'5si' '5so' '6sn' '6sw' '7o' '7i'}, {'58' '82i' '82o' '6t'}, ...
    {'8n' '8w' '9n' '9w' '91w'}, ...
    {'91n'}};
%}

%{
strat_groups = {{'37' '7i' '7o'}  {'6sn' '6sw' '6t'}  {'51' '5si' '5so' '5ti' '63'}, {'4i' '4on' '4ow' '5to' }, ...
    {'1ni' '1no' '1ws'}, {'1wt' '25' '2an' '2aw' '2i' '2o' '3i' '3o'}, {'27' '28'}, ...
    { '72n' '72w' '82wo'}, { '81i' '83'     '91n'}, ...
    {'58' '81o' '82n' '82wi'   }, {'8w' '91w' '9m' '9n'}, ...
    };
%}

%%{
strat_groups = {{'37' '7i' '7o'}  {'6sn' '6sw' '6t'}  {'51' '5si' '5so' '5ti' '63'}, {'4i' '4on' '4ow' '5to' }, ...
    {'1ni' '1no' '1ws' '1wt' '2o' }, {'25' '2an'  '3i'  '27' '28' }, {'2i' '2aw'  '3o' '3x'}, ...
    { '72n' '72w'}, {'81i' '82wo' '82n' '82wi' }, ...
    {     '91n' '91w' '81o'   '83' '58' }, ...
    {  '9m' '9n' '8w'}, ...
    };
%}


figcol = figcol + 1;
startrow = 0;
endrow = 0;
for group = strat_groups
    group = group{1};
    startrow = endrow;
    cells = get_cell_info(cell_info, group);
    nsubtypes = length(unique({cells.type}));
    endrow = startrow + nsubtypes;
    %endrow = startrow + length(group);
    group = reorderedTypes(startrow+1:endrow);  % make sure about order consistency
    subplot(nfigrows, nfigcols, nfigcols * [startrow:endrow-1] + figcol);
    cell_info_plot_strat(cell_info, group, [], 0, 0, 1)
    ax = gca();
    %ax.Color = 'none';
    %ax.Visible = 'off';
    xlim([0 100]);
    xticklabel = ax.XTickLabel;
    ax.XTickLabel = [];
    ax.YTick = [];
    lh = legend(repmat({''}, 1, nsubtypes), 'Location', 'westoutside');
    %legend()
    %legend('boxoff')
    %lh.Position
end
ax.XTickLabel = xticklabel;

figcol = figcol + 1;
%for ctype = reorderedTypes(:).'
    %ctype = ctype{1};
for k = 1:nrow
    ctype = reorderedTypes{k};
    cells = get_cell_info(cell_info, ctype);
    idx = ismember(cell_dict_j(:,2), [cells.cell_id]);
    ca_ids = cell_dict_j(idx,1);
    ca = roi_sums_means_flatten(:,ca_ids);
    ca = mean(ca, 2);
    ca = ca(:, 1);

    n_em = length(cells);
    n_ca = length(ca_ids);
    fprintf('%s \t %d \t %d \n', ctype, n_em, n_ca);


    subplot(nfigrows, nfigcols, nfigcols * [k-1] + figcol);
    plot(ca);
    ax = gca();
    %ax.Color = 'none';
    %ax.Visible = 'off';
    %xlim([0 100]);
    ax.XTick = t1t2_to_ti(9,16); %grid on;
    xticklabel = ax.XTickLabel;
    ax.XTickLabel = [];
    ax.YTick = [];
    ylabel(ctype, 'Rotation', 0, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
end
%ax.XTickLabel = xticklabel;



figure;
ca_groups = {
    'off sustained' {'1ni' '1no' '1ws' '1wt' '25' '27' '2an' '2aw' '2i' '2o' '3i' '3o'   '28' }
    'off transient' {'4i' '4on' '4ow' '5to' }
    'on-off transient' {'51' '5si' '5so' '5ti' '5to' '63'}
    'on transient' {'6sn' '6sw' '6t'}
    %'on sustained' {'28'    '58' '72n' '72w' '81i' '81o' '82n' '82wi' '82wo' '83' '8n' '8w' '91n' '91w' '9m' '9n' '9w'}
    %'on sustained' {    '58' '72n' '72w' '81i' '81o' '82n' '82wi' '82wo' '83' '8n' '8w' '91n' '91w' '9m' '9n' '9w'}
    'on sustained' {    '58' '72n' '72w' '81i' '81o' '82n' '82wi' '82wo' '83'  '8w' '91n' '91w' '9m' '9n' }
};
ca_groups = cell2table(ca_groups, 'VariableNames', {'name', 'types'});

n_groups = 5; %length(ca_groups);
plot_grouped_ca(cell_info, cell_dict_j, roi_sums_means_flatten, ca_groups.types);
legend(ca_groups.name);
title('Calcium');
types = {'1wt'  '4ow' '6sw' '8w'};
figure;plot_grouped_ca(cell_info, cell_dict_j, roi_sums_means_flatten, types)

figure;
for k = 1:n_groups
    %group = strat_groups.types.'
    subplot(n_groups, 1, k)
    plot_grouped_ca(cell_info, cell_dict_j, roi_sums_means_flatten, ca_groups.types{k});
    title(ca_groups.name{k})
end


%cell_info_get_ca 
%roi_sums_means_flatten(:,ind)

