function fig_2(cell_info, cell_dict_j, roi_sums_means_flatten)

leaves = ...
{
    '1ni'   1   2    02 -1  -1  -1  -1  -1
    '1no'   1   2    01 -1  -1  -1  -1  -1
    '1ws'   1   1   -1  -1  -1  -1  -1  -1
    '1wt'   2   1   1   1   -1  -1  -1  -1
    '25'    2   1   2   -1  -1  -1  -1  -1
    '27'    2   2   1    01 -1  -1  -1  -1
    '28'    2   2   1    02 -1  -1  -1  -1
    '2an'   2   1   1   2   1    01 -1  -1
    '2aw'   2   1   1   2   1    02 -1  -1
    '2i'    2   1   1   2   2   2   2    01
    '2o'    2   1   1   2   2   2   1   -1
    '37'    3   3   -1  -1  -1  -1  -1  -1
    '3i'    2   1   1   2   2   1   -1  -1
    '3o'    2   1   1   2   2   2   2    02
    '4i'    3   1   1   2   2   -1  -1  -1
    '4on'   3   1   1   2   1    01 -1  -1
    '4ow'   3   1   1   2   1    02 -1  -1
    '51'    3   1   2    02 -1  -1  -1  -1
    '58'    4   2   1    02 -1  -1  -1  -1
    '5si'   4   1   1    01 -1  -1  -1  -1
    '5so'   4   1   1    02 -1  -1  -1  -1
    '5ti'   3   1   2    01 -1  -1  -1  -1
    '5to'   3   1   1   1   -1  -1  -1  -1
    '63'    3   2   -1  -1  -1  -1  -1  -1
    '6sn'   4   1   2   1    01 -1  -1  -1
    '6sw'   4   1   2   1    02 -1  -1  -1
    '6t'    4   2   1    01 -1  -1  -1  -1
    '72n'   3   4   -1  -1  -1  -1  -1  -1
    '72w'   2   2   2    01 -1  -1  -1  -1
    '7i'    4   1   3   -1  -1  -1  -1  -1
    '7o'    4   1   2   2   -1  -1  -1  -1
    '81i'    2   2   2    03 -1  -1  -1  -1
    '82wi'   4   2   2    01 -1  -1  -1  -1
    '81o'   4   2   2    02 -1  -1  -1  -1
    '82wo'   2   2   2    02 -1  -1  -1  -1
    '83'    3   5   -1  -1  -1  -1  -1  -1
    '8n'    5   2    02 -1  -1  -1  -1  -1
    '8w'    5   2    03 -1  -1  -1  -1  -1
    '91n'   10  -1  -1  -1  -1  -1  -1  -1
    '91w'   5   1   -1  -1  -1  -1  -1  -1
    '9n'    5   2    01 -1  -1  -1  -1  -1
    '9w'    5   2    04 -1  -1  -1  -1  -1
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
subplot(nfigrows, nfigcols, nfigcols * [0:nfigrows-1] + 1);

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

strat_groups = {{'1ni' '1no' '1ws'}, {'1wt' '2an' '2aw' '2i' '2o'  '3i' '3o',  '25'}, ...
    {'27' '28' '72w' '81' '82w'}, {'4i' '4on' '4ow' '5ti' '5to' '51'}, {'63', '37', '72n', '83'}, ...
    {'5si' '5so' '6sn' '6sw' '7o' '7i'}, {'58' '82i' '82o' '6t'}, ...
    {'8n' '8w' '9n' '9w' '91w'}, ...
    {'91n'}};

startrow = 0;
endrow = 0;
for group = strat_groups
    group = group{1};
    startrow = endrow;
    endrow = startrow + length(group);
    group = reorderedTypes(startrow+1:endrow);  % make sure about order consistency
    subplot(nfigrows, nfigcols, nfigcols * [startrow:endrow-1] + 2);
    cell_info_plot_strat(cell_info, group, [], 0, 0, 1)
    ax = gca();
    %ax.Color = 'none';
    %ax.Visible = 'off';
    xlim([0 100]);
    xticklabel = ax.XTickLabel;
    ax.XTickLabel = [];
    ax.YTick = [];
    lh = legend(repmat({''}, 1, length(group)), 'Location', 'westoutside');
    %legend()
    %legend('boxoff')
    %lh.Position
end
ax.XTickLabel = xticklabel;

%for ctype = reorderedTypes(:).'
    %ctype = ctype{1};
for k = 1:nrow
    ctype = reorderedTypes(k);
    cells = get_cell_info(cell_info, ctype);
    idx = ismember(cell_dict_j(:,2), [cells.cell_id]);
    ca_ids = cell_dict_j(idx,1);
    ca = roi_sums_means_flatten(:,ca_ids);
    ca = mean(ca, 2);
    ca = ca(:, 1);


    subplot(nfigrows, nfigcols, nfigcols * [k-1] + 3);
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

%cell_info_get_ca 
%roi_sums_means_flatten(:,ind)

