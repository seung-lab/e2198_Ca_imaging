function cluster_linkage(cell_info, types, method)

suson = { '2an' '2aw' '2i' '2o' '3i' '3o'};
suson = {'1wt' '2an' '2aw' '2i' '2o' '3i' '3o'};
suson = {'1wt' '2aw' '2i' '2o' '3i' '3o'};
suson = {  '2aw' };
suson = {'1wt' '2i' '2o' '3i' '3o'};
suson = {  '2aw' '2i' '2o' '3i' '3o'};
%types = suson;

suson = get_cell_info(cell_info, types);
susonstrat = cat(2, suson.strat_nrml);
susonstrat = susonstrat(:, 2:2:end);

%{
metric = 'naive';
tree = linkage(susonstrat.');
%}

%%{

if ~exist('method', 'var')
method = 'median';
method = 'average';
method = 'centroid';
method = 'ward';
method = 'complete';
%method = 'weighted';
%method = 'single';
end

metric = 'spearman';
metric = 'cosine';
metric = 'correlation';

tree = linkage(susonstrat.', method, metric);
%}

%labels = strcat({suson.type}.', {sprintf('\t_')}, cellstr(num2str([suson.cell_id].')));
labels = strcat(cellstr(num2str([suson.cell_id].')), {'  '}, {suson.type}.');
%labels = {suson.type}.';
figure; dendrogram(tree, Inf, 'Orientation', 'left', 'Labels', labels)
%{ no way to display a tab '\t' in figures...
ax = gca();
ax.TickLabelInterpreter = 'none';
%}

title([strjoin(types)  '    '  method  '    '  metric])
