function contact_heat_map(contingency_mat, groups1, groups2, labels1, labels2)

if ~exist('labels1', 'var')
	labels1 = {};
	labels2 = {};
end
if ~iscell(labels1)
	labels1 = cellstr(num2str(labels1(:)));
end
if ~iscell(labels2)
	labels2 = cellstr(num2str(labels2(:)));
end


[~,idx] = sort(groups1);
contingency_mat = contingency_mat(idx, :);

[~,idx] = sort(groups2);
contingency_mat = contingency_mat(:, idx);

figure;
ax = imagesc(contingency_mat);
ax = gca();

assert(iscategorical(groups1) && iscategorical(groups2))

bothticks = {};
bothlabels = {};
for group = {groups1 groups2; labels1 labels2}
	group = group{1};

	[counts, categories] = histcounts(group, unique(group));
	% class(categories{1})  % char
	counts
	ticks = cumsum(counts);
	ticks = [ticks - counts / 2; ticks];
	ticks = [0; ticks(:)];
	ticks = ticks + 0.5;  % imagesc's grid cells are centered on the coord values
	ticklabels = cell(size(ticks));
	ticklabels(2:2:end) = cellstr(char(categories));

	bothticks{end+1} = ticks;
	bothlabels{end+1} = ticklabels;
end
ax.XTick = bothticks{2};
ax.YTick = bothticks{1};
ax.XTickLabel = bothlabels{2};
ax.YTickLabel = bothlabels{1};
ax.TickDir = 'out';
%axis equal
set(gca,'dataAspectRatio',[1 1 1])
colorbar()
