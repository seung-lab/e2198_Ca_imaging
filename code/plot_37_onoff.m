%function plot_type_onoff(typefit)

%gc_types = cell_info_typedef_gc();
%gc_types = {gc_types.name};

%typefit = fit_type_means(roi_sums_xcond_typemeans);

overlay_id = 0;
textsize = 13;

xx=[];
yy=[];

[~, tuning_onoff] = tuning_from_fit(coeffs16{3,2});
tuning_onoff_sum = sum(tuning_onoff, 2);
tuning_onoff_max = max(tuning_onoff, [], 2);
tuning_onoff_sum = squeeze(tuning_onoff_sum);
tuning_onoff_max = squeeze(tuning_onoff_max);
%tuning_onoff_sum = tuning_onoff_max;
ca_on_off_ratio = tuning_onoff_sum(1,:)./tuning_onoff_sum(2,:);
%ca_on_off_ratio = tuning_onoff_max(1,:)./tuning_onoff_max(2,:);
ca_on_off_index = (tuning_onoff_sum(1,:) - tuning_onoff_sum(2,:)) ./  (tuning_onoff_sum(1,:) + tuning_onoff_sum(2,:));

figure;
for celltype = {'37'}

	cells = get_ca_cell_info(cell_dict_j, cell_info, celltype);

	yy = [];
	for elem = cells(:).'
		ca_id = cell_dict_j(elem.cell_id==cell_dict_j(:,2), 1);
		%{
		% average of tuning across conditions
		ca_id = get_ca_ids(cell_dict_j, cell_info, elem.cell_id);
		yy(end+1) = ca_on_off_index(ca_id);
		%}

		%%{
		% re-fitting using average trace across conditions
		method = 'exp-exp'; debug_fit=0;
		yactual = roi_sums_xcondmeans(:,ca_id);
		[params, cost, ispoor] = fit_single_onoff(yactual, method, debug_fit);
		%ind_onoff = [4 5]';
		%onoff = params(ind_onoff);
		[~,onoff] = tuning_from_fit(params);
		yy(end+1) = (onoff(1)-onoff(2)) / (onoff(1)+onoff(2));
		%}
	end
	%ca_ids = look_up(cell_dict_j(:,[2 1]), [cells.cell_id]);
	%yy = ca_on_off_index(ca_ids);

	stat1 = cell_info_get_strat_property(cells, 'on');
	stat2 = cell_info_get_strat_property(cells, 'off');
	stat = (stat1-stat2) ./ (stat1+stat2);
	%stat = mean(stat);

	%stat = log10(mean(stat1)/mean(stat2));
	%xx(end+1) = stat;
	%yy(end+1) = log10(onoff(1)/onoff(2))

	xx = stat;

	labels = num2str([cells.cell_id].');

	%{
	exclude = [26162 26103 26158];	% makes r=0.15 p=0.48
	include = ~ismember([cells.cell_id], exclude);
	xx = xx(include);
	yy = yy(include);
	labels = labels(include,:);
	%}
end
length(xx)
h = scatter(xx, yy, 'filled', 'LineWidth', 2); %, 'filled','MarkerFaceColor',kolor,'MarkerEdgeColor',kolor);
hold on
grey = 0.6*[1 1 1];
plot([-1 1], [-1 1], '-.k', 'Color', grey);
%plot([-1 1], [-1 1], '-.k', 'Color', grey, 'LineWidth', 1);

if overlay_id
	text(xx+0.01,yy, labels, 'FontSize', textsize);
	figure_size_x2(2);
else
	figure_size_x2(1.5);
end
ax = gca();
ax.FontSize = 15;

[R, P] = corrcoef(xx, yy)

xtext = 0.3;

%{
xx = xx(:);
yy = yy(:);
bb = [ones(length(xx),1), xx] \ yy;
plot([-1 1], bb(1)+bb(2)*[-1 1], '-.k', 'Color', grey);
text(xtext, bb(1)+bb(2)*xtext+0.01, sprintf('   r = %.2f \n p = %.2f', R(2), P(2)), ...
	 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', textsize);
%}

%{
xx = yy(:);
yy = xx(:);
bb = [ones(length(xx),1), xx] \ yy;
plot(bb(1)+bb(2)*[-1 1], [-1 1], '-.k', 'Color', grey);
%}



xlabel('(Inner-Outer)/(Inner-Outer) in dendritic arbor');
xlabel('$Inner - Outer \over Inner - Outer$ in dendritic arbor', 'Interpreter', 'LaTex');
xlabel('Inner-Outer index')
ylabel('(On-Off)/(On+Off) in calcium');
ylabel('$On-Off \over On+Off$ in calcium', 'Interpreter', 'LaTex');
ylabel('On-Off index')

axis equal;
ax = gca();
xlim([-1 1]);
ylim([-1 1]);
ax.YTick = [-1:0.5:1];

title(sprintf('37, max, r = %.2g  p = %.2g', R(2), P(2)))
title(sprintf('37, mean, r = %.2g  p = %.2g', R(2), P(2)))
title(sprintf('37, fit to mean, fit params, r = %.2g  p = %.2g', R(2), P(2)))
title(sprintf('37, fit to mean, r = %.2g  p = %.2g', R(2), P(2)))
