%function plot_type_onoff(typefit)

%gc_types = cell_info_typedef_gc();
%gc_types = {gc_types.name};

%typefit = fit_type_means(roi_sums_xcond_typemeans);

overlay_id = 1;
textsize = 13;

ind_onoff = [4 5]';
xx=[];
yy=[];

figure;
for celltype = roi_sums_xcond_typemeans.Properties.RowNames(:).'  % gc_types(:).'
	%celltype

	%{
	if strcmp(celltype, '27') || strcmp(celltype, '1ws')
		continue;
	end
	% Using coeffs:
	% R = 0.8421
	% P = 7.8493e-13
	% with 27 and 1ws: R = 0.7230   P = 1.3902e-08
	% Using peaks:
	% R = 0.84
	% P = 8.6e-13
	% with 27 and 1ws: R = 0.7222   P = 1.46e-08
	%}

	cells = get_ca_cell_info(cell_dict_j, cell_info, celltype);
	cells = get_cell_info(cell_info, celltype);

	params = typefit{3}{celltype, 'params'};
	onoff = params(ind_onoff);
	[~,onoff] = tuning_from_fit(params);

	stat1 = cell_info_get_strat_property(cells, 'on');
	stat2 = cell_info_get_strat_property(cells, 'off');
	stat = (stat1-stat2) ./ (stat1+stat2);
	stat = mean(stat);

	%stat = log10(mean(stat1)/mean(stat2));
	%xx(end+1) = stat;
	%yy(end+1) = log10(onoff(1)/onoff(2))

	xx(end+1) = stat;
	yy(end+1) = (onoff(1)-onoff(2)) / (onoff(1)+onoff(2));
end
h = scatter(xx, yy, 'LineWidth', 2); %, 'filled','MarkerFaceColor',kolor,'MarkerEdgeColor',kolor);
hold on
grey = 0.6*[1 1 1];
plot([-1 1], [-1 1], '-.k', 'Color', grey);
if overlay_id
	text(xx+0.02,yy, roi_sums_xcond_typemeans.Properties.RowNames(:), 'FontSize', textsize);
end
figure_size_x2(1.5);
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


%{
xlabel('(Inner-Outer)/(Inner-Outer) in dendritic arbor');
xlabel('$Inner - Outer \over Inner - Outer$ in dendritic arbor', 'Interpreter', 'LaTex');
ylabel('(On-Off)/(On+Off) in calcium');
ylabel('$On-Off \over On+Off$ in calcium', 'Interpreter', 'LaTex');
%}
xlabel('Inner-Outer index')
ylabel('On-Off index')

axis equal;
ax = gca();
xlim([-1 1]);
ylim([-1 1]);
ax.YTick = [-1:0.5:1];

title(sprintf('r = %.2g  p = %.2g', R(2), P(2)))
