function cartoon_sustainedness(cell_info, cell_dict_j, roi_sums_means_flatten)

style = 1;
hat_style = 2;
slab_style = 2;
with_size_adjustment_slab = 0;
mark_measurement_time = 1;
layout = 1;

figure;
plot_grouped_ca(cell_info, cell_dict_j, roi_sums_means_flatten, {1}, [],0, 'on',1.7, 0);
legend off
title ''
ax = gca();

hold on
frametime = 0.128;


baseframe = 15 - 1;   % -1 for 0 frame adjustment

if with_size_adjustment_slab
	switch slab_style
	case 1
		%plot([1 1]/frametime, [0 0.8], '--k')
		plot([baseframe baseframe], [0 0.8], '--k')
		plot([baseframe baseframe]-1.6, [0 0.8], '--k')
	case 2

		height = 1.3;
		plot([1 1]/frametime, [0 1.3], '--k')
		plot([baseframe baseframe], [0 1.3], '--k')
		plot([baseframe baseframe]-1.6, [0 0.8], '--k')
		%draw_horiz_caliper();
	end
	theDot = [13, 0.724];
else
	theDot = [14, 0.631];
end


%plot(theDot(1), theDot(2), 'x', 'MarkerSize', 10)
plot(theDot(1), theDot(2), '.', 'MarkerSize', 20)


if mark_measurement_time
	if layout==1
		draw_horiz_caliper(0.55, 0.6, 0.05);
	elseif layout==2
		%draw_horiz_caliper(-theDot(2)/5*1, -theDot(2)/5*2, 0);
		draw_horiz_caliper(-0.15, -0.25, -0.05);
	end
end


ax.YTick = [0 1];
% hat bar
if hat_style == 1
	%plot([0 2/frametime], [1 1], '--k')
	plot([0 2/frametime], [1 1], '-.', 'Color', [1 1 1]*0.8)
elseif hat_style == 2
	plot([0 2/frametime], [1 1], '--', 'Color', [1 1 1]*0.8)
elseif hat_style == 3
	plot([0 2/frametime], [1 1], ':', 'Color', [1 1 1]*0.3)
else
ax.YGrid = 'on';
end


if with_size_adjustment_slab || layout==2
	dimlineoffset = 5;
	extline_far = 6;
	extline_near = 0.6;
elseif mark_measurement_time && layout==1
	dimlineoffset = -1.8;
	extline_far = -3;
	extline_near = -0.6;
else
	dimlineoffset = -1.8;
	extline_far = -3;
	extline_near = -0.6;
end
%textloc = 
draw_vert_caliper(dimlineoffset, extline_far, extline_near);

ax = gca();
ax.FontSize = 13;

%function fig_draw_dimension_caliper(anchors, xx, yy)
function draw_vert_caliper(dimlineoffset, extline_far, extline_near)

	xx = theDot(1);
	yy = theDot(2);

	x1 = xx+dimlineoffset;

	% dim/arrow line
	plot([xx xx] + dimlineoffset, [0 yy], 'k');
	% ext line
	%plot([xx+extline_near xx+extline_far], [yy yy; 0 0].', 'k');
	plot([xx+extline_near xx+extline_far], [yy yy].', 'k');

	text(x1+0.5, yy/2, sprintf('Sustainedness \nindex (On)'), 'FontSize', 15, 'VerticalAlignment', 'middle');
	%plot([xx-0.3 xx+0.1], [0 0], 'k');

	% arrows
	offset = 0.03;
	plot(x1, offset, 'v', x1, yy-offset, '^', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

end

function draw_horiz_caliper(dimlineoffset, extline_far, extline_near)
	xx = theDot(1);
	yy = theDot(2);

	y1 = theDot(2)+dimlineoffset; % 1.25;
	y2 = theDot(2)+extline_far;
	y0 = theDot(2)+extline_near;
	%yt = y1 + 0.05; % y text
	
	% ext line
	plot([xx xx], [y0 y2], 'k')
	plot([1 1]/frametime, [y0 y2], 'k')

	% dim/arrow line
	plot([1/frametime xx], [1 1]*y1, 'k');

	% arrows
	offset = 0.3;
	plot(1/frametime+offset, y1, '<', xx-offset, y1, '>', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');


	text(xx-3, y1+0.05, '0.8 s', 'FontSize', 15, 'HorizontalAlignment', 'center');

end

end