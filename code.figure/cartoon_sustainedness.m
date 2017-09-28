function ax = cartoon_sustainedness(cell_info, cell_dict_j, roi_sums_means_flatten)

style = 3;
hat_style = 4;
slab_style = 2;
with_size_adjustment_slab = 0;
layout = 1;

on_and_off = 1;

xgap = 0.6;
ygap = 0.05;

xdelay_tick = 0;
xdelay_tickline = 0;
xcaliper = 0;
ycaliper = 0;
y_index_tick = 0;
y_index_line = 0;

if on_and_off
	switch style
		case 0
			%nothing
		case 1
			xcaliper = 1;
			ycaliper = 1;
		case 2
			xdelay_tick = 1;
			xdelay_tickline = 1;
			y_index_tick = 1;
			y_index_line = 1;
		case 3	% good one
			y_index_tick = 1;
			y_index_line = 1;
			xcaliper = 1;
		case {4, 5, 6}
			y_index_line = 1;
			ycaliper = 1;
			xdelay_tick = 1;
		otherwise
			body
	end
else
	xcaliper = 1;
	ycaliper = 1;
end


figure;
plot_grouped_ca(cell_info, cell_dict_j, roi_sums_means_flatten, {1}, [],0, 'on',1.7, 0);
legend off
title ''
ax = gca();

hold on
frametime = 0.128;


if on_and_off && 0  % arh affects plot colors...
	yyaxis right
	ylim([0 1.7])

	yyaxis left
end



baseframe = 15 - 1;   % -1 for 0 frame adjustment

if with_size_adjustment_slab  % ON only
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

if on_and_off
	theDot = [theDot; 22, 0.8605];
	peak = 1.582;
end


%plot(theDot(1), theDot(2), 'x', 'MarkerSize', 10)
plot(theDot(:,1), theDot(:,2), '.', 'MarkerSize', 20)



if xcaliper
	if layout==1
		draw_horiz_caliper(1, theDot(1,1), theDot(1,2), 0.55, 0.6, 0.05);
	elseif layout==2
		%draw_horiz_caliper(-theDot(2)/5*1, -theDot(2)/5*2, 0);
		draw_horiz_caliper(1, theDot(1,1), theDot(1,2), -0.15, -0.25, -0.05);
	end

	if on_and_off
		draw_horiz_caliper(2, theDot(2,1), theDot(2,2), -0.4, -0.6, -0.05, -0.3);
	end
end


ax.YTick = [0 1];
% hat bar
hat_styles = {
	'-.', 'Color', [1 1 1]*0.8
	'--', 'Color', [1 1 1]*0.8
	':', 'Color', [1 1 1]*0.3
	'--', 'Color', [1 1 1]*0.5
};
if hat_style
	plot([0 2/frametime], [1 1], hat_styles{hat_style,:})
else
	ax.YGrid = 'on';
end

if on_and_off && hat_style
	plot([2 4]/frametime, peak*[1 1], hat_styles{hat_style,:});
end


if on_and_off
	yyaxis right
	ylim([0 1.7])
	ax.YAxis(2).Color = 'k';

	if ~y_index_tick
		ax.YAxis(2).TickValues = [0 peak];
		ax.YAxis(2).TickLabels = {0, 1};
	else
		ax.YAxis(1).TickValues = [0 theDot(1,2) 1];
		ax.YAxis(2).TickValues = [0 theDot(2,2) peak];
		ontext = 'On';
		%offtext = sprintf('Off sustainedness');
		offtext = 'Off';
		if 1
			ontext = 'On sustainedness';
			offtext = 'Off sustainedness';
		end
		%ax.YTickLabels = {0, offtext, 1};
		%ax.YTickLabels = {'0', offtext, '1'};
		ax.YAxis(2).TickLabels = {'0', offtext, '1'};

		ax.YAxis(1).TickLabels = {'0', ontext, '1'};
	end
	yyaxis left


	if xdelay_tickline
		xx = [theDot(1,1) theDot(1,1)];
		yy = [0 theDot(1,2)-ygap];
		line(xx, yy, 'LineStyle', hat_styles{hat_style,:})
		xx = [theDot(2,1) theDot(2,1)];
		yy = [0 theDot(2,2)-ygap];
		line(xx, yy, 'LineStyle', hat_styles{hat_style,:})
	end

	if y_index_line
		xx = [0 theDot(1,1)-xgap];
		yy = [theDot(1,2) theDot(1,2)];
		line(xx, yy, 'LineStyle', hat_styles{hat_style,:})
		xx = [theDot(2,1)+xgap 4/0.128];
		yy = [theDot(2,2) theDot(2,2)];
		line(xx, yy, 'LineStyle', hat_styles{hat_style,:})
		%plot([xx xx] + dimlineoffset, [0 yy], '-k');
	end

end

if xdelay_tick
	%ax.XTick = sort([[1:4]/frametime theDot(:,1).']);
	if on_and_off
		ax.XTick = sort([1:4 1.8 2.8]/frametime);
	else
		ax.XTick = sort([1:4 1.8]/frametime);
	end
	%ax.XTickLabels = ax.XTick * frametime;	% char array causing labels to be left-aligned...
	ax.XTickLabels = cellstr(num2str(ax.XTick.' * frametime));
else
	ax.XTick = [1 2]/frametime;	 % to be PS'ed
	xlabel('Time')
end

if xdelay_tickline
end


textside = 1;  % right side. 0 = left side
if with_size_adjustment_slab || layout==2
	dimlineoffset = 5;
	extline_far = 6;
	extline_near = 0.6;
elseif xcaliper && layout==1
	dimlineoffset = -1.8;
	extline_far = -3;
	extline_near = -0.6;
elseif style == 4
	dimlineoffset = -5.8;
	extline_far = -7;
	extline_near = -4.6;
elseif style == 5
	dimlineoffset = -2.8;
	extline_far = -4;
	extline_near = -1.6;
	% dimlineoffset = -1.8;
	% extline_far = -3;
	% extline_near = -0.6;
	textside = 0;
elseif style==6
	dimlineoffset = -0;
	extline_far = -0;
	extline_near = -0;
	textside = 0;
else
	dimlineoffset = -1.8;
	extline_far = -3;
	extline_near = -0.6;
end

if ycaliper
	draw_vert_caliper(theDot(1,1), theDot(1,2), dimlineoffset, extline_far, extline_near, ...
		sprintf('Sustainedness \nindex (On)'), textside);

	if on_and_off
		if style==5
			dimlineoffset = 3.8;
			extline_far = 5;
			extline_near = 2.6;
		elseif style==6
			dimlineoffset = -0;
			extline_far = -0;
			extline_near = -0;
			textside = 1;
		end
		draw_vert_caliper(theDot(2,1), theDot(2,2), dimlineoffset, extline_far, extline_near, ...
			sprintf('Sustainedness \nindex (Off)'), textside);
	end
end


ax.FontSize = 13;
if style == 3
	figure_size_x2([1.33 1]);
end


%function fig_draw_dimension_caliper(anchors, xx, yy)
function draw_vert_caliper(xx, yy, dimlineoffset, extline_far, extline_near, txt, textside)

	x1 = xx+dimlineoffset;
	switch textside
	case 0  % text on left
		textoffset = -0.5;
		sideoption = {'HorizontalAlignment', 'right'};
	case 1
		textoffset = 0.5;
		sideoption = {'HorizontalAlignment', 'left'};
	end

	% dim/arrow line
	plot([xx xx] + dimlineoffset, [0 yy], '-k');
	% ext line
	%plot([xx+extline_near xx+extline_far], [yy yy; 0 0].', 'k');
	plot([xx+extline_near xx+extline_far], [yy yy].', '-k');

	text(x1+textoffset, yy/2, txt, 'FontSize', 15, 'VerticalAlignment', 'middle', sideoption{:});
	%plot([xx-0.3 xx+0.1], [0 0], 'k');

	% arrows
	offset = 0.03;
	plot(x1, offset, 'v', x1, yy-offset, '^', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

end

function draw_horiz_caliper(left, xx, yy, dimlineoffset, extline_far, extline_near, extline_near_left)

	if ~exist('extline_near_left', 'var')
		extline_near_left = extline_near;
	end

	y1 = yy+dimlineoffset; % 1.25;
	y2 = yy+extline_far;
	y0 = yy+extline_near;
	%yt = y1 + 0.05; % y text
	
	% ext line
	plot([xx xx], [y0 y2], '-k')
	plot(left*[1 1]/frametime, [yy+extline_near_left y2], '-k')

	% dim/arrow line
	plot([left/frametime xx], [1 1]*y1, '-k');

	% arrows
	offset = 0.3;
	plot(left/frametime+offset, y1, '<', xx-offset, y1, '>', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');


	text(xx-3, y1+0.05, '0.8 s', 'FontSize', 15, 'HorizontalAlignment', 'center');

end

end