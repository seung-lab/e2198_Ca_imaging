function fig_ca_summary_brief(roi_sums_means_flatten, tuning, order, ca_ids)

%load('guassiansubstraction.shift=5.mat')
%load('coeffs16.20160819.mat')
load('coeffs16.20160822.mat')


% fit
figure;
n = length(ca_ids);
ncol = 3;

for ii=1:n
	ca_id = ca_ids(ii);

subplot(n, ncol, (ii-1)*ncol+[2 3])
t_frame = 0.128;
t_frame = 4/31;

params = coeffs16{3,2}(ca_id, :);
yfit = fit16{3,2}(:, ca_id);
yfit = yfit.';
y = roi_sums_means_flatten(:, ca_id);
t = 1:31*8;
tt = t*t_frame;
%tt = t;


plot(tt, y);
hold on;
ax = gca();
%plot(tt, yfit, 'LineWidth', 1);
ColorOrderIndex = ax.ColorOrderIndex;
plot(tt(8:end), yfit(8:end), 'LineWidth', 1);
ax.ColorOrderIndex = ColorOrderIndex;
plot([tt(1:8)], [yfit(1:8)], '--', 'LineWidth', 1);
ax.ColorOrderIndex = ColorOrderIndex;
plot([tt+tt(end)], [yfit], '--', 'LineWidth', 1);

plot([-100 100], [params(3) params(3)], '--')
xlim([0 33])
ax.YTick = [-10:10:20];
%ax.XTick = [0:4:32];
ax.XGrid = 'on';
%ax.XMinorTick = 'on';
ax.XTick = sort([1:4:32 2:4:32]);
ax.XTickLabel = [];

figure_size_x2([2 1]);
%applystyle(ax);
%xlabel('time (s)')
%title(ca_id)
ylabel(ca_id)


% tuning
subplot(n, ncol, (ii-1)*ncol+1)
polar_tuning2(tuning(:,:,ca_id),order([7:8 1:6]),2);   % rotated to final coord sys
%applystyle(gca);
%title(ca_id)

end %for
end % func


function applystyle(ax)
	fontsize = 13;

	ax.FontSize = fontsize;
end


