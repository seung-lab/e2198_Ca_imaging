function fig_ca_typefit(roi_sums_xcond_typemeans, typefit, order, typename)


if ~exist('typename', 'var')

typename = '37v';
typename = '37d';
end

%load('coeffs16.20160822.mat');


t_frame = 0.128;
t_frame = 4/31;

params = typefit{1,3}{typename, 'params'};
y = roi_sums_xcond_typemeans{typename,:};
t = 1:31;
tt = t*t_frame;

method = 'exp-exp';
ff = @(x, xdata) (ftrial(x(1:end-2), method, x(end-1:end), 1));
[yfit, CAg] = ff(params);
C = CAg(1);
%CAg(2:end,:) = CAg(2:end,:) + CAg(1);
%CAg(2:end,1:31) = CAg(2:end,1:31) + CAg(2:end,32:end);
n = (size(CAg,1)-1)/2;
for k = 3:2*n+1
	CAg(k,:) = CAg(k,:) + CAg(k-1,:);
end
%CAg(2:end,1:31) = CAg(2:end,1:31) + repmat(CAg(end,32:end), n, 1);
CAg(1,:) = 0;
CAg(:,1:31*n) = CAg(:,1:31*n) + repmat(CAg(end,1+31*n:end), 2*n+1, 1);
CAg = CAg + C;


ti = params(end-1:end);
%{
[~, g_raw] = generate_alphas(params(1:2), n, ti, method);
%parts = 
g_raw
CAg(2:end,:) = CAg(2:end,:) + CAg(1);
%}

% fit
figure;
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


tt2 = [1:size(CAg,2)]*t_frame;
%plot(tt2, CAg.', '--', 'LineWidth', 1);
for k = 1:2*n
	kolor = ax.ColorOrder(k,:);
	kolor = kolor*0.8;
	fill([tt2 flip(tt2)], [CAg(k, :) flip(CAg(k+1, :))], kolor, 'FaceAlpha', 0.5, 'LineStyle', 'none');%, 'EdgeAlpha', 0);
	fill([tt flip(tt)], [CAg(k, 31*n+1:end) flip(CAg(k+1, 31*n+1:end))], kolor, 'FaceAlpha', 0.5, 'LineStyle', 'none');%, 'EdgeAlpha', 0);
end

%fill([t_e*t_frame flip(t_e*t_frame)], [e1, flip(e2)], kolor, 'FaceAlpha', 0.5, 'LineStyle', 'none');%, 'EdgeAlpha', 0);


plot([-100 100], [params(3) params(3)], '--')
xlim([0 5])
ax.YTick = [-10:10:20];
ax.XTick = [0:4:32];
ax.XGrid = 'on';
ax.XMinorTick = 'on';
ax.XAxis.MinorTickValues = sort([ti*t_frame 0 4]);
ax.XMinorGrid = 'on';

figure_size_x2([2 1]);
xlabel('time (s)')
title(typename)



end %func


function g = exp_exp(t, ti, tau)
		t_ti = t - ti;  t_ti(t_ti<0) = 0;
		g = exp(-t_ti/tau(1)) - exp(-t_ti/tau(2));
		g(t_ti<=0) = 0;
end
