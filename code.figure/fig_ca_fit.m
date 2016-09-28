function fig_ca_fit(roi_sums_means_flatten, tuning, order, ca_id, coeffs16)

if ~exist('ca_id', 'var')
% 90001 / 508
ca_id = 508;
end

%generate_alphas
%load('coeffs16.20160505.mat');

%if ~exist('coeffs16', 'var')
%end
%load('guassiansubstraction.shift=5.mat')

load('coeffs16.20160822.mat');


t_frame = 0.128;
t_frame = 4/31;

params = coeffs16{3,2}(ca_id, :);
yfit = fit16{3,2}(:, ca_id);
yfit = yfit.';
y = roi_sums_means_flatten(:, ca_id);
t = 1:31*8;
tt = t*t_frame;
%offset = coeffs16{3,2}(ca_id, 1:2);

plot_fit=1;
if plot_fit

% exp-exp function
figure;
tau = params(1:2);
%{
mag = tuning_from_fit(params);
assert(isequal(mag, tuning(:,:,ca_id)))
mag ./ (params(end-32+1:2:end-16) + params(end-32+2:2:end-16))
%}
step_size = 1;
step_size = 0.2;
t_ti = [-20 -1e-8 0:step_size:50];
g = exp(-t_ti/tau(1)) - exp(-t_ti/tau(2));
g(t_ti<0) = 0;
[maxi, ind] = max(g)

t_e = [0:step_size:50];
e1 = exp(-t_e/tau(1)); e1(t_e<0) = 1;
e2 = exp(-t_e/tau(2)); e2(t_e<0) = 1;
e2 = 1-e2;

plot(t_ti*t_frame, g, 'LineWidth', 1.5);
hold on
plot(t_e*t_frame, [e1; e2], '--', 'LineWidth', 1);

xx = (t_ti(ind)+5)*t_frame;
plot([xx xx], [0 maxi], 'k');
plot([xx-0.3 xx+0.15], [maxi maxi; 0 0].', 'k');
text(xx+0.05, maxi-0.1, 'R', 'FontSize', 15);
%plot([xx-0.3 xx+0.1], [0 0], 'k');
offset = 0.02;
plot(xx, offset, 'v', xx, maxi-offset, '^', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');


xlim([-1, 3])
ylim([-0.3, 1.2])
ax = gca();
ax.YTick = [0 1];
ax.XTick = [0];
applystyle(ax);
xlabel('t', 'FontAngle', 'italic');
grid on
ax.GridLineStyle = '--';
legend({'A\cdotg(t)' 'A\cdote^{-t/\tau_1}' 'A\cdot(1-e^{-t/\tau_2})'})
legend({'A\cdotg(t)' 'A\cdot(1-exp(-t/\tau_1))' 'A\cdot(1-exp(-t/\tau_2))'})
legend({'A\cdotg(t)' 'A\cdot(1-exp(-t/tau1))' 'A\cdot(1-exp(-t/tau2))'})
%legend({'A g(t)' 'A exp(-t/tau1)' 'A (1-exp(-t/tau2))'})
legend({'A*g(t)' 'A*exp(-t/tau1)' 'A*(1-exp(-t/tau2))'})
legend({'$A\cdot g(t)$,  $g(t) = e^{-t/\tau_1}-e^{-t/\tau_2}$' '$A\cdot e^{-t/\tau_1}$' '$A\cdot (1-e^{-t/\tau_2})$'}, 'Interpreter', 'LaTex');
legend({'$A\cdot (e^{-t/\tau_1}-e^{-t/\tau_2})$' '$A\cdot e^{-t/\tau_1}$' '$A\cdot (1-e^{-t/\tau_2})$'}, 'Interpreter', 'LaTex');
legend({'$A\cdot g(t)$' '$A\cdot e^{-t/\tau_1}$' '$A\cdot (1-e^{-t/\tau_2})$' '$g(t) = e^{-t/\tau_1}-e^{-t/\tau_2}$'}, 'Interpreter', 'LaTex');
ax.YTickLabel{2} = 'A';

%% exp-exp function v2
figure;
legends = {'g(t)' 'exp(-t/tau1)' 'exp(-t/tau2)'};
legends = {'g(t)' 'exp(-t/\tau_1)' 'e^{-t/\tau_2}'};
legends = {'$g(t)$' 'exp$(-t/\tau_1)$' '$e^{-t/\tau_2}$'};
legends = {'$g(t)$' '$e^{-t/\tau_1}$' '$e^{-t/\tau_2}$'};

t_e = [0:step_size:50];
e1 = exp(-t_e/tau(1)); e1(t_e<0) = 1;
e2 = exp(-t_e/tau(2)); e2(t_e<0) = 1;

subplot(2,1,2)
plot(t_ti*t_frame, g, 'LineWidth', 1.5);
legend({'$A \cdot g(t)$'}, 'Interpreter', 'LaTex');

subplot(2,1,1)
ax = gca;
kolor = ax.ColorOrder(1,:);
kolor = kolor*0.8;
fill([t_e*t_frame flip(t_e*t_frame)], [e1, flip(e2)], kolor, 'FaceAlpha', 0.5, 'LineStyle', 'none');%, 'EdgeAlpha', 0);
fill([t_e*t_frame flip(t_e*t_frame)], [e1, flip(e2)], kolor+(1-kolor)/2,  'LineStyle', 'none');%, 'EdgeAlpha', 0);
hold on;
plot(t_e*t_frame, [e1; e2], '--', 'LineWidth', 1.5);
%legend(legends, 'FontAngle', 'italic');
legend(legends, 'Interpreter', 'LaTex');

subplot(2,1,2)
hold on
xx = (t_ti(ind)+5)*t_frame;
plot([xx xx], [0 maxi], 'k');
plot([xx-0.3 xx+0.15], [maxi maxi; 0 0].', 'k');
text(xx+0.05, maxi-0.15, 'R', 'FontSize', 15);
%plot([xx-0.3 xx+0.1], [0 0], 'k');
offset = 0.06;
plot(xx, offset, 'v', xx, maxi-offset, '^', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');


for subp = 1:2
subplot(2,1,subp)
xlim([-1, 3])
ylim([-0.3, 1.2])
ax = gca();
ax.YTick = [0 1];
ax.XTick = [0];
ax.GridLineStyle = '--';
applystyle(ax);
grid on
end
subplot(2,1,2)
xlabel('t', 'FontAngle', 'italic');
ax.YTickLabel{2} = 'A';

end % if plot_fit

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

plot([-100 100], [params(3) params(3)], '--')
xlim([0 33])
ax.YTick = [-10:10:20];
ax.XTick = [0:4:32];
ax.XGrid = 'on';
ax.XMinorTick = 'on';
ax.XAxis.MinorTickValues = sort([1:4:32 2:4:32]);

figure_size_x2([2 1]);
applystyle(ax);
xlabel('time (s)')
title(ca_id)

% fit v2
figure;
%lines = [];
clear lines;
method = 'exp-exp';
ff = @(x, xdata) (ftrial(x(1:end-16), method, x(end-16+1:end)));
%ff = @(x, xdata) (ftrial(x(1:end-2), method, x(end-1:end), 1));
[yfit, CAg] = ff(params);
ti = params(end-16+1:end);
C = CAg(1);
n = (size(CAg,1)-1)/2;
for k = 3:2*n+1
	CAg(k,:) = CAg(k,:) + CAg(k-1,:);
end
CAg(1,:) = 0;
CAg(:,1:31*n) = CAg(:,1:31*n) + repmat(CAg(end,1+31*n:end), 2*n+1, 1);
CAg = CAg + C;
hold on
tt2 = [1:size(CAg,2)]*t_frame;
%plot(tt2, CAg.', '--', 'LineWidth', 1);
ax.ColorOrder = [ax.ColorOrder; 0 0 0]; % make it 8

%plot(tt, y, 'k', 'LineWidth', 1);
for k = 1:2*n
	kolor = ax.ColorOrder(1+mod(k,2),:);
	%kolor = kolor;
	kolor = {'r', 'k'};
	kolor = {[0.3 0.3 0.3], 'k'};
	kolor = kolor{1+mod(k-1,2)};
	kolor = ax.ColorOrder(1+mod(k,8),:);
	alphav = 0.7;
	tx = floor(ti(k));
	len = 30;
	%plot(tt2(tx:tx+len), CAg(1+k, tx:tx+len), 'LineWidth', 2);
	%%{
	fill([tt2 flip(tt2)], [CAg(k, :) flip(CAg(k+1, :))], kolor, 'FaceAlpha', alphav, 'LineStyle', 'none');%, 'EdgeAlpha', 0);
	%CAg(k, 31*n+1:end)
	fill([tt flip(tt)], [CAg(k, 31*n+1:end) flip(CAg(k+1, 31*n+1:end))], kolor, 'FaceAlpha', alphav, 'LineStyle', 'none');%, 'EdgeAlpha', 0);
	%}
end
lines(1) = plot(tt, y, 'k', 'LineWidth', 2);
%plot(tt, y, '.');
%plot(tt, yfit, '-', 'LineWidth', 1);

kolor = ax.ColorOrder(3,:);
t1 =  floor(ti(1));
%%{
lines(2) = plot(tt(t1:end), yfit(t1:end), 'LineWidth', 2, 'Color', kolor);
plot([tt(1:t1)], [yfit(1:t1)], '--', 'LineWidth', 2, 'Color', kolor);
plot([tt(end) tt+tt(end)], [yfit(end) yfit.'], '--', 'LineWidth', 2, 'Color', kolor);
%}

lines(3) = plot([-100 100], [C C], '--k');
text
text1toN = num2str([1:2*n].');
texts = [repmat('A_{', 2*n, 1) text1toN repmat('}*g_{', 2*n, 1) text1toN repmat('}', 2*n, 1)];
%text(ti*t_frame, y(ceil(ti)+1)+0.8, texts);

%kolor = ax.ColorOrder(1+mod(k,2),:);
xlim([0 32])
ax = gca;
ax.XTick = sort([ti*t_frame]);

ax.XTickLabels = [repmat('$t_{', 2*n, 1) num2str([1:2*n].') repmat('}$', 2*n, 1)]; % {'t_1' 't_2'}
ax.XGrid = 'on';
ax.GridLineStyle = ':';

ax.XAxis.MinorTickValues = sort([0:4:32]);
ax.XMinorTick = 'on';
ax.XMinorGrid = 'on';
ax.MinorGridLineStyle = '-';

ax.GridColor = 0.1*[1 1 1];
ax.MinorGridColor = 0.15*[1 1 1];
ax.GridAlpha = 0.25;
ax.MinorGridAlpha = 0.15;

[ax.YTick, idx] = sort([C -10:10:20]);
ax.YTickLabels{idx(1)} = 'C';
%ax.YAxis.MinorTickValues = C;
%ax.YMinorTick = 'on';

figure_size_x2([2 1]);
applystyle(ax);
xlabel('t', 'FontAngle', 'italic');
ax.XAxis.TickLabelInterpreter = 'LaTex';
title(ca_id)
legend(lines, {'averaged trial' 'fit' 'baseline from fit'}, 'Location', 'northwest')



% tuning
figure;polar_tuning2(tuning(:,:,ca_id),order([7:8 1:6]),2);   % rotated to final coord sys
applystyle(gca);
title(ca_id)

end


function applystyle(ax)
	fontsize = 12;

	ax.FontSize = fontsize;
end
