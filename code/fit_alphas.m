tau = 3.0803;
t1 = 8.9468; t2 = 16.5457;
% g(t) = t * exp(-t/tau)
% y = c + A1 * g(t - t1) + A2 * g(t - t2) + A2 * g(t - t2 + 31)
%   = c + A1 * g1 + A2 * g2

t = 1:31;
t_t1 = t - t1;  t_t1(t_t1<0) = 0;
t_t2 = t - t2;  t_t2(t_t2<0) = 0;
t_t2_ = t - t2 + 31;

g1 = t_t1 .* exp(-t_t1/tau);
g2 = t_t2 .* exp(-t_t2/tau) + t_t2_ .* exp(-t_t2_/tau);

c0 = ones(1, 31);

% yfit.' = [c A1 A2] * [c0; g1; g2]
% [c A1 A2] = yactual.' / [c0; g1; g2];

% yfit = [c0; g1; g2].' * [c A1 A2].'
% [c A1 A2].' =  [c0; g1; g2].' \ yactual;
divider = [c0; g1; g2];
%denominator


roi_sums_all_reshape = reshape(roi_sums_all, 31, 8, 5, 634);

roi_sums_all_reshape = roi_sums_all_reshape(:, :, 2:5, :);
roi_sums_means = squeeze(mean(roi_sums_all_reshape, 3));


% c, A1, A2
coeffs = reshape(roi_sums_means, 31, 8 * 634).' / divider;
coeffs = reshape(roi_sums_means, 31, 8 * 634).' / divider;
roifit = (coeffs * divider).';

[ordered, order] = sort(str2num(char(angles)));

coeffs_reshape = reshape(coeffs.', 3, 8, 634);
coeffs_ordered = coeffs_reshape(:, order, :);


ind = 90;

figdir = '~/dev/e2198_Ca_imaging/cell_summary';
for ind = c5ti
	responses = coeffs_ordered(:, :, ind);
	tmp = [responses(:, 1) responses];	%BUG: end, not 1
	range = (0:8)*pi/4;
	%range = range.';
	figure;
	polar(range, tmp(2, :), '');
	hold on;
	polar(range, tmp(3, :), '');
	polar(range, tmp(1, :), '--');
	title(num2str(ind));

	print(gcf, '-r300', sprintf('%s\\%s.png',figdir,num2str(ind)), '-dpng');
end