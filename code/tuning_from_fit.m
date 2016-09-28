function [tuning, tuning_onoff] = tuning_from_fit(coeffs, tuning_method)
% 1*8*634, 2*8*634

% exp or exp-exp

%g = generate_alphas(x(1:end-1-2*n), n, ti, method);
n_coeffs = size(coeffs, 2);
if n_coeffs==34
	fit_method = 'exp';
	n_conds = 8;  % 8 directions
elseif n_coeffs==35
	fit_method = 'exp-exp';
	n_conds = 8;  % 8 directions
elseif n_coeffs==7
	fit_method = 'exp-exp';
	n_conds = 1;  % averaged to single on-off
else
	error('unknown fit method for input fit')
end
n_peaks = 2 * n_conds;  % on and off

n_rois = size(coeffs, 1);

ti = coeffs(:, end-n_peaks+1:end);  % 2 or 16 for each ROI
As = coeffs(:, end-n_peaks-n_peaks+1:end-n_peaks);  % ditto
tau = coeffs(:, 1:2);  % if exp, only first column is valid

start_times = repmat(0:31:(n_conds-1)*31, 2, 1);
ti = ti - repmat(start_times(:).', n_rois, 1);  % such that t=0 means start of the condition

ti = ti.';
As = As.';
%tau

% we don't care about condition or cell now. (all time data have been aligned to the start of the condition)
ti = ti(:);
As = As(:);

t = 1:31;
t_ti = repmat(t, n_rois * n_peaks, 1) - repmat(ti, 1, 31);

tau = reshape(tau, 1, n_rois, 2);
tau = repmat(tau, n_peaks, 31, 1);
tau = reshape(tau, n_peaks * n_rois, 31, 2);


zero = t_ti<0;
t_ti(zero) = 0;

switch fit_method
case {'exp'}
		g = exp(-t_ti./tau(:,:,1));
		g(zero) = 0;
case 'exp-exp'
		g = exp(-t_ti./tau(:,:,1)) - exp(-t_ti./tau(:,:,2));
otherwise
	error('lllllllll')
end

gmax = max(g, [], 2) .* As;
garea = sum(g, 2) .* As;

tuning_onoff = reshape(gmax, 2, n_conds, n_rois);
tuning = sum(tuning_onoff, 1);

%{
[tuning, tuning_onoff] = tuning_from_fit(coeffs{3,2});
%}
