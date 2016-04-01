tau = 3.0803;
t1 = 8.9468; t2 = 16.5457;
% g(t) = t * exp(-t/tau)
% y = c + A1 * g(t - t1) + A2 * g(t - t2) + A2 * g(t - t2 + 31)
%   = c + A1 * g1 + A2 * g2

t = 1:31*8;
t1s = (0:7)*31 + t1;
t2s = (0:7)*31 + t2;
ti = [t1s; t2s];
ti = ti(:);


g = zeros(8*2, length(t));
for k = 1:8*2
t_ti = t - ti(k);  t_ti(t_ti<0) = 0;
g(k,:) = t_ti .* exp(-t_ti/tau);
end

% assuming last peak causes the initial downslope at the beginning
t_t2_ = t - t2 + 31;
g(end,:) = g(end,:) + t_t2_ .* exp(-t_t2_/tau);


c0 = ones(1, 31*8);

% yfit.' = [c A1 A2 ...] * [c0; g1; g2; ...]
% [c A1 A2] = yactual.' / [c0; g1; g2];

divider = [c0; g];
%denominator


% c, A1, A2, ...
%{
clip = [floor(t1)-1:t1+3, floor(t2)-1:t2+3];
clip2 = reshape(repmat(0:31:31*7, 10, 1), 1, []) + reshape(repmat(clip, 1, 8), 1, []);
clipped = roi_sums_means(clip, :, :);
clippeddivider = divider(:, clip2);
coeffs16clipped = reshape(clipped, 10 * 8, 634).' / clippeddivider;	% 634x17
fit16 = (coeffs16clipped * divider).';
%}

%{
coeffs16 = reshape(roi_sums_means, 31 * 8, 634).' / divider;	% 634x17
fit16 = (coeffs16 * divider).';
%}

[ordered, order] = sort(str2num(char(angles)));

%coeffs16_reshape = reshape(coeffs16.', 1+16, 634);
%coeffs16_ordered = coeffs16_reshape(:, order, :);

%%{
%coeffs16tau = [NaN(634, 1) NaN(size(coeffs16))];
%coeffs16expdiff = [NaN(634, 2) NaN(size(coeffs16))];
%fit16tau = NaN(31*8, 634);
%fit16expdiff = NaN(31*8, 634);
method_list = {'exp', 'alpha', 'exp-exp'}; %, 'exp_w_ti'};
fit_list = {'tau+As', 'tau+As+ti'};
%{
fit16 = cell(0);
for method = 1:length(method_list)
	coeffs16{method} = NaN(634, 1+17);
	fit16{method} = NaN(31*8, 634);
end
coeffs16{3} = NaN(634, 2+17);  % 'exp-exp'
%}

for fit_index = 1:2
if fit_index==2
	optimize_ti = true;
else
	optimize_ti = false;
end

for ind = 1:n_rois %74 %48 %41 %19 %5 %409 %3 %double checking 73 %67 %61 %22 %41  %1:634  %[346 90] %c5ti

	ind

	for method_id = [1 3] %1:length(method_list) % {'exp', 'alpha', 'exp-exp'}

		method = method_list{method_id};

		if ~strcmp(method, 'exp')
			% fit to cross condition(direction) mean first
			initialguess = containers.Map();
			initialguess('exp-exp') = [4 1  -3 10 10  9 16]; % exp-exp
			initialguess('alpha') = [4  -5  5 5  9 16];
			guess = initialguess(method);
			%ff = @(x, xdata) (ftrial(x(1:4), 'alpha', x(5:6), 1));
			%ff = @(x, xdata) (ftrial(x(1:5), 'exp-exp', x(6:7), 1));
			ff = @(x, xdata) (ftrial(x(1:end-2), method, x(end-1:end), 1));
			yactual = roi_sums_xcondmeans(:,ind);
			[guess, cost] = lsqcurvefit(ff, guess, [], yactual);
			yfit = ff(guess);
			%guess
			%cost
			%figure; plot([yactual, yfit])
			newguess = [guess(1:end-4) repmat(guess(end-3:end-2), 1, 8)];
			ti = guess(end-1:end);
			guess = newguess;
			lb = zeros(size(guess));
			lb(end-16) = -100;	% vertical intercept

		else  % exp
			guess = [30 -50 reshape(onoff(:,:,ind),1,16)];
			ti = maxdiffs_trial(:, ind) + 0.5;
			lb = [0 -100 zeros(1, 16)];
		end

		%guess = [tau 1 coeffs16(ind, :)*2];
		%guess = [tau 1 coeffs16(ind, 1) coeffs16(ind, 2:end)*2];
		%guess = [tau 0.1 coeffs16(ind, 1) ones(1, 16)*10]; % exp-exp
		%guess = [tau coeffs16(ind, :)];	% alpha
		%guess = [30 -50 coeffs16(ind, 2:end)];	 % exp
		%guess = [30 -50 reshape(onoff(:,:,ind),1,16)];
		%guess = ones(1, 18);


		yactual = roi_sums_means(:, :, ind);
		yactual = yactual(:);

		%[coeffs16tau(ind,:), yfit, cost] = fminsearch_16(guess, yactual);
		%[coeffs16tau(ind,:), yfit, cost] = lsqsearch(@ftrial, guess, yactual);


		if optimize_ti && strcmp(method, 'exp')

			ti_expected = floor(t1t2_to_ti(t1, t2)) + 0.5; %8.5 16.5 ...
			%ti_to_optimize = find((ti < ti_expected-1) | (ti > ti_expected+2));
			ti_to_optimize = 1:16;
			ti(ti_to_optimize) = ti_expected(ti_to_optimize);

			opts = optimset('display','off');
			ti_orig = ti;
			for k = ti_to_optimize(:).'
				ti_best = [];
				cost_best = 1e10;
				guess_best = [];
				for offset= -1:2 %[-1 1 2]
					ti(k) = ti_orig(k) + offset;
					ff = @(x, xdata) (ftrial(x, method, ti));	% need to do this to use the up to date "ti"
					[guess, cost] = lsqcurvefit(ff, guess, [], yactual, lb, [], opts);
					if cost < cost_best
						ti_best = ti;
						cost_best = cost;
						guess_best = guess;
					end
				end
				ti = ti_best;
			end
			coeffs16{method_id,fit_index}(ind,:) = [guess_best ti_best(:).'];
			ff = @(x, xdata) (ftrial(x, method, ti_best));	% need to do this to use the up to date "ti"
		else
			%ff = @(x, xdata) (ftrial(x, 'alpha'));
			%ff = @(x, xdata) (ftrial(x, 'exp-exp'));
			ff = @(x, xdata) (ftrial(x, method, ti));

			%[guess, cost] = lsqcurvefit(ff, guess, [], yactual);
			[guess_best, cost_best] = lsqcurvefit(ff, guess, [], yactual, lb);
			%coeffs16tau(ind,:) = guess;
			%coeffs16expdiff(ind,:) = guess;

			%cost_best

			if optimize_ti
				ti_expected = t1t2_to_ti(ti(1), ti(2));
				guess2 = [guess_best, ti_expected.'];
				%guess2 = [guess, ti_expected.'];
				ub = Inf(size(lb));
				ub = [ub ti_expected.'+2];
				lb = [lb ti_expected.'-1];
				ff = @(x, xdata) (ftrial(x(1:end-16), method, x(end-16+1:end)));
				[guess_best, cost_best] = lsqcurvefit(ff, guess2, [], yactual, lb);
			end
			coeffs16{method_id,fit_index}(ind,:) = guess_best;
		end
		yfit = ff(guess_best);

		cost_best
		%fit16tau(:,ind) = yfit;
		%fit16expdiff(:,ind) = yfit;
		fit16{method_id,fit_index}(:,ind) = yfit;
		fit16_cost{method_id,fit_index}(ind) = cost_best;

	end

end % for fit_index = 1:2

end
%}


%figure;plot([roi_sums_means_flatten(:,ind) fit16{1}(:,ind) fit16{3}(:,ind)])

