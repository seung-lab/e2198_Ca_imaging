function [params, cost, ispoor] = fit_single_onoff(yactual, method, debug_fit)


		ti_lb = [7.8 15.6] - 1;
		ti_ub = [7.8 15.6] + 1;
		ti_ub = [7.8 15.6] + 4;   % 82wo - extremely delayed rise / suppress by contrast?

		if 1 || ~strcmp(method, 'exp')

			%yactual = roi_sums_xcond_typemeans{celltype,:}(:);

			%{
			As = mean(onoff(:,:,ind),2);
			As = mean(onoff2(:,:,ind),2);
			%}
			As = find_max_min_diffs(yactual)
			initialguess = containers.Map();
			initialguess('exp-exp') = [4 1  -3 10 10  9 16]; % exp-exp
			initialguess('alpha') = [4  -5  5 5  9 16];
			initialguess('exp-exp') = [4 1  -3  As(1) As(2)  9 16]; % exp-exp
			%initialguess('exp-exp') = [14 2  -3  As(1) As(2)  9 16]; % exp-exp
			initialguess('exp-exp') = [15 1  -3  As(1) As(2)  9 16]; % exp-exp
			initialguess('exp-exp') = [4 1  -3  As(1) As(2)  7.8 15.6]; % exp-exp
			initialguess('exp') = [30  -5  5 5  9 16];
			initialguess('exp') = [30  -15  As(1) As(2)  9 16];
			guess = initialguess(method);
			lb = -Inf(size(guess));
			lb(end-1:end) = guess(end-1:end)-2;	% constrain time
			lb(end-1:end) = ti_lb;	% constrain time

			%ms: new for single-onoff fit
			lb(end-3:end-2) = 0; %constrain A.    '8w' '4ow' '82wi'  '1ws'

			% not used
			%%{
			ub = Inf(size(guess));
			ub(end-1:end) = ti_ub;	% constrain time  % 296. 1ni 4on typemean
			%}

			%ff = @(x, xdata) (ftrial(x(1:4), 'alpha', x(5:6), 1));
			%ff = @(x, xdata) (ftrial(x(1:5), 'exp-exp', x(6:7), 1));
			ff = @(x, xdata) (ftrial(x(1:end-2), method, x(end-1:end), 1));
			%lb=[];
			%ub=[];
			%[guess, cost] = lsqcurvefit(ff, guess, [], yactual);
			%options = optimoptions('lsqcurvefit', 'FunctionTolerance', 1e-9, 'MaxFunctionEvaluations', 10000);
			%[guess, cost] = lsqcurvefit(ff, guess, [], yactual, lb, [], options);
			%[guess, cost] = lsqcurvefit(ff, guess, [], yactual, lb);
			[guess, cost] = lsqcurvefit(ff, guess, [], yactual, lb, ub);
			yfit = ff(guess);
			%guess
			%cost
			cost_thres = (max(yactual)-min(yactual))/5
			if cost_thres>1 		% WARNING: THIS IS NOT PROPER for renormalized [0~1] traces..
				cost_thres = 2 * cost_thres.^2;
			else
				cost_thres = 2;		% WARNING: THIS IS NOT PROPER for renormalized [0~1] traces..
			end
			if cost>50 || cost>cost_thres
				warning(sprintf('poor single peak fit:  cost = %g, thres = %g', cost, cost_thres))
				ispoor = 1;
			else
				ispoor = 0;
			end

			if debug_fit
				figure; plot([yactual, yfit])
			end
		end

params = guess;
