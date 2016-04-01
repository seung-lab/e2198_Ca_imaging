function [x, val] = fminsearch_expo(y_actual)

	opts = optimset('MaxFunEvals',50000, 'MaxIter',10000);

	% 3 independant exp fit
	guess = [-1 4  5 3 3  9 16 -15]
	% >>    -2.2377    3.0902    4.2131    3.5008  1.6123e+10    8.9416   16.5448  -88.0457

	% 2 exp fit
	guess = [-1 4  5 3  9 16]
	%    -2.0722    2.9869    4.2363    3.6055    8.9821   16.5538	% 2 exp without initial slope
	%    -2.2217    3.0803    4.2113    3.5119    8.9468   16.5457	% 2+1 exp with initial slope
	% tau = 3.0803
	% t1 = 8.9468, t2 = 16.5457
	% g(t) = t * exp(-t/tau)
	% y = c + A1 * g(t - t1) + A2 * g(t - t2) + A2 * g(t - t2 + 31)

	guess = [-5 4  5 5  9 16]

	% 2*(exp - exp)
	guess = [-1   4 1   5 3  9 16]
	guess = [-3   4 1   10 10  9 16]
	%-2.7218    6.9411    0.9888    8.5424    6.2784    9.3993   16.7412   cost =    0.5372
	%guess = [-3   6 4   10 17  9 16] or
	%guess = [-2   10 2   10 6  9 16]
	% -2.6667    6.1073    1.5379   10.8618    8.1619    8.9445   16.5903  cost = 0.8679 %fminsearch only, or lsq with the next guess below
	%guess = [-3   4 1   4 4]%  8.9*1e5 16.5*1e5]
	% -2.6494    6.0481    1.5665   10.9964    8.2793   cost =0.8845
	%guess = [-13   4 1   10 10  9.2 16.5]

	[x, val] = fminsearch(@cost, guess, opts);

	%guess = [-5 4  5 3  9 16]
	ff = @(x, xdata) (f(x));
	[x, val] = lsqcurvefit(ff, guess, [], y_actual);
	%{
	[guess, cost] = lsqcurvefit(ff, guess, [], yactual);
	coeffs16tau(ind,:) = guess;
	yfit = ff(guess);
	%}

	function c = cost(x)
		%ypred = sum_exp(x(1), x(2:3), x(4:5), x(6:7), 1:length(y_actual));
		%ypred = sum_exp(x(1), x(2:3), x(4:5), x(6:7), [0:0.1:60]);
		ypred = f(x);
		c = (y_actual - ypred).^2; % least mean squares
		%c = abs(y_actual - ypred); % mean abs error
		c = sum(c(:));
	end

	function ypred = f(x)
		% 3 exp, independant
		%ypred = sum_exp_fix_tau(x(1), x(2), x(3:5), x(6:8), 1:length(y_actual));

		%{
		% 2 exp
		ypred = sum_exp_fix_tau(x(1), x(2), x(3:4), x(5:6), 1:length(y_actual));

		% +1 exp same as 2nd
		ypred = ypred + sum_exp_fix_tau(0, x(2), x(4), x(6)-length(y_actual), 1:length(y_actual));
		%}

		%%{
		% 2*(exp - exp)
		
		%g = generate_alphas(x(2:3), 1); %, x(6:7)/1e5);
		g = generate_alphas(x(2:3), 1, x(6:7), 'exp-exp');
		As = x(4:5);

		ypred = x(1) + (As * g).';
		%}

		%{
		% 2 alpha + 2 alpha
		g = generate_alphas(x(2), 1, x(5:6), 'alpha');
		As = x(3:4);

		ypred = x(1) + (As * g).';
		%}
	end

	n = length(y_actual);
	ypred = f(x);
	figure; plot(1:n, [y_actual, ypred])

end





% dx + x = delta(t)
% t * exp(-t)