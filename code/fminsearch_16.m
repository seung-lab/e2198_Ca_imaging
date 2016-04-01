function [x, yfit, cost] = fminsearch_16(guess, y_actual)

	opts = optimset('MaxFunEvals',50000, 'MaxIter',30000);

	[x, cost] = fminsearch(@costf, guess, opts);

	function c = costf(x)
		%ypred = sum_exp(x(1), x(2:3), x(4:5), x(6:7), 1:length(y_actual));
		%ypred = sum_exp(x(1), x(2:3), x(4:5), x(6:7), [0:0.1:60]);
		ypred = f(x);
		c = (y_actual - ypred).^2; % least mean squares
		%c = abs(y_actual - ypred); % mean abs error
		c = sum(c(:));
	end

	function ypred = f(x)
		ti = [9 16  41 49  70 79  102 110  133 142  164 172  196 204  226 234];
		g = generate_alphas(x(1), 8, ti);
		%g = generate_alphas(x(1));
		coeffs16 = x(2:end);
		%g = generate_alphas(x(1:2));
		%coeffs16 = x(3:end);

		c0 = ones(1, 31*8);

		ypred = (coeffs16 * [c0; g]).';
	end

	n = length(y_actual);
	yfit = f(x);

end