function [x, yfit, cost, residual] = lsqsearch(func, guess, y_actual)

	f = @(x) (func(x)-y_actual);
	ff = @(x, xdata) (func(x));
	opts = optimset('MaxFunEvals',50000, 'MaxIter',30000);

	[x, resnorm,residual] = lsqnonlin(f, guess); %, opts);
	%[x, resnorm,residual] = lsqnonlin(f, guess, [], [], opts);
	%{
	xx =x;
	[x, resnorm,residual] = lsqcurvefit(ff, guess, [], y_actual); %, opts);
	isequal(xx,x)
	%}

	yfit = func(x);
	%sum((yfit-y_actual).^2);
	cost = resnorm;
end