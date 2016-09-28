function [ypred, CAg] = ftrial(x, varargin)

	nvarargin = length(varargin);
	ti = [];
	%90
	%ti = [9 16  41 49  70 79  102 110  133 142  164 172  196 204  226 234]+0.5;
	%346
	%ti = [9 17  41 49  71 79  103 111  134 141  164 173  196 204  226 234]+0.5;
	optargs = {'exp', ti, 8};
	optargs(1:nvarargin) = varargin;
	[method, ti, n] = optargs{:};

		g = generate_alphas(x(1:end-1-2*n), n, ti, method);
		coeffs16 = x(end-2*n:end);
	%{
	switch method
	case 'exp'
		g = generate_alphas(x(1), n, ti, method);
		coeffs16 = x(2:end);
	case 'alpha'
		g = generate_alphas(x(1), n, ti, method);
		coeffs16 = x(2:end);
	case 'exp-exp'
		g = generate_alphas(x(1:2), n, ti, method);
		coeffs16 = x(3:end);
	otherwise
		error('unknown method')
	end
	%}

	c0 = ones(1, 31*n);

	ypred = (coeffs16 * [c0; g]).';
	if nargout>1
		[~, g_raw] = generate_alphas(x(1:end-1-2*n), n, ti, method);
		c0 = ones(1, size(g_raw,2));
		CAg = repmat(coeffs16.', 1, length(c0)) .* [c0; g_raw];
	end
end