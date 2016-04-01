function ti = t1t2_to_ti(t1, t2, varargin)  % n=8


	nvarargin = length(varargin);
	optargs = {8};
	optargs(1:nvarargin) = varargin;
	[n] = optargs{:};

	t1s = (0:n-1)*31 + t1;
	t2s = (0:n-1)*31 + t2;
	ti = [t1s; t2s];
	ti = ti(:);

end