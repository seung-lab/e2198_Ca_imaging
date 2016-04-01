function out = vali(varargin)

    out = {};
	for x = varargin
        x = x{1};
		if any(x)
			out{end+1} = x;
		else
			%out{end+1} = [0 0 0];
            out{end+1} = [NaN NaN NaN];
            %out{end+1} = 1e3 * [1.5493   -0.1904   -1.0847];
		end
	end
end