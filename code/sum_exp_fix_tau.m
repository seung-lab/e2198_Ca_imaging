function y = sum_exp_fix_tau(c, tau, A, t_i, t)
% f = c + sum A_i * g((t-t_i)/tau)
	t = t(:).';	 %column vec
	tt = zeros(length(t), length(t_i));
	for k = 1:length(t_i)
		tt(:, k) = (t - t_i(k));
	end
	exps = arrayfun(@alpha, tt);
	y = c + exps * A(:);

	function y = alpha(t)
		if t>=0
			y = t * exp(-t/tau);
		else
			y = 0;
		end
	end
end
