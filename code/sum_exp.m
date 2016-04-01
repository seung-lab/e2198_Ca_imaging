function y = sum_exp(c, A, t_i, tau_i, t)
% f = c + sum A_i * g((t-t_i)/tau_i)
	t = t(:).';	 %column vec
	tt = zeros(length(t), length(t_i));
	for k = 1:length(t_i)
		tt(:, k) = (t - t_i(k)) / tau_i(k);
	end
	exps = arrayfun(@expdelta, tt);
	y = c + exps * A(:);
end


function y = expdelta(t)
	if t>=0
		y = exp(-t);
	else
		y = 0;
	end
end
