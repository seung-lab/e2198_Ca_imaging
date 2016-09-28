function [onoff, onoffmaxmins]= find_max_min_diffs(roi_sums_means)

	%roi_sums_means:  31*8*634, or 31*8 or 31*1
	t1 = 8.9468; t2 = 16.5457;
	t1 = 9; t2 = 16;

	dims = size(roi_sums_means);
	if length(dims)<3
		dims(3) = 1;
	end
	dims(1) = 2;

	onoff = NaN(dims);
	for k = 1:dims(3)   % 634
		for a = 1:dims(2)   % 8
			%offset = 31*(a-1);
			on = roi_sums_means(t1-5:t1+4, a, k);
			off = roi_sums_means(t2-2:t2+8, a, k);
			%{
			on = max(on) - min(on);
			off = max(off) - min(off);
			%}
			% Find max amount of increase (no decrease)
			on = max(max(tril(  repmat(on, 1, length(on)) - repmat(on.', length(on), 1)  )));
			off = max(max(tril(  repmat(off, 1, length(off)) - repmat(off.', length(off), 1)  )));
			onoff(:, a, k) = [on; off];
		end
	end
end
