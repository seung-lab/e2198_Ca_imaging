function [onoff, onoffmaxmins]= find_max_min_diffs(roi_sums_means)

	%roi_sums_means:  31*8*634
	t1 = 8.9468; t2 = 16.5457;
	t1 = 9; t2 = 16;

	onoff = NaN(2, 8, 634);
	for k = 1:634
		for a = 1:8
			%offset = 31*(a-1);
			on = roi_sums_means(t1-5:t1+4, a, k);
			off = roi_sums_means(t2-2:t2+8, a, k);
			on = max(on) - min(on);
			off = max(off) - min(off);
			onoff(:, a, k) = [on; off];
		end
	end
end
