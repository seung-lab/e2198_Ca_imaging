function [onoff, onoffmaxmins]= find_diff_peaks(roi_sums_means)

	%roi_sums_means:  31*8*634
	t1 = 8.9468; t2 = 16.5457;
	t1 = 9; t2 = 16;

% too broad, e.g. #605
% 273, 21, how to fix?
% 41, 3rd cond, how?or any?
	onoff = NaN(2, 8, 634);
	for k = 1:634
		for a = 1:8
			%offset = 31*(a-1);
			on = roi_sums_means(t1-5:t1+4, a, k);
			off = roi_sums_means(t2-2:t2+8, a, k);
			on = diff(on);
			off = diff(off);
			[val, ind] = max(on);
			on = ind-1 + t1-5 ;
			[val, ind] = max(off);
			off = ind-1 + t2-2 ;
			onoff(:, a, k) = [on; off];
		end
	end
end
