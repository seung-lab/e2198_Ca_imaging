function [ret] = ca_quality_control(roi_sums_all_reshape, framerange, method, threshold)
% returns: quality values, or idx for elements that passed.


	%cellwise_raw = roi_sums_all_reshape(:, : ,:, ca_ids);
	cellwise_raw = roi_sums_all_reshape;
	cellwise_raw = cellwise_raw(framerange, :, :, :);

	[frames, dirs, trials, cells] = size(cellwise_raw);
	assert(dirs==8 && trials==4 && frames<=31)


	switch method
	case 'none'
		ret = ones(cells,1);
		threshold = 0;
	case 'snr'
		cellwise_raw = reshape(cellwise_raw, frames*dirs, trials, cells);
		snr = var(mean(cellwise_raw, 2), [], 1) ./ mean(var(cellwise_raw,[],2), 1);
		ret = snr(:);
	case {'qi', 'Qi', 'quality index'}	% Baden 2016
		cellwise_raw = reshape(cellwise_raw, frames*dirs, trials, cells);
		qi = var(mean(cellwise_raw, 2), [], 1) ./ mean(var(cellwise_raw,[],1), 2);
		ret = qi(:);
	case {'mean_max-min', 'magnitude'}
		cellwise_xcondmeans = squeeze(mean(reshape(cellwise_raw, frames, dirs*trials, cells), 2));
		magnitude = max(cellwise_xcondmeans, [], 1) - min(cellwise_xcondmeans, [], 1);
		ret = magnitude(:);
	end

	if exist('threshold', 'var')
		ret = find(ret > threshold);
	end
