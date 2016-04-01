%{

roi_sums_all;
n_rois;
%n_rois = size(roi_sums_all, 2); %634

nconds; % 8 directions
frames_per_condition = 31;	  % t_pre_frames+t_stim_frames+t_post_frames
angles = {stim_struct.condnames};

roi_centers;
roi_borders;


run cell_mapping_verified.m
%omni_id_from_roi = NaN(n_rois, 1);
cell_dict;

%}


roi_sums_all_reshape = reshape(roi_sums_all, 31, 8, 5, 634);
roi_sums_all_reshape = roi_sums_all_reshape(:, :, 2:5, :);  % disgard first trial

roi_sums_means = squeeze(mean(roi_sums_all_reshape, 3));	% 31*8*634
roi_sums_means_flatten = reshape(roi_sums_means, 31*8, 634);

roi_sums_xcondmeans = squeeze(mean(roi_sums_means, 2));	% 31*634

xcellmeans4 = mean(reshape(roi_sums_means, 31*8,[]),2);
xcellmeans = mean(reshape(roi_sums_all, 31*8,[]),2);
ymean = mean(reshape(roi_sums_all, 31, []), 2);

% figure;plot()
% hold on; plot(ymean)


onoff = find_max_min_diffs(roi_sums_means);
maxdiffs = find_diff_peaks2(roi_sums_means);
maxdiffs_trial = NaN(size(maxdiffs));
for a = 1:8
	maxdiffs_trial(:, a, :) = maxdiffs(:, a, :) + 31*(a-1);
end
maxdiffs_trial = reshape(maxdiffs_trial, 2*8, 634);