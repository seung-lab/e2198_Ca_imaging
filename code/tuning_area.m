function [roi_sums_area, roi_sums_area_fakeonoff] = tuning_area(roi_sums_means)

% roi_sums_means_flatten
% roi_sums_means  % 31*8*634
% e.g. #090 #190 #282 in the 5ti group

baselines = min(roi_sums_means(1:20, :, :), [], 1);
rebased = roi_sums_means - repmat(baselines, size(roi_sums_means,1), 1, 1);

%find( sum(sum(rebased<0), 2) )	% cell IDs with <0 responses

% TODO: Only sum over after 9th frame?
%rebased(rebased<0) = 0;  % get rid of <0 values
roi_sums_area = sum(rebased, 1);
roi_sums_area_fakeonoff = [roi_sums_area; roi_sums_area];

ca_id = 090;
ca_id = 051;
figure; plot(reshape(rebased(:,:, ca_id), [], 1))
h = gca();
h.XTick = sort([9:31:320 16:31:320 32:31:320]);
grid on

%{
[ordered, order] = sort(str2num(char(angles)));
ca_id = 090;
figure;polar_tuning2(roi_sums_area_fakeonoff(:,:,ca_id), order);
%}
