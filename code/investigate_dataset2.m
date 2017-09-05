data_fn = [
    '011_p1_g6_DSbars_200um_export.mat';
    '012_p1_g6_DSbars_200um_export.mat';
    '013_p1_g6_DSbars_200um_export.mat';
    '014_p1_g6_DSbars_200um_export.mat';
    '015_p1_g6_DSbars_200um_export.mat';
    '016_p1_g6_DSbars_200um_export.mat';
    '017_p1_g6_DSbars_200um_export.mat';
    '018_p1_g6_DSbars_200um_export.mat';
    '019_p1_g6_DSbars_200um_export.mat';];
%{
data_fn = ['035_p2_g7_DSbars_200um_export.mat';
           '036_p2_g7_DSbars_200um_export.mat';
           '037_p2_g7_DSbars_200um_export.mat';
           '038_p2_g7_DSbars_200um_export.mat';
           '039_p2_g7_DSbars_200um_export.mat';
           '040_p2_g7_DSbars_200um_export.mat';
           '041_p2_g7_DSbars_200um_export.mat';
           '042_p2_g7_DSbars_200um_export.mat';
           '034_p2_g7_DSbars_200um_export.mat'];
%}

roi_tile = [];
for count = 1:size(data_fn,1)
       
    load(data_fn(count,:));
    [q,w] = find(sum(roi_sums,1)~=0);

    roi_tile(w) = count;
end



figure;

ncols = 5;
nrows = 3;
%startId = 214; % tile 4
%startId = 365;  % tile 6
%startId = 494;  % tile 8
startId = 562;  % tile 9
startId = 1;  % tile 1

stepsize = 4;
tile = 9;
w = find(roi_tile == tile);
for ii = 0:14

	subplot(nrows, ncols, 1+ii);
	%id = startId + ii*stepsize;
	id = w(1+ii*stepsize);
	%plot([roi_sums_all(:,id) imgaussfilt(roi_sums_all(:,id), 8/0.128)])
	plot([roi_sums_all(:,id) imgaussfilt(roi_sums_all(:,id), 20/0.128, 'Padding', 'symmetric')])
	title([id roi_tile(id)])
	xlim([0,1240])
	set(gca,'XTick', 0:31*8:1240);
	ax=gca;
	ax.XGrid = 'on';
end

% e2198:
% 166 in 8
