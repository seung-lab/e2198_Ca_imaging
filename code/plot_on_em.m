
%{
%em_img = imread('em_gcl.png');

figure;
hold off
axes = imshow(em_img);

hold on
xx = (roi_centers(:,1) - 50) * 2230/450; % - 300;
yy = (roi_centers(:,2) -15) * 2500/515; % - 400;
tform = [   ]

em_slice_coords = tform .* 

plot(xx, yy, '*');
%}
