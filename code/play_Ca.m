function play_Ca()
% actual frame rate is about 8 fps (128ms/frame).

rawdata_fn = ['034_p2_g7_DSbars_200um.cfd';
           '035_p2_g7_DSbars_200um.cfd';
           '036_p2_g7_DSbars_200um.cfd';
           '037_p2_g7_DSbars_200um.cfd';
           '038_p2_g7_DSbars_200um.cfd';
           '039_p2_g7_DSbars_200um.cfd';
           '040_p2_g7_DSbars_200um.cfd';
           '041_p2_g7_DSbars_200um.cfd';
           '042_p2_g7_DSbars_200um.cfd'];

col = [2 3 3 2 1 1 1 2 3];
row = [2 2 3 3 3 2 1 1 1];

% Read raw CFD file
[dx,dy,channels,frames,ch1,ch2,t_scan,t_frame,ch1_t_stim] = read_cfd(rawdata_fn(1,:));
ch1 = permute(ch1,[2 1 3]);

frames = 1325;	% longest
bytes_stim_info = 3;
dx = dx-bytes_stim_info;
alltiles = uint8(ones(3*dy+3, 3*dx+3, frames)) * 150;

%{
tifff = 'AVG_043_p2_g7_ZoomedOut_center.tif';
tiff = imread(tifff);
tiffc = tiff(:,4:end);  % first 3 columns black
tiff2 = tiffc(1:2:end,1:2:end);
%}

for ii = 1:9
	[dx,dy,channels,frames,ch1,ch2,t_scan,t_frame,ch1_t_stim] = read_cfd(rawdata_fn(ii,:));
	size(ch1)
	ch1 = permute(ch1,[2 1 3]);
	max(ch1(:))
  dx = dx-bytes_stim_info;
  ch1 = ch1(:,bytes_stim_info+1:end,:);

%{
  maxi = max(ch1, [], 3);
  maxi = maxi(:,4:end);  % first 3 columns black
  xc = xcorr2(double(tiff2), double(maxi));
  figure; imshow(xc); set(gca,'CLim',[0 max(xc(:))])
%}

	scaling = 1;
	%alltiles(dy*row(ii)-dy+1 : dy*row(ii), dx*col(ii)-dx+1 : dx*col(ii), 1:frames) = ch1*(300/180);
	alltiles(dy*row(ii)-dy+1 + row(ii) : dy*row(ii) + row(ii), dx*col(ii)-dx+1 + col(ii) : dx*col(ii) + col(ii), 1:frames) = ch1*scaling;
  alltiles(dy*row(ii)-dy+1 + row(ii) : dy*row(ii) + row(ii), dx*col(ii)-dx+1 + col(ii) : dx*col(ii) + col(ii), frames+1:end) = 0;
end

%alltiles = permute(alltiles,[2 1 3]);

scaling = 500/180;
%implay(alltiles * scaling);

%%{
v = VideoWriter('newfile');
open(v)

%writeVideo(v,rand(300,300,1,5))
%class(alltiles * scaling)
%alltiles = ch1;
sz = size(alltiles);
sz = [sz(1:2) 1 sz(3)]
writeVideo(v,reshape(alltiles * scaling, sz))
close(v)
%}


maxi = max(alltiles, [], 3);
max(maxi(:))
figure; imshow(maxi*1.2);
%imwrite(maxi, 'max_tiles.tiff')
%set(gca,'CLim',[0 40])

end
%}


%{
ch1 = ch1(20:29,11:21,:);
  scaling = 500/180;
  ch1 = ch1 * scaling+70;
  implay(ch1);
  %}
