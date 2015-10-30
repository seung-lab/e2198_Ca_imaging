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
alltiles = uint8(ones(3*dy+3, 3*dx+3, frames)) * 200;

for ii = 1:9
	[dx,dy,channels,frames,ch1,ch2,t_scan,t_frame,ch1_t_stim] = read_cfd(rawdata_fn(ii,:));
	size(ch1)
	ch1 = permute(ch1,[2 1 3]);
	%max(ch1(:))
	scaling = 300/180;
	%alltiles(dy*row(ii)-dy+1 : dy*row(ii), dx*col(ii)-dx+1 : dx*col(ii), 1:frames) = ch1*(300/180);
	alltiles(dy*row(ii)-dy+1 + row(ii) : dy*row(ii) + row(ii), dx*col(ii)-dx+1 + col(ii) : dx*col(ii) + col(ii), 1:frames) = ch1*scaling;
end

%alltiles = permute(alltiles,[2 1 3]);

implay(alltiles);

end