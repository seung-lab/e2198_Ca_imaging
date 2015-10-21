function [dx,dy,channels,frames,ch1,ch2,t_scan,t_frame,ch1_t_stim] = f(fname)

info = imfinfo(fname);

dx = info.Width;
dy = info.Height;
channels = 1;
frames = length(info);

data = imread(fname);
size(data)
data = rot90(flipud(data),-1);
ch1 = reshape(data,[dx dy frames]);
ch2 = [];

t_scan = [0:frames*dy-1]*2; %2ms per line
ch1_t_stim = reshape(ch1(1,:,:),[dy*frames 1]);
t_frame = [0:frames-1]*2*dy;


