function [dx,dy,channels,frames,ch1,ch2,t_scan,t_frame,ch1_t_stim] = f(fname)

fid = fopen(fname);
header = fread(fid,768/2,'uint16');
dx = header(141);
dy = header(142);
channels = header(149);
frames = header(147);

if(channels>1)
    channels = 2;
    data = fread(fid,dx*dy*frames*channels,'*uint8');
    ch1 = reshape(data(1:channels:end),[dx dy frames]);
    ch2 = reshape(data(2:channels:end),[dx dy frames]);
else
    channels = 1;
    data = fread(fid,dx*dy*frames*channels,'*uint8');
    ch1 = reshape(data,[dx dy frames]);
    ch2 = [];
end

t_scan = [0:frames*dy-1]*2; %2ms per line
ch1_t_stim = reshape(ch1(1,:,:),[dy*frames 1]);
t_frame = [0:frames-1]*2*dy;

% ch1=ch2;
ch1(1:3,:,:) = NaN;

fclose(fid);

