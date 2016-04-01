fid = fopen('../spipe/ws.chunks/0/0/0/.seg');
seg = fread(fid, [5*63*63], 'uint32=>uint16');
seg = reshape(seg, [63 63 5]);
fclose(fid);
figure; imshow(seg(:,:,2))
colormap([0 0 0; rand(2^16,3)])
