%function show_Ca()
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
[dx,dy,channels,frames,ch1,ch2,t_scan,t_frame,ch1_t_stim] = read_cfd(rawdata_fn(4,:));
ch1 = permute(ch1,[2 1 3]);

corrx = zeros(size(ch1,1), size(ch1,2));
corry = zeros(size(ch1,1), size(ch1,2));
meanall = mean(ch1, 3);
%max(ch1(:))
%max(meanall(:))
%min(meanall(:))
figure; imshow(uint8(meanall*300/180)); title('mean');

mini = min(ch1,[],3);
figure; imshow(mini*(300/180));
 title('min');
%subtract = 

%%{
for ii = 1:size(ch1, 1)-1
  for jj = 1:size(ch1, 2)-1
    a = ch1(ii,jj,:); b = ch1(ii+1,jj,:);
    corry(ii,jj) = corr(double(a(:)), double(b(:)));
    a = ch1(ii,jj,:); b = ch1(ii,jj+1,:);
    corrx(ii,jj) = corr(double(a(:)), double(b(:)));
  end
end
corry(isnan(corry)) = 0;
corrx(isnan(corrx)) = 0;

img = [ones(1, size(corry, 2)); corry(1:end-1,:)];
affx = zeros(size(img)-1);
affy = zeros(size(img)-1);
for ii = 1:size(img, 1)-1
  for jj = 1:size(img, 2)-1
    affx(ii,jj) = min(img(ii,jj), img(ii,jj+1));
    affy(ii,jj) = min(img(ii,jj), img(ii+1,jj));
  end
end
corryaffx = affx;
corryaffy = affy;

%}
%img = corry;
%img = double(mini) / double(max(mini(:)));
img = meanall / max(meanall(:));
affx = zeros(size(img)-1);
affy = zeros(size(img)-1);
affz = ones(size(img)-1);
for ii = 1:size(img, 1)-1
  for jj = 1:size(img, 2)-1
    affx(ii,jj) = min(img(ii,jj), img(ii,jj+1));
    affy(ii,jj) = min(img(ii,jj), img(ii+1,jj));
  end
end
affx = max(0, affx);
affy = max(0, affy);
affx = affx + corry(1:end-1,1:end-1);
affy = affy + corry(1:end-1,1:end-1);
%affx = affx + corryaffx;
%affy = affy + corryaffy;
affx = affx / max(affx(:));
affy = affy / max(affy(:));
%affy = corry(1:end-1,1:end-1);
affx = affy;

figure; imshow(affx); title('xtmp')
figure; imshow(affy); title('ytmp')
%aff = cat(4, affz, affy, affx);
aff = cat(4, affx, affy, affz);
aff = repmat(aff, [1 1 5 1]);
filename = 'aff_corr.h5';
h5create(filename,'/main',[size(img)-1, 5, 3], 'Datatype', 'single')
h5write(filename,'/main',single(aff))


%{
figure; imshow(corrx); title('in x');
figure; imshow(corry); title('in y');
corrx = watershed(-corrx(:,4:end));
corry = watershed(-corry(:,4:end));
figure; imshow(uint8(corrx)); title('in x');
figure; imshow(uint8(corry)); title('in y');
%}
%{
ch1 = ch1(20:29,11:21,:);
  scaling = 500/180;
  ch1 = ch1 * scaling+70;
  implay(ch1);
  %}

%end