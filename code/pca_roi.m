%pca

[dx,dy,channels,frames,ch1,ch2,t_scan,t_frame,ch1_t_stim] = read_cfd(rawdata_fn(4,:));
ch1 = permute(ch1,[2 1 3]);

X = permute(ch1,[3 1 2]);
sz = size(X);
X = reshape(X, [sz(1) sz(2)*sz(3)]);
coeff = pca(double(X));

max(coeff(:))
coeff = coeff / max(coeff(:)) / 0.3;
size(X)
size(coeff)
comp_img = reshape(coeff, sz(2), sz(3), []);
implay(comp_img)