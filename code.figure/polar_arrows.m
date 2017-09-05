function lines = polar_arrows(theta, rho, bidir, varargin)
%e.g. 
%b = ca_dsos_allrois(get_ca_ids(cell_dict_j, cell_info, '82wi'), :);
%figure;polar_arrows(b.os_theta, b.os_r./b.r_mean, 2);

if 1    % to final "standard" coord
	theta = pi/2+theta;
end

if ~exist('bidir', 'var')
	bidir = 0;
end

if bidir
	tt = insert_val3(theta, theta);
	tt = insert_val3(tt, theta, 2);
	rr = insert_val3(0*rho, rho);
	rr = insert_val3(rr, -rho, 2);
	lines = polarplot(tt, rr, varargin{:});
else
	lines = polarplot(insert_val(0*theta, theta), insert_val(0*rho, rho), varargin{:});
end

function out = insert_val(in, to_insert)
	d = ndims(in);
	in = cat(d+1, in, to_insert);
	in = shiftdim(in, d);
	sz = size(in);
	in = reshape(in, sz(1)*sz(2), deal(sz(3:end)));
	out = in;

function out = insert_val2(in, to_insert)
	out = [shiftdim(-1,in); shiftdim(-1,to_insert)];
	% TODO

function out = insert_val3(in1, in2, pitch1, pitch2)
	if ~exist('pitch1', 'var')
		pitch1 = 1;
	end
	if ~exist('pitch2', 'var')
		pitch2 = 1;
	end
	sz1 = size(in1);
	in1 = reshape(in1, pitch1, [], deal(sz1(2:end)));
	sz2 = size(in2);
	in2 = reshape(in2, pitch2, [], deal(sz2(2:end)));
	out = [in1; in2];
	sz = size(out);
	out = reshape(out, sz(1)*sz(2), deal(sz(3:end)));
	