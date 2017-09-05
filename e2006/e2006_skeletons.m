function [cell_info, A, strat_profile, skels] = e2006_skeletons(A)
	% input: optional, the transform to apply

dirpath = '~/dev/e2198-gc-analysis/e2006/Source-Version/';

if 0 && ~exist('A', 'var')
	search_pattern = [dirpath 'amf*.swc'];
else
	%search_pattern = [dirpath 'g*.swc'];
	search_pattern = [dirpath '*.swc'];
end

search_result=dir(search_pattern);
if isempty(search_result)
    error('cell not found');
end


figure; subplot(6,6,1)
cell_info = [];
ii = 0;
for file = search_result(:).'
	ii = ii+1;

	fname = file.name;
	cname = fname(1:end-4);
	c = [];
	c.annotation = cname;
	if strncmp(cname, 'amfon', 5)
		%c.type = 'on sac';
		c.type = 'ON SAC';
		c.class = 2;
	elseif strncmp(cname, 'amfoff', 5)
		%c.type = 'off sac';
		c.type = 'OFF SAC';
		c.class = 2;
	elseif strncmp(cname, 'gao', 3)
		c.type = 'dsgc';
		c.class = 1;
	else
		c.type = 'e2006gc';
		c.class = 1;
	end
	a = dlmread(fname, ' ', 4, 2);
	c.skel = a(:, 1:3);   % maybe change coord to match e2198? wait it already matches
	c.skel = [c.skel ones(size(c.skel(:,1)))];  % make it homogeneous coordinates
	%id = -ii - 10000;
	id = ii + 10000;
	c.cell_id = id;

	cell_info = [cell_info; c];

	if ii > 36
		continue;
	end
	subplot(6,6,ii);
	plot(a(:,1), a(:,2), '.', 'MarkerSize', 1)
	xlim([70, 120])
end

a = vertcat(cell_info.skel);

if ~exist('A', 'var')
onsacs = get_cell_info(cell_info, 'ON SAC');
onsacs = vertcat(onsacs.skel);
onsacs = onsacs(onsacs(:,1)<110, :);  % basic cropping
offsacs = get_cell_info(cell_info, 'OFF SAC');
offsacs = vertcat(offsacs.skel);
offsacs = offsacs(offsacs(:,1)>75, :);  % basic cropping

A = affine_e2006(onsacs, offsacs)
end

a = a * A;

cell_info_new = [];
for elem = cell_info(:).'
	elem.flatskel = elem.skel * A;
	x = [-20:120].';
	strat = histcounts(elem.flatskel(:,1), x);
	strat = strat(:);
	x = x(1:end-1)+0.5;
	elem.strat_unrml = [x strat];
	elem.strat_nrml = [x strat/sum(strat)];

	%{
	% doesn't do normalization.. wait actually does
	elem.cell_id = elem.cell_id + 1e7;
	strat_nrml = cell_info_bin_strat(elem, 1);
	strat_nrml = strat_nrml{elem.cell_id};
	assert(isequal(x, strat_nrml(:,1)))
	elem.strat_nrml = strat_nrml;
	elem.cell_id = elem.cell_id - 1e7;
	%}

	cell_info_new = [cell_info_new; elem];
end
cell_info = cell_info_new;

figure; plot(a(:,1), a(:,2), '.', 'MarkerSize', 0.1)

%figure;plot3(a(:,1), a(:,2), a(:,3), '.')
