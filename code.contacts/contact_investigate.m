function contact_stat = contact_investigate(cell_info, cellquery, feature, varargin)

self = @contact_investigate;

%{
if ~exist(args, 'var')
	args = {};
elseif ~iscell(args)
	args = {args};
end
%}


cells = get_cell_info(cell_info, cellquery);
cells = get_cell_info_table(cells);
cells.soma_mip2 = nan(size(cells,1),3);
%class(cells.soma_coords_warped_mip2_zscaled)
for ii = 1:size(cells,1)
	%	cells.soma_mip2(ii) = nan(1,3);
	%else
	if ~isempty(cells.soma_coords_warped_mip2_zscaled(ii))
		cells.soma_mip2(ii,:) = cells.soma_coords_warped_mip2_zscaled(ii,:);
		cells.soma_mip2(ii,2) = cells.soma_mip2(ii,2) *16.5/23;
	end
end

%basepath = '~/dev/e2198_Ca_imaging/data/contacts/raw';
basepath = 'contacts/raw2d_600-1167';
basepath = 'contacts/raw3d_445-1167';
surfacearea_file_path = [basepath '/surface_area.mat'];
load(surfacearea_file_path);  % var: 'surfacearea'

contacts = {};
%flat = [];

cellids = [];
%for c = cells(:).'
for cell_id = cells.cell_id(:).'
	%{
	search_pattern=sprintf('%s/%d.mat',basepath,c.cell_id);
	search_result=dir(search_pattern);
	if isempty(search_result)
	    error('cell not found');
	end
	file_path=sprintf('%s/%s',dir_sac2gc,search_result.name);
	%}

	file_path = sprintf('%s/%d.mat', basepath, cell_id);
	if ~exist(file_path, 'file')
		display(sprintf('cell not found: %d \n', cell_id));
		continue
	end
	tmp = load(file_path);  % var: 'contacts'
	contacts{cell_id} = tmp.contacts;
	cellids(end+1) = cell_id;
end
display(sprintf('cells:  \n  %s', num2str(cellids)))

flat = horzcat(contacts{cellids});  %note: turns out horzcat(contacts{:}) is much slower
flat0 = flat;
flat = [];
surface_area = [];
contact_stat = [];
for cell_id = cells.cell_id(:).'
	con = contacts{cell_id};
	con_count = size(con,2);
	flat = [flat [repmat(cell_id, 1, con_count); con]];

	contact_stat = [contact_stat [surfacearea(:, cell_id == surfacearea(1,:)); con_count]];	% cell_id, surface_area, contact_area, contact_area_of_interest
end
%class(contact_stat)	% int64
contact_stat = double(contact_stat);	%avoids all sorts of problems including silent conversion of double to int when accessing ranges in a table
contact_stat = array2table(contact_stat.', 'VariableNames', {'cell1_id','surface_area','contact_area'}); %,'percent_contact'});

	flat = flat.';
	if size(flat, 2) > 4  %3d
		flat = array2table(flat, 'VariableNames', {'cell1_id','cell_id','x','y','z'});
	else
		flat = array2table(flat, 'VariableNames', {'cell1_id','cell_id','x','y'});
		flat.z = double(flat.z);
	end
	% Getting around this error: Integers can only be combined with integers of the same class, or scalar doubles.
	flat.x = double(flat.x);
	flat.y = double(flat.y);
npoints = size(flat0, 2);

xyz = [2; 3; 1];
ijk = [3; 1; 2];
resolution = [16.5; 16.5; 23];
resolution_xy = resolution(xyz(1:2));

fig = figure;
somacolor = 'k';
	if npoints > 1e5	% for speed
		markertype = '.';
	else
		markertype = 'o';
	end

switch feature
case {'sacdir'}
	nvarargin = length(varargin);
	optargs = {3};
	optargs(1:nvarargin) = varargin;
	[onoff] = optargs{:};

	sac_somas = 'sac_soma_centers_m2_warped.csv';
	sac_somas = csvread(sac_somas);
	sac_somas = array2table(sac_somas, 'VariableNames', {'cell_id','sac_i','sac_j','sac_k'}); % x,y,z in omni coord directions
	sac_somas = sac_somas(:, {'cell_id','sac_j','sac_k'}); % light axis removed
	
	onsacs = struct2table(get_cell_info(cell_info, 'ON SAC'));
	onsacs.onoff(:,1) = 1;
	offsacs = struct2table(get_cell_info(cell_info, 'OFF SAC'));
	offsacs.onoff(:,1) = 2;
	sacs = [onsacs; offsacs];
	sacs = sacs(:, {'cell_id', 'onoff'});
	sac_somas = join(sac_somas, sacs);

	[flat, ~, ib] = outerjoin(flat, sac_somas, 'Type', 'left'); % left join
	flat.sacdir = rad2deg(atan2( (flat.y-flat.sac_k)*resolution_xy(2), (flat.x-flat.sac_j)*resolution_xy(1) ));
	%flat.sacdir(isnan(flat.sacdir)) = -40;
	flat = flat(~isnan(flat.sacdir), :);

	% on or off
	flat = flat(~~bitand(flat.onoff, onoff), :);

	flat.Properties.VariableNames{'sacdir'} = 'value';

	scatter([],[]); hold on; % initiate view as 2D plot
	scatter(cells.soma_mip2(:,1), cells.soma_mip2(:,2), somacolor, 'o');
	colormap('hsv');
	%scatter3(flat.x, flat.y, flat.value, [1], flat.value, markertype);
	%scatter3(flat.x, flat.y, flat.value, [1], flat.value, 'filled', markertype);
	scatter(flat.x, flat.y, [1], flat.value, markertype);
	hold off;

	colorbar();
	title(cellquery)

case {''}
	nvarargin = length(varargin);
	optargs = {1};
	optargs(1:nvarargin) = varargin;
	[colorschemeIdx] = optargs{:};

	colorscheme = {'default', 'colorcube'};
	colorscheme = colorscheme{colorschemeIdx};

	somacolor = 'r';

	%cellstr(char(cell_info.type))
	%allcelltypes = categories({cell_info.type});
	
	%uniqueIds = unique(flat0(1,:));
	%uniquepartners = get_cell_info(cell_info, uniqueIds);

	if 1 % version 1
	cell_info_table = get_cell_info_table(cell_info);
	%unique(cell_info_table.type)  % empty string converts <undefined>
	cell_info_table.type = categorical(cell_info_table.type);

	%length(find(isnan(double(cell_info_table.type))))	%  <undefined> converts to NaN
	typevalue = double(cell_info_table.type);	%  <undefined> converts to NaN
	typevalue(isnan(typevalue)) = 0;
	cell_info_table.value = cell_info_table.class*100 + typevalue;
	%type_table = cell_info_table(:, {'cell_id', 'class', 'type'});
	value_table = cell_info_table(:, {'cell_id', 'value'});

	[flat, ~, ib] = outerjoin(flat, value_table, 'Type', 'left'); % left join
	% Add value for cells not found
	%length(find(flat.value(~ib)))
	%flat.value(~ib) = min(flat.value)-100;
	flat.value(~ib) = 0;
	flat.value(flat.value<100) = 100;
	flat.value(flat.value>=400) = 80; % 
	%length(find(isnan(flat.value)))
	
	%scatter(flat.x, flat.y, [], flat.value);
	plot([],[]); hold on; % initiate view as 2D plot
	scatter(cells.soma_mip2(:,1), cells.soma_mip2(:,2), somacolor, 'o');

	colormap(colorscheme);
	scatter3(flat.x, flat.y, flat.value, [], flat.value, markertype);
	%scatter3(flat.x, flat.y, flat.value, [], categorical(flat.value));
	hold off
	colorbar();
	title(cellquery)
	
	else % version 0, very slow

	%%{
	partners = get_cell_info(cell_info, flat0(1,:), [], true);
	%partners = struct2table(partners);
	partners = get_cell_info_table(partners);
	partners.type = categorical(partners.type);
	%length(find(partners.class<0))

	warning('bug here. contains <undefined>/NaN. see the other version')
	value = partners.class*100 + double(partners.type);

	%plot(flat0(2,:), flat0(3,:), '.');
	%scatter(flat0(2,:), flat0(3,:), 1, partners.class);
	scatter(flat0(2,:), flat0(3,:), [], value);
	%plot([],[]); hold on; % initiate view as 2D plot
	%scatter3(flat0(2,:), flat0(3,:), value, [], value);
	%}
	end

case {'self'}	% within type contacts
	%close(fig);
	%self(cell_info, cellquery, '--', cellquery, varargin{:});
	%return;
	pool = get_cell_info_table(cell_info, cellquery, varargin{:});

	bg = flat;	% hack to plot the full cells in the background
	bg.color = categorical(bg.cell1_id);

	flat = innerjoin(flat, pool);

	colormixing = 0;
	colors = distinguishable_colors(length(categories(bg.color)));
	colors0 = colors;
	%colors = colors + 0.5; colors(colors>1) = 1;
	colors = 0.5*colors + 0.5;
	if ~colormixing
		% No color mixing:
		colors = distinguishable_colors(1+length(categories(bg.color)));
		colors = colors([1:3 5:end], :);  % 4th is close to black

		colors = 0.5*colors + 0.5;
	else
		% Color mixing:
		flat.cell_id = categorical(flat.cell_id);
		flat.cell1_id = categorical(flat.cell1_id);	% these two sets should be identical
		flat.value = colors0(flat.cell1_id, :) + colors0(flat.cell_id, :);
		flat.value = flat.value./2;
		flat.value(flat.value>1) = 1;
	end
	
	colormap(colors);

	hold on
	scatter(bg.x, bg.y, [1], bg.color, '.');
	scatter(cells.soma_mip2(:,1), cells.soma_mip2(:,2), somacolor, 'h');
	
	if ~colormixing
		%scatter(flat.x, flat.y, [30], colors(end, :), 'o');  %markertype
		%scatter(flat.x, flat.y, [30], 'k', 'o');  %markertype
		scatter(flat.x, flat.y, [1], 'k', markertype);
	else
		scatter(flat.x, flat.y, [30], flat.value, 'o');  %markertype
	end
	hold off;

	title(cellquery)

	%xlim([0 5376]);ylim([0,3456])  % 5376,3456
	xlim([0 6000]);ylim([0,3500])


case {'bc', 'gc', 'gc ac', 'ac', '--'}

	if 0	% filter
		include = flat.z < ipl_to_mip2warpedi(0.5);
		flat = flat(include,:);
	end

	switch feature
	case '--'
		nvarargin = length(varargin);
		optargs = {{}, {}};
		optargs(1:nvarargin) = varargin;
		[types, vararg] = optargs{:};

		pool = get_cell_info_table(cell_info, types, vararg{:});
	case 'bc'
		pool = struct2table(get_cell_info(cell_info, {'bc', 'xbc', 'rbc'}, true)); % or use class 3
	case 'gc'
		pool = get_cell_info_table(cell_info, 1); % class 1
	case 'ac'
		pool = get_cell_info_table(cell_info, 2); % class 2
		pool.type(cellfun(@isempty,pool.type)) = {'unknown?'};
	case 'gc ac'
		pool = [get_cell_info(cell_info, 1); get_cell_info(cell_info, 2)];
		pool = struct2table(pool);
		pool.type(cellfun(@isempty,pool.type)) = {'unknown?'};
		%pool = struct2table(get_cell_info(cell_info, 2));
		%pool.soma_coord = [];
		%pool.is_cutoff = [];
		%pool2 = struct2table(get_cell_info(cell_info, 1));
		%pool = [pool2; pool];

		%pool = [pool; struct2table(get_cell_info(cell_info, 1))];
	end

	pool.type = categorical(pool.type);
	typesinpool = categories(pool.type(1));
	ntypes = length(typesinpool)
	pool = pool(:, {'cell_id', 'type'});

	bg = flat;	% hack to plot the full cells in the background
	flat = innerjoin(flat, pool);
	flat = flat(~isundefined(flat.type), :);
	flat.Properties.VariableNames{'type'} = 'value';

	scatter([],[]); hold on; % initiate view as 2D plot
	scatter(bg.x, bg.y, [1], 0.8*[1 1 1], '.');
	scatter(cells.soma_mip2(:,1), cells.soma_mip2(:,2), somacolor, 'h');

	colormap(distinguishable_colors(ntypes));
	caxis([1, ntypes+0.1]) % prevent rescaling of colormap
	%scatter3(flat.x, flat.y, flat.value, [1], flat.value, markertype);
	%scatter3(flat.x, flat.y, flat.value, [1], flat.value, 'filled', markertype);
	scatter(flat.x, flat.y, [1], flat.value, markertype);
	hold off;

	colorbar('Ticks',1:ntypes, 'TickLabels', typesinpool);
	title(cellquery)

	%xlim([0 5376]);ylim([0,3456])  % 5376,3456
	xlim([0 6000]);ylim([0,3500])

otherwise
	error('invalid input: feature');
end
%axis equal
daspect([10 10/23*16.5 1])

realnpoints = size(flat,1)

counts = grpstats(flat, {'cell1_id'}, {'numel'}, 'DataVars',{'cell1_id'});
counts.numel_cell1_id = []; % numel_cell1_id or GroupCount
contact_stat = join(counts, contact_stat);
if size(contact_stat,1) > 1
contact_stat{'Total', :} = [0 sum(contact_stat{:,2:end})];
end
% total contact : total surface
contact_stat.percent_contact = 100 * contact_stat.contact_area ./ double(contact_stat.surface_area); % originally int64
% contact of interest : total surface
contact_stat.percent_shown = 100 * contact_stat.GroupCount ./ double(contact_stat.surface_area); % originally int64
contact_stat.percent_of_contact = 100 * contact_stat.GroupCount ./ double(contact_stat.contact_area); % originally int64
if size(contact_stat,1) > 1
%summary(contact_stat)
tmp = contact_stat('Total', :);
contact_stat('Total', :) = [];
tmp{'Median', :} = [zeros(1,4) median(contact_stat{:,5:end})];	% excluding last row of 'Total'
tmp{'Mean', :} = [zeros(1,4) mean(contact_stat{:,5:end})];	% excluding last row of 'Total'
tmp{'Std', :} = [zeros(1,4) std(contact_stat{:,5:end})];	% excluding last row of 'Total'
% Doesn't work (switch row order for Mean and Total), row names are not transferred/assigned
%contact_stat(end-2:end, :) = contact_stat(end:-1:end-1, :);

contact_stat = vertcat(contact_stat, tmp(end:-1:1, :));
contact_stat(end-3:end, :)

ratio = hull_overlap(cell_info, cellquery);
contact_stat{'Total',end-2:end} / ratio
end
