
%{
function cell_stat = cell_info_get_type_property(cell_info, property_type_arg, use_normalized_strat, varargin)

if ~iscell(property_type_arg)
    property_type_arg = {property_type_arg};
end
property_type = property_type_arg{1};
property_arg = property_type_arg(2:end);
%}
%property_type = 'ca_sus_ind_on'; type_list = types_ca_has_on;
%property_type = 'ca_sus_ind_off'; type_list = types_ca_only_off;


%gc_types = list_types(cell_info);

get_super_types
%ca_groups = ca_groups2;
ca_groups = ca_groups4;
ca_groups = ca_groups_alphas;
n_groups = size(ca_groups,1);
%N = length(n_groups);

colors=distinguishable_colors(n_groups);

roi_sums_xcond_typemeans = table();
roi_sums_xcond_typemeans_rescale = table();
type_stat = table();


cellwise_value = 1; % cellwise sustainedness or cellwise risetime
cellwise_size_correction = 0;	% technically implies cellwise_value... TODO: check
compute_risetime = 0;
x_is_peak = 1;
filter_small_response = 3;
filter_small_response = -1;  % snr based, w/o size_correction
filter_small_response = -0.6;  % snr based (Qi), w/o size_correction
filter_small_response = -0.9;  % snr based
filter_small_response = -0.5;  % snr based (Qi)
%filter_small_response = 0;
%draw_indiviual = 1;   filter_small_response = 0;
size_correction = 1;
%size_correction = 0;
frametime_correction = 1;  % effective only when size_correction = 1. account for the fact that frame 15 and 23 are actually slightly different in term of timing
%frametime_correction = 0;
draw_indiviual = 0;
%draw_indiviual = 1;
anova_separate_onoff = 0;

vertlayout = 0
vary_marker_size = 1;


final_eps_figs = 1;
if final_eps_figs  % overrides
	if 0
	elseif 0    % ext fig
		size_correction = 1;
		vary_marker_size = 0;
	elseif 1	% main fig
		size_correction = 0;
	elseif 1  % typewise traces for ext fig. also remember to adjust the super groups
		size_correction = 0;
	end
end


figure;
legends = {};
legendentries = [];

anova_values = [];
anova_groups = [];

%{
fake_cell_info, fake_cell_dict_j, fake_roi_sums_means_flatten
filtered_cell_info, rescaled_roi_sums_means_flatten
or filtered_cell_info,  normalization_max / value_to_scale_to_1 but not on off...
ah could also have just modified the "cells"
%}
fake_cell_info = cell_info([]);
fake_cell_info(1).type = 'to be set';
fake_cell_info(1).cell_id = 1e6;
fake_cell_dict_j = [1 1e6];
fake_roi_sums_means_flatten = []; %sin(1:31);
filtered_cell_info = {[], []};  % on, off
cellwise_scaling = {nan(1,n_rois), nan(1,n_rois)};  % on, off

typecounts = table(cell(0,1), zeros(0,1), 'VariableNames', {'typename', 'count'});
all_values = [];
frameoffsets = [];
radii_all = [];
%ca_groups.
%for group = ca_groups.name
for k = 1:n_groups
    %group = strat_groups.types.'
	ca_groups.name{k}
    typelist = ca_groups.types{k};
    if ca_groups.has_on(k)
    	normalization = 'on';
    	property_type = 'ca_sus_ind_on';
    else
    	normalization = 'off';
    	property_type = 'ca_sus_ind_off';
    end

    group_values = [];
	for celltype = typelist(:).'
		celltype = celltype{1};

		[ca_ids, cells] = get_ca_ids(cell_dict_j, cell_info, celltype, false);
		cellwise_xcondmeans = roi_sums_xcondmeans(:, ca_ids);
		cellwise_raw = roi_sums_all_reshape(:, : ,:, ca_ids);
		roi_sums_xcond_typemeans{celltype, :} = mean(cellwise_xcondmeans, 2).';
		cellwise_min = min(cellwise_xcondmeans, [], 1);
		cellwise_rescale = cellwise_xcondmeans - repmat(cellwise_min, 31, 1);
		if 0 && filter_small_response
			cellwise_max = max(cellwise_rescale, [], 1);
			include = cellwise_max>filter_small_response;
			cellwise_rescale = cellwise_rescale(:, include);
			ca_ids = ca_ids(include);
			cells = cells(include);
		end
		roi_sums_xcond_typemeans_rescale{celltype, :} = roi_sums_xcond_typemeans{celltype, :} - min(roi_sums_xcond_typemeans{celltype, :});
		roi_sums_xcond_typemeans_rescale{celltype, :} = roi_sums_xcond_typemeans_rescale{celltype, :} / max(roi_sums_xcond_typemeans_rescale{celltype, :});

		if size_correction
			radii = cell_info_get_strat_property(cells, {'maxdiameter' 'um'}) / 2;
			if ~cellwise_size_correction
				radii = max(radii);  % taking the biggest cell for now
			end
			frameoffset = -1.6 + radii/(1000*0.128); %1000um/s * 0.128s/frame  % -1.6 is taken such that the max frameoffset is still negative
		end

		switch property_type
		case {'ca_sus_ind_on'}
			[maxi, maxind] = max(roi_sums_xcond_typemeans_rescale{celltype, 7:15});
			maxind = maxind + 6;
			%{
			while diff(roi_sums_xcond_typemeans_rescale{celltype, maxind:maxind+1})>0 ...
				 || diff(roi_sums_xcond_typemeans_rescale{celltype, maxind:2:maxind+2})>0
				disp(celltype)
				maxind = maxind+1
				maxi = roi_sums_xcond_typemeans_rescale{celltype, maxind};
			end
			frame = maxind + 3;
			%}
			frame = 15;
			if filter_small_response < 0 
				cellwise_raw = cellwise_raw(4:15, :, :, :);
				cellwise_raw = reshape(cellwise_raw, (15-4+1)*8, 4, []);
				% cellwise_raw = cellwise_raw(1:15, :, :, :);
				% cellwise_raw = reshape(cellwise_raw, (15-1+1)*8, 4, []);
				if filter_small_response <= -0.8  % snr
					snr = var(mean(cellwise_raw, 2), [], 1) ./ mean(var(cellwise_raw,[],2), 1);
				else % Qi. -1 < filter_small_response < 0
					snr = var(mean(cellwise_raw, 2), [], 1) ./ mean(var(cellwise_raw,[],1), 2);
				end
				%snr = reshape(var(reshape(cellwise_raw,[],7),[],1), 1,1,7) - mean(var(cellwise_raw,[],1), 2) - var(mean(cellwise_raw, 2), [], 1) - mean(var(cellwise_raw,[],2), 1)
				include = snr > -filter_small_response;
				cellwise_rescale = cellwise_rescale(:, include);
				cells = cells(include);
				ca_ids = ca_ids(include);
				if cellwise_size_correction
					frameoffset = frameoffset(include);
				end
			end
			if cellwise_value
				[~, maxind] = max(cellwise_rescale(7:15,:), [], 1);
				maxind = maxind + 6;
			end

			cellwise_rescale_copy = cellwise_rescale;
			% if cellwise_value || size_correction
			% 	frame = repmat(frame, size(cells));
			% end
			if 0 && size_correction
				cellwise_max = max(cellwise_rescale(7:14,:), [], 1);
			elseif 1 && size_correction
				frame = frame + frameoffset;
				if ~cellwise_size_correction
					frame = repmat(frame, size(cells));
				end
				frameceil = ceil(frame);
				%tmp = cellwise_rescale;
				
				cellwise_max = zeros(1, length(frame));
				for ii = 1:length(cells)
					cellwise_rescale(frameceil,ii) = interp1(cellwise_rescale(:,ii), frame(ii));
					cellwise_max(ii) = max(cellwise_rescale(7:frameceil,ii));
				end
			else
				cellwise_max = max(cellwise_rescale(7:15,:), [], 1);
			end
			%{
			cellwise_rescale = cellwise_rescale ./ repmat(cellwise_max, 31, 1);
			if filter_small_response > 0
				cellwise_rescale = cellwise_rescale(:, cellwise_max>filter_small_response);
			end
			typemean = mean(cellwise_rescale, 2);
			typestd = std(cellwise_rescale, [], 2); %std(cellwise_rescale, 0, 2);

			roi_sums_xcond_typemeans_rescale{celltype, :} = roi_sums_xcond_typemeans_rescale{celltype, :} / maxi;
			
			%stat_ca
			tmp = roi_sums_xcond_typemeans_rescale{celltype, :};
			tmp = typemean;
			%sus_idx = (1-tmp(13))*sign(tmp(13)-tmp(12));
			%sus_idx = (1-tmp(14))*sign(tmp(14)-tmp(13));
			%sus_idx = (1-tmp(15))*sign(tmp(15)-tmp(14));
			sus_idx = tmp(frame);
			%}

			%cells = get_ca_cell_info(cell_dict_j, cell_info, celltype);
		stat = cell_info_get_strat_property(cells, 'sus_on-trans_on');
		stat1 = cell_info_get_strat_property(cells, 'sus_on');
		stat2 = cell_info_get_strat_property(cells, 'trans_on');
		stat = (stat1-stat2) ./ (stat1+stat2);
		%stat = mean(stat1-stat2) ./ mean(stat1+stat2); % only slight visible difference for 27,28

			%yy2 = maxind - 1/0.128;
			yy2 = maxind*0.128 - 1;

		case {'ca_sus_ind_off'}
			% off
			[maxi, maxind] = max(roi_sums_xcond_typemeans_rescale{celltype, 15:end});
			maxind = maxind+14;
			%off
			frame = 24;
			frame = 23;

			if filter_small_response < 0 
				cellwise_raw = cellwise_raw(15:end, :, :, :);
				cellwise_raw = reshape(cellwise_raw, (31-15+1)*8, 4, []);
				% cellwise_raw = cellwise_raw(15:26, :, :, :);
				% cellwise_raw = reshape(cellwise_raw, (26-15+1)*8, 4, []);
				if filter_small_response <= -0.8  % snr
					snr = var(mean(cellwise_raw, 2), [], 1) ./ mean(var(cellwise_raw,[],2), 1);
				else % Qi. -1 < filter_small_response < 0
					snr = var(mean(cellwise_raw, 2), [], 1) ./ mean(var(cellwise_raw,[],1), 2);
				end
				include = snr > -filter_small_response;
				cellwise_rescale = cellwise_rescale(:, include);
				cells = cells(include);
				ca_ids = ca_ids(include);
				if cellwise_size_correction
					frameoffset = frameoffset(include);
				end
			end
			if cellwise_value
				[~, maxind] = max(cellwise_rescale(15:end,:), [], 1);
				maxind = maxind + 14;
			end

			cellwise_rescale_copy = cellwise_rescale;
			% if cellwise_value || size_correction
			% 	frame = repmat(frame, size(cells));
			% end
			if 1 && size_correction
				if frametime_correction
					frame = 15 + 1/0.128; % instead of 23
				end
				frame = frame + frameoffset;
				if ~cellwise_size_correction
					frame = repmat(frame, size(cells));
				end
				frameceil = ceil(frame);
				%tmp = cellwise_rescale;
				
				cellwise_max = zeros(1, length(frame));
				for ii = 1:length(cells)
					cellwise_rescale(frameceil,ii) = interp1(cellwise_rescale(:,ii), frame(ii));
					cellwise_max(ii) = max(cellwise_rescale(15:frameceil,ii));
				end
			else
				cellwise_max = max(cellwise_rescale(15:end,:), [], 1);
			end

			%cells = get_ca_cell_info(cell_dict_j, cell_info, celltype);
			stat1 = cell_info_get_strat_property(cells, 'sus_off');
			stat2 = cell_info_get_strat_property(cells, 'trans_off');
			stat = (stat1-stat2) ./ (stat1+stat2);
			%stat = cell_info_get_strat_property(cells, 'sus');

			%yy2 = maxind - 2/0.128;
			yy2 = maxind*0.128 - 2;

		end  %switch

		cellwise_rescale = cellwise_rescale ./ repmat(cellwise_max, 31, 1);
		cellwise_rescale_copy = cellwise_rescale_copy ./ repmat(cellwise_max, 31, 1);
		if filter_small_response > 0
			include = cellwise_max>filter_small_response;
			cellwise_rescale = cellwise_rescale(:, include);
			cellwise_rescale_copy = cellwise_rescale_copy(:, include);
			cells = cells(include);
			ca_ids = ca_ids(include);
			cellwise_max = cellwise_max(ca_ids);
			if size_correction
			frame = frame(include);
			end
			frameceil = frameceil(include);
		end

		% done with all/any filtering
    	on_or_off = 2-logical(ca_groups.has_on(k));
		filtered_cell_info{on_or_off} = [filtered_cell_info{on_or_off}; cells];
		cellwise_scaling{on_or_off}(ca_ids) = cellwise_max;


		%typemean = mean(cellwise_rescale, 2);
		%typestd = std(cellwise_rescale, [], 2); %std(cellwise_rescale, 0, 2);

		if ~size_correction
			roi_sums_xcond_typemeans_rescale{celltype, :} = roi_sums_xcond_typemeans_rescale{celltype, :} / maxi;

			%stat_ca
			if ~cellwise_value %old: ~draw_indiviual
				tmp = roi_sums_xcond_typemeans_rescale{celltype, :};
				%tmp = typemean;
				sus_idx = tmp(frame);
			else
				sus_idx = cellwise_rescale(frame,:).';
			end

		else % size_correction
			if 0
				radii = cell_info_get_strat_property(cells, {'maxdiameter' 'um'})/2;
				frameoffset = -1 + radii/(1000*0.128); %1000um/s * 0.128s/frame
				frame = frame + frameoffset;
				%max(radii)
				%max(frame)
				sus_idx = zeros(size(frame));
				for ii = 1:length(cells)
					sus_idx(ii) = interp1(cellwise_rescale(:,ii), frame(ii));
				end
			else
				sus_idx = zeros(size(frame));
				for ii = 1:length(cells)
					sus_idx(ii) = cellwise_rescale(frameceil(ii),ii);
				end
			end

			frameoffsets = [frameoffsets; frameoffset repmat(ca_groups.has_on(k),size(frameoffset))];
			radii_all = [radii_all; radii];
		end

		sus_idx_mean = mean(sus_idx);
		sus_idx_std = std(sus_idx);


		%{
		if tmp(frame) > tmp(frame-1)
			warning('aaaaaa')
			celltype
			maxind
			tmp(frame-1:frame+1)
		end
		%}

		if x_is_peak
			stat = cell_info_get_strat_property(cells, {'peak', normalization, 'smooth'});
			if celltype(1) == '5'	% hack. Always use global peak for 5xx types
				stat = cell_info_get_strat_property(cells, {'peak', [], 'smooth'});
			end
			stat = stat/100; % 0 to 1
		end

		if ~draw_indiviual && 0
		%stat = median(stat);
		stat = mean(stat);
		end

		xx = stat;
		yy = sus_idx;
		if compute_risetime
			yy = yy2.';
		end
		typecount = length(yy);
		typecounts(end+1,:) = {{celltype}, typecount};
		if strcmp(celltype, '63') || strncmp(celltype, '2a', 2)
			disp(celltype)
			disp(typecount)
		end
		yy_mean = sus_idx_mean;
		yy_std = sus_idx_std;
		%size(yy)
		yy_mean = mean(yy);
		yy_std = std(yy);
		xx_mean = mean(xx);
		anova_values = [anova_values; yy];
		if anova_separate_onoff
			groupname = [ca_groups.name{k} ' (' normalization ')'];
			anova_groups = [anova_groups; repmat({groupname},size(yy))];
		else
			anova_groups = [anova_groups; repmat(ca_groups.name(k),size(yy))];
		end

		if ~draw_indiviual && 0
			yy = yy_mean;
		group_values(end+1,:) = [xx yy];
		else
			%yy = cellwise_rescale(frame,:);
			group_values(end+1:end+length(xx),:) = [xx yy(:)];
		end

		%type_stat{celltype, {'sus_idx', 'sus_on'}} = [sus_idx, stat];

		%k = size(type_stat, 1);
		%k = size(roi_sums_xcond_typemeans_rescale, 1);
		kolor = colors(mod(k,n_groups)+1,:);
		%h = scatter(xx, yy, 'filled','MarkerFaceColor',kolor,'MarkerEdgeColor',kolor);
		%h = scatter(xx, yy, 90, ca_groups.style1{k}{:}, 'filled');
		%h = scatter(xx, yy, 90, ca_groups.style2{k}{:}, 'filled');

		if draw_indiviual
			%size(xx)
			%size(yy)
			scatter(xx, yy, 50, ca_groups.style3{k}{:});
			%scatter(xx, yy);
		else
				if any(strcmp(ca_groups.style3{k}, 'filled'))
					border = {'MarkerEdgeColor', [1 1 1]};
				end
			if vary_marker_size
				h = scatter(xx_mean, yy_mean, 30*typecount, ca_groups.style3{k}{:}, border{:});
			else
		h = scatter(xx_mean, yy_mean, 90, ca_groups.style3{k}{:}, border{:});
			end

			if 1 %&& typecount<10 % draw_dots
				hold on
				scatter(xx_mean*yy./yy, yy, 7, ca_groups.style3{k}{:});
			end

			%errorbar([xx], yy, [yy_std])
			line([xx_mean,xx_mean], [yy_mean-yy_std yy_mean+yy_std], 'Color', ca_groups.style3{k}{end});

			fontsize = 13;
			if 1 && vary_marker_size
				textoffset = sqrt(typecount) * 0.009;
			else
				textoffset = 0.0105;
			end
			if x_is_peak
				textoffset = textoffset * 0.9/2;
			end
			if vertlayout
				text(xx_mean+textoffset,yy_mean, typename2displayname(celltype),  'VerticalAlignment', 'cap', 'FontSize', fontsize); %'top', 'HorizontalAlignment', 'center',
			else
				text(xx_mean+textoffset,yy_mean, typename2displayname(celltype),'FontSize', fontsize);
			end
			hold on
		end
		
		%scatter(repmat(sus_idx,size(stat)),stat, 'filled','MarkerFaceColor',colors(k,:),'MarkerEdgeColor',colors(k,:));
		hold on;
	end

	all_values = [all_values; group_values];
	[R, P] = corrcoef(group_values(:,1), group_values(:,2));
	if length(R)>1
		disp([R(2), P(2)])
	else
		disp('   not enough points to compute correlation')
	end
	legendentries(end+1) = h;
	legends{end+1} = [ca_groups.name{k} ' (' normalization ')'];
end %for k = 1:n_groups
if ~draw_indiviual
%legend(legendentries, ca_groups.name);
legend(legendentries, legends);
end
if x_is_peak
	if vertlayout
		view(90,90);
		figure_size_x2([1 2]);
	else
		figure_size_x2([2 1]);
	end
else
plot([0 0], [-1 1], '--k');
figure_size_x2();
end
set(gca, 'FontSize', 13)

[R, P] = corrcoef(all_values(:,1), all_values(:,2))
[~, tmp] = max(typecounts.count);
%typecounts(tmp,:)
largestdottext = [',  max(n) =' table_row_to_text(typecounts(tmp,:))]

if ~compute_risetime
xlabel 'Ca'
ylabel 'strat '
%xlim([-0.8 0])
ylim([-0.8 0]+1)

switch property_type
case {'ca_sus_ind_on'}
	xlabel 'Ca @ 15'
	ylabel(sprintf('strat\n (sus_on-trans_on)/(sus_on+trans_on)'), 'Interpreter', 'none')
	title('on sustainedness')
case {'ca_sus_ind_off'}
	xlabel 'Ca @ 23'
	ylabel(sprintf('strat\n (sus_off-trans_off)/(sus_off+trans_off)'), 'Interpreter', 'none')
	title('off sustainedness')
end %switch

ylabel 'Ca @ 15/23'
ylabel 'Sustainedness index'
xlabel(sprintf('strat\n (sus_on-trans_on)/(sus_on+trans_on)\n or \n (sus_off-trans_off)/(sus_off+trans_off)'), 'Interpreter', 'none')
title(['sustainedness' largestdottext])
else % compute_risetime
	xlabel(sprintf('strat\n (sus_on-trans_on)/(sus_on+trans_on)\n or \n (sus_off-trans_off)/(sus_off+trans_off)'), 'Interpreter', 'none')
	ylabel('rise time in #frames')
	title('rise time in #frames')
	ylabel('rise time in seconds')
	title('rise time')
end

if x_is_peak
	%xlabel('on/off peak strat location')
	xlabel('IPL depth of peak stratification (inner or outer)')
end
if size_correction && ~compute_risetime
	title('sustainedness adjusted by cell size, with caveat')
	ylabel('Ca @  15or23 -1 + radii/(1000*0.128)')
	%ylim([0 1.3])
	if cellwise_size_correction
		title(['sustainedness adjusted by cell size' largestdottext])
	else
		title(['sustainedness adjusted by cell size of type' largestdottext])
	end
	ylabel('Ca @  15or23 - 1.6 + radii/(1000*0.128)')
	ylabel('Sustainedness index')
end
if 0 && compute_risetime
	%title('risetime in #frames')
	title('risetime in seconds')
end

asdfsadfsa
%legend(type_stat.Properties.RowNames)
%legend(roi_sums_xcond_typemeans_rescale.Properties.RowNames)
hold off;
%scatter(type_stat{:,1},type_stat{:,2}, 'filled','MarkerFaceColor',colors(k,:),'MarkerEdgeColor',colors(k,:));


if draw_indiviual
	figure
	h = histogram(all_values(:,1),50);
	%h.fewerbins()
end

[p, t, stats] = anova1(anova_values, anova_groups);
drawnow; set(gca, 'XTickLabelRotation', 70)
%xtickangle(90) % r2016b

figure;[c,m,h,nms] = multcompare(stats, 'Alpha', 0.005);  %  p <= 0.0024

anova_data = table();
anova_data.value = anova_values;
anova_data.group = categorical(anova_groups, unique(anova_groups, 'stable'));

% grpstats(anova_data, )
% 'mean'	Mean
% 'sem'	Standard error of the mean
% 'numel'	Count, or number, of non-NaN elements
drawnow
figure; %grpstats(anova_values, anova_groups, 0.01);
grpstats(anova_values, anova_groups, 0.01);
if vertlayout
	view(90,90);
	figure_size_x2([1 2]);
else
	%drawnow;
	set(gca, 'XTickLabelRotation', 70)
end

%grpstats(anova_values, anova_groups, {'mean', 'std', 'sem', 'meanci'})
st = grpstats(anova_data, 'group', {'mean', 'std', 'sem', 'meanci'}, 'Alpha', 0.01);
st.mean_plus_std = [st.mean_value - st.std_value, st.mean_value + st.std_value]
figure
% plot(st.mean_value, st.group, 'o')
% hold on
% plot(st.mean_plus_std.', repmat(st.group.',2,1))
% axis ij  % set(gca,'Ydir','reverse')
% set(gca, 'YLim', [0 6], 'YTick', 1:5, 'YTickLabel', categories(anova_data.group))
errorbar(st.mean_value, st.std_value, 'x');
if vertlayout
	view(90,90)	%, 'Xdir', 'reverse')
	figure_size_x2([1 2]);
else
	set(gca, 'XTickLabelRotation', 70)
end
set(gca, 'XLim', [0 6], 'XTick', 1:5, 'XTickLabel', categories(anova_data.group))
if compute_risetime
	%title('risetime: mean +/- std, in #frames')
	title('risetime: mean +/- std')
else
	title('sustainedness: mean +/- std')
end

figure;
if compute_risetime
	jitter = 0.5;
else
	jitter = 0;
end
cats = 1:size(st,1);
if 0
	boxplot(anova_values, anova_groups);
	hold on
	errorbar([1:5]-0.4, st.mean_value, st.sem_value, 'x') %, 'MarkerSize', 20);
elseif 0
	boxplot(anova_values, anova_groups, 'Whisker', 0, 'Colors', 0.8*[1 1 1], 'Symbol', '.', 'BoxStyle', 'filled', 'Widths', 0.2); %, 'Orientation', 'horizontal'); 
	hold on
	errorbar([1:5]-0.2, st.mean_value, st.sem_value, 'x') %, 'MarkerSize', 20);
elseif 1
	%boxplot(anova_values, anova_groups, 'Whisker', 0, 'Colors', 0.8*[1 1 1], 'Symbol', ''); %, 'Orientation', 'horizontal');
	boxplot(anova_values, anova_groups, 'Colors', 0.8*[1 1 1], 'Symbol', '.', 'Jitter', jitter); %, 'Orientation', 'horizontal');
	hold on
	errorbar(cats, st.mean_value, st.sem_value, 'x') %, 'MarkerSize', 20);
	%errorbar(cats, st.mean_value, st.sem_value, '.', 'MarkerSize', 20, 'CapSize', 20)
elseif 0
	co = get(groot, 'defaultAxesColorOrder');  % blue, red, yellow/orange, purple, green
	co = co + (1-co)/2;
	%co = 1 - (1-co)/3;
	boxplot(anova_values, anova_groups, 'Colors', co, 'Symbol', '.', 'Jitter', jitter); %, 'Orientation', 'horizontal');
	hold on
	errorbar([1:5], st.mean_value, st.sem_value, 'xk'); %, 'CapSize', 7); %, 'LineWidth', 1) %, 'MarkerSize', 20);
end
if vertlayout
	view(90,90)	%, 'Xdir', 'reverse')
	%figure_size_x2([1 2]);
else
	set(gca, 'XTickLabelRotation', 70)
end
set(gca, 'XLim', [0 6], 'XTick', 1:5, 'XTickLabel', categories(anova_data.group))
ylabel('Sustainedness index')
if compute_risetime
	%title('risetime: mean +/- sem, in #frames')
	title('risetime: mean +/- sem')
else
	title('sustainedness: mean +/- sem')
end



n_groups = size(ca_groups,1);
figure;
for k = 1:n_groups
    %group = strat_groups.types.'
    subplot(n_groups, 1, k)
    typelist = ca_groups.types{k};
    on_or_off = 2-logical(ca_groups.has_on(k));
    if ca_groups.has_on(k)
    	normalization = 'on';
    else
    	normalization = 'off';
    end
    if 1 && cellwise_value
    	plot_grouped_ca(filtered_cell_info{on_or_off}, cell_dict_j, roi_sums_means_flatten, typelist, [], 0, normalization, 2.2, 1, typelist, cellwise_scaling{on_or_off});
    else
    plot_grouped_ca(cell_info, cell_dict_j, roi_sums_means_flatten, typelist, [], 0, normalization, 2.2);
	end
    ylabel(ca_groups.name{k}, 'Rotation', 0, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    l = legend();
    %legend('show', 'Location', 'eastoutside');
    %l.Position(3) = 0.15;
    if ~final_eps_figs
    	l.Location = 'eastoutside';
	end
    ax = gca();
    xtext = ax.XLabel;
    ax.XLabel = [];
    xticklabel = ax.XTickLabel;
    ax.XTickLabel = [];
end
ax.XTickLabel = xticklabel;
ax.XLabel = xtext;
return

ca_classical_groups = {
	types_oodsgc 'on'
	types_ondsgc 'on'
	types_alpha 'full'	%
	types_alpha_all 'full'	%
};
for k = 1:size(ca_classical_groups,1)
    [typelist, normalization] = ca_classical_groups{k,:};
    figure; plot_grouped_ca(cell_info, cell_dict_j, roi_sums_means_flatten, typelist, [], 0, normalization);
end
