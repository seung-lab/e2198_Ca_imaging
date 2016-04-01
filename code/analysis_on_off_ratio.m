

% 20239


show_on_offs = 0;
show_taus = 1;


on_off_ratios = [];
on_off_ratios_bytype = containers.Map();
sustain_indices = [];
sustain_indices_bytype = containers.Map();

excluded_cell_ids = [];

%% Ca ratio

coeffs16_reshape = reshape(coeffs16{1,2}(:,3:end-16).', 2, 8, 634);
coeffs16_reshape_sum = sum(coeffs16_reshape, 2);
coeffs16_reshape_max = max(coeffs16_reshape, [], 2);
coeffs16_reshape_sum = squeeze(coeffs16_reshape_sum);
coeffs16_reshape_max = squeeze(coeffs16_reshape_max);
ca_on_off_ratio = coeffs16_reshape_sum(1,:)./coeffs16_reshape_sum(2,:);
%ca_on_off_ratio = coeffs16_reshape_max(1,:)./coeffs16_reshape_max(2,:);
ca_taus = coeffs16{1,2}(:,1);

%% preferred directions
load('oodsgc_pref_axes.mat', 'oodsgc_axes');
% '''
% columns 2 and 3 represent the center of mass of contactivity. Columns 4 and 5 are
% the max of the off-layer connectivity over direction and Columns 6 and 7 are the same for ON
% they should all agree pretty well
% ''' - Matt
% Matt vs Kevin:
% x-   =  90
% x+   =  270
% y+   =  180
oodsgc_pref_dir_group = ...
 (oodsgc_axes(:,2) > abs(oodsgc_axes(:,3))) * 1 + ...
 (oodsgc_axes(:,2) < -abs(oodsgc_axes(:,3))) * 2 + ...
 (oodsgc_axes(:,3) >= abs(oodsgc_axes(:,2))) * 3 + ...
 (oodsgc_axes(:,3) <= -abs(oodsgc_axes(:,2))) * 4;
oodsgc_pref_dir_group = repmat(oodsgc_pref_dir_group, [1 3]);

%%{
oodsgc_kevin_x = -oodsgc_axes(:,[3 5 7]);
oodsgc_kevin_y = -oodsgc_axes(:,[2 4 6]);
oodsgc_pref_dir = atan2(oodsgc_kevin_y, oodsgc_kevin_x);
oodsgc_pref_dir_group = mod( round((oodsgc_pref_dir + 2*pi) / (pi/4)) - 1, 8 ) + 1;
oodsgc_pref_dir_column = order(oodsgc_pref_dir_group);
%}


%% contact ratio
load('ooDSGC_vericose_contact_summary.mat', 'dsgc_vericose_contact_summary');
% 'columns: 1- ooDSGC ID, 2- OFF contact, 3- ON contact' - Matt
contact_on_off_ratio = dsgc_vericose_contact_summary(:,3) ./ dsgc_vericose_contact_summary(:,2);


%% strat ratio
%cells = get_cell_info(cell_info, 20239); %'37');
cells = get_cell_info(cell_info, '37');
cells = cell_info;
%{
cells = [cells; get_cell_info(cell_info, '7i')];
cells = [cells; get_cell_info(cell_info, '7o')];
cells = [cells; get_cell_info(cell_info, '1wt')];
cells = [cells; get_cell_info(cell_info, '4ow')];
cells = [cells; get_cell_info(cell_info, '6sw')];
cells = [cells; get_cell_info(cell_info, '8w')];
%}
filter_low_confidence = 0;
for cell_info_elem = cells.'

	%strat = cell_info_elem.strat_nrml;
	strat = cell_info_elem.strat_unrml;
	binsums = strat(:,2);

	% 45% is eyeballed. 47% or 48% may be better, see 5ti(?)
	off_vol = sum(binsums( strat(:,1)>10 & strat(:,1)<45 ));
	on_vol = sum(binsums( strat(:,1)<80 & strat(:,1)>45 ));
	%off_vol = sum(binsums( strat(:,1)>0 & strat(:,1)<45 ));
	%on_vol = sum(binsums( strat(:,1)<90 & strat(:,1)>45 ));
	strat_ratio = on_vol/off_vol;


	sustained_vol = sum(binsums( strat(:,1)>10 & strat(:,1)<28 )) + sum(binsums( strat(:,1)>62 & strat(:,1)<80 ));
	transient_vol = sum(binsums( strat(:,1)>28 & strat(:,1)<62 ));
	%{
	sustained_vol = sum(binsums( strat(:,1)>0 & strat(:,1)<20 )) + sum(binsums( strat(:,1)>70 & strat(:,1)<80 ));
	transient_vol = sum(binsums( strat(:,1)>36 & strat(:,1)<54 ));
	%}
	strat_sustain_index = sustained_vol / transient_vol;
	%strat_sustain_index = log(1 / transient_vol);
	%strat_sustain_index = log(sustained_vol / transient_vol);
	%transience_index

	%{
	%TODO: 37 is fine. but need mod for others where the tow largest peaks don't necessary uniquely match the on off peaks.
	binsums = smooth(binsums);
	[pks,locs, w] = findpeaks(binsums, 'SortStr', 'descend');

	%for ii locs(1) < locs(2)
	%area = sum(binsums()
	% estimated area
	pks = pks.*w;

	% bins are actually in descending order in our data
	if locs(1) < locs(2)
		pks = pks(2:-1:1);
	else
		pks = pks(1:2);
	end
	strat_ratio = pks(2)/pks(1);
	%}

	omni_id = cell_info_elem.cell_id;
	ca_id = cell_dict(cell_dict(:,1)==omni_id, 2);

	% skip outliers
	if omni_id == 20239 && show_on_offs
		excluded_cell_ids(end+1) = omni_id;
		continue;
	end

	if isempty(ca_id) || length(ca_id)>1

		% directly skip for now to simplify computing the correlation
		excluded_cell_ids(end+1) = omni_id;
		continue;

		ca_id = NaN;
		ca_ratio = NaN;
	else
		ca_ratio = ca_on_off_ratio(ca_id);

		%%{
		filter_aggressively = 0;
		if filter_aggressively
			sum_min = 25;
			max_min = 7;
			if coeffs16_reshape_sum(1,ca_id)<sum_min && coeffs16_reshape_sum(2,ca_id)<sum_min
				ca_id
				excluded_cell_ids(end+1) = omni_id;
				continue
			end
			if coeffs16_reshape_max(1,ca_id)<max_min && coeffs16_reshape_max(2,ca_id)<max_min
				ca_id
				%coeffs16_reshape_max(:,ca_id)
				excluded_cell_ids(end+1) = omni_id;
				continue
			end
		else
			sum_min = 15;
			max_min = 4;
		end
		if coeffs16_reshape_sum(1,ca_id)<15 && coeffs16_reshape_sum(2,ca_id)<15
			ca_id
			excluded_cell_ids(end+1) = omni_id;
			continue
		end
		if coeffs16_reshape_max(1,ca_id)<4 && coeffs16_reshape_max(2,ca_id)<4
			ca_id
			%coeffs16_reshape_max(:,ca_id)
			excluded_cell_ids(end+1) = omni_id;
			continue
		end
		%}
	end

	%dir_group = oodsgc_pref_dir_group( oodsgc_axes(:,1)==omni_id );
	dir_group = oodsgc_pref_dir_group( oodsgc_axes(:,1)==omni_id, :) / 8;
	if isempty(dir_group)
		%dir_group = 0;
		dir_group = [0 0 0];
	end
	if omni_id == 26103 % mismatch between Matt and Kevin
		%excluded_cell_ids(end+1) = omni_id;
		%continue;
		%dir_group = 5;
		dir_group(1) = -dir_group(1); % opposite
	end
	if omni_id == 26032 && filter_low_confidence % mismatch between Matt and Kevin
		excluded_cell_ids(end+1) = omni_id;
		continue;
		%dir_group = 5;
	end

	%{
	% ca_ratio on the prefered direction only
	column = oodsgc_pref_dir_column( oodsgc_axes(:,1)==omni_id, : );
	%[omni_id coeffs16_reshape(1, column, ca_id) coeffs16_reshape(2, column, ca_id)]
	% one dir
	ca_ratio = coeffs16_reshape(1, column(1), ca_id) ./ coeffs16_reshape(2, column(1), ca_id);
	% separate on/off dir
	ca_ratio = coeffs16_reshape(1, column(3), ca_id) ./ coeffs16_reshape(2, column(2), ca_id);
	if coeffs16_reshape(1, column(3), ca_id) < 5 && filter_low_confidence
		omni_id
		excluded_cell_ids(end+1) = omni_id;
		continue;
	end
	%}

	contact_ratio = contact_on_off_ratio(dsgc_vericose_contact_summary(:,1)==omni_id);
	if isempty(contact_ratio)
		contact_ratio = 0;
	end

	% skip fitting problems or division by 0s (no off ca/strat)
	if show_on_offs
		if ca_ratio > 100 || isnan(ca_ratio)
			excluded_cell_ids(end+1) = omni_id;
			continue;
		end
		if isnan(strat_ratio)
			omni_id
			[on_vol off_vol]
			excluded_cell_ids(end+1) = omni_id;
			continue;
		end
		if isinf(strat_ratio)
			strat_ratio = 100;
			excluded_cell_ids(end+1) = omni_id;
			continue;
		end
	end

	if cell_info_elem.class ~= 1   % not Ganglion cell
		excluded_cell_ids(end+1) = omni_id;
		continue;
	end

	if show_taus && isinf(strat_sustain_index)
		strat_sustain_index = 1e4;
		excluded_cell_ids(end+1) = omni_id;
		continue;
	end
	if show_taus
		if 0 && omni_id == 20075	
			excluded_cell_ids(end+1) = omni_id;
			continue;
		end
	end

	% skip fitting problems.
	if ca_taus(ca_id)<1		% e.g. 26070
		excluded_cell_ids(end+1) = omni_id;
		continue;
	end

	sustain_indices(end+1, :) = [omni_id ca_id strat_sustain_index contact_ratio ca_taus(ca_id)];


	%on_off_ratios(end+1, :) = [omni_id ca_id strat_ratio ca_ratio dir_group];
	on_off_ratios(end+1, :) = [omni_id ca_id strat_ratio contact_ratio ca_ratio dir_group];
	ctype = cell_info_elem.type;
	if ~ischar(ctype)
		ctype = '';
		excluded_cell_ids(end+1) = omni_id;
		continue;  % skip unidentified cell
	end
	if on_off_ratios_bytype.isKey(ctype)
		tmp = on_off_ratios_bytype(ctype);
		tmp2 = sustain_indices_bytype(ctype);
	else
		tmp = [];
		tmp2 = [];
	end
	tmp(end+1, :) = on_off_ratios(end, :);
	tmp2(end+1, :) = sustain_indices(end, :);
	on_off_ratios_bytype(ctype) = tmp;
	sustain_indices_bytype(ctype) = tmp2;

end

on_off_ratios_type = [];
sustain_indices_type = [];
for ctype = on_off_ratios_bytype.keys()
	ctype = ctype{1};
	on_off_ratios_type(end+1,:) = mean(on_off_ratios_bytype(ctype), 1);
	sustain_indices_type(end+1,:) = mean(sustain_indices_bytype(ctype), 1);
	tmp = sustain_indices_bytype(ctype);
	sustain_indices_type(end,4) = std(tmp(:,5));
	if strcmp(ctype, '8w') || strcmp(ctype, '9w')
		sustain_indices_bytype(ctype)
	end
end

if show_on_offs
	[R, P] = corrcoef(on_off_ratios(:,3), on_off_ratios(:,5))
	%[R, P] = corrcoef(on_off_ratios(:,4), on_off_ratios(:,5))
	[R, P] = corrcoef(log(on_off_ratios(:,3)), log(on_off_ratios(:,5)))
	figure;
	%plot(on_off_ratios(:,3), on_off_ratios(:,4), '.');
	scatter(on_off_ratios(:,3), on_off_ratios(:,5), [], on_off_ratios(:,6), 'filled');
	%scatter(on_off_ratios(:,3), on_off_ratios(:,5), [], on_off_ratios(:,6:end), 'filled');
	%loglog(on_off_ratios(:,3), on_off_ratios(:,5), '.');
	%text(on_off_ratios(:,3), on_off_ratios(:,5), num2str(on_off_ratios(:,1:2), '%d-%d'));
	%axis equal;
	%loglog(on_off_ratios_type(:,3), on_off_ratios_type(:,5), '.');
	[R, P] = corrcoef(log(on_off_ratios_type(:,3)), log(on_off_ratios_type(:,5)))
	%text(on_off_ratios_type(:,3), on_off_ratios_type(:,5), char(on_off_ratios_bytype.keys()));
	xlabel('strat'); ylabel('Ca');

	%{
	figure;
	scatter(on_off_ratios(:,4), on_off_ratios(:,5), [], on_off_ratios(:,6), 'filled');
	axis equal;
	xlabel('contact'); ylabel('Ca');
	%}
end

if show_taus
	[R, P] = corrcoef(sustain_indices(:,3), sustain_indices(:,5))
	%[R, P] = corrcoef(log(sustain_indices(:,3)), log(sustain_indices(:,5)))
	figure;
	scatter(sustain_indices(:,3), sustain_indices(:,5), [], on_off_ratios(:,6), 'filled');

	%loglog(sustain_indices_type(:,3), sustain_indices_type(:,5), '.');
	semilogx(sustain_indices_type(:,3), sustain_indices_type(:,5), '.');
	[R, P] = corrcoef(sustain_indices_type(:,3), sustain_indices_type(:,5))
	text(sustain_indices_type(:,3), sustain_indices_type(:,5), char(sustain_indices_bytype.keys()));
	colors = ismember(sustain_indices_bytype.keys(), {'6sw', '8w', '4ow', '1wt'});
	%colors = [colors.' zeros(length(colors), 2)];
	names = char(sustain_indices_bytype.keys());
	text(sustain_indices_type(colors,3), sustain_indices_type(colors,5), names(colors,:), 'Color', 'red');
	xlabel('strat'); ylabel('tau');
end



% Other:
%26158 is 37?