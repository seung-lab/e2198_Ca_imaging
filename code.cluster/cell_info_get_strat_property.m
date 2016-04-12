function cell_stat = cell_info_get_strat_property(cell_info, property_type, use_normalized_strat, varargin)
% Note the cell_info input here is expected to be the already filtered subset of interest, rather than the full array of all cells.
%   	# TODO: perhaps I should just add the celltype argument and do filtering here.
% varargin: used for certain properties, such as 'corr' against a certain type of known strat profile.

N = length(cell_info);
cell_stat=zeros(N, 1);

if ~exist('use_normalized_strat', 'var')
	use_normalized_strat = true;
end

%strat=cell_info_bin_strat(cell_info,binstep);

if use_normalized_strat
	stratname = 'strat_nrml';
else
	stratname = 'strat_unrml';
end


%TODO: this for loop can be made into an array operation instead

for k = 1:N
	if use_normalized_strat
		strat = cell_info(k).strat_nrml;
    	binwidth = abs(strat(2,1) - strat(1,1));
    	strat(:,2) = strat(:,2) * binwidth;
	else
		strat = cell_info(k).strat_unrml;
	end
	s=strat(:,2);
    x=strat(:,1);

	switch property_type
	case {'sus_strat_vol', 'cum_sus', 'sus'}
        cell_stat(k) = cell_info_get_strat_property(cell_info(k), 'sus_on') + cell_info_get_strat_property(cell_info(k), 'sus_off');

	case {'trans_strat_vol', 'cum_trans', 'trans'}
        cell_stat(k) = cell_info_get_strat_property(cell_info(k), 'trans_on') + cell_info_get_strat_property(cell_info(k), 'trans_off');

	case {'cum_on', 'on'}
        cell_stat(k) = sum(s(x>45));

	case {'cum_off', 'off'}
        cell_stat(k) = sum(s(x<45));

	case {'cum_trans_on', 'trans_on'}
        cell_stat(k) = sum(s(x>45 & x<62));

	case {'cum_sus_on', 'sus_on'}
        cell_stat(k) = sum(s(x>62));

	case {'cum_trans_off', 'trans_off'}
        cell_stat(k) = sum(s(x<45 & x>28));

	case {'cum_trans_on', 'sus_off'}
        cell_stat(k) = sum(s(x<28));

   	case {'on-off'}
        cell_stat(k) = cell_info_get_strat_property(cell_info(k), 'on') - cell_info_get_strat_property(cell_info(k), 'off');

   	case {'sus-trans'}
        cell_stat(k) = cell_info_get_strat_property(cell_info(k), 'sus') - cell_info_get_strat_property(cell_info(k), 'trans');

   	case {'sus_on-trans_on'}
        cell_stat(k) = cell_info_get_strat_property(cell_info(k), 'sus_on') - cell_info_get_strat_property(cell_info(k), 'trans_on');

   	case {'sus_off-trans_off'}
        cell_stat(k) = cell_info_get_strat_property(cell_info(k), 'sus_off') - cell_info_get_strat_property(cell_info(k), 'trans_off');

   	case {'sus_on-sus_off'}
        cell_stat(k) = cell_info_get_strat_property(cell_info(k), 'sus_on') - cell_info_get_strat_property(cell_info(k), 'sus_off');

	case {'corr', 'sac_corr'}
		sac = varargin{1};
		cell_stat(k) = s.'*sac / sqrt(sum(s.^2) * sum(sac.^2));

	case {'corr_unrml'}
		sac = varargin{1};
		cell_stat(k) = s.'*sac;% / sqrt(sum(s.^2) * sum(sac.^2));

    % case {'sus/tran'}
    % 	cell_stat(k) = cell_info_get_strat_property(cell_info())
    %     if at
    %         cell_stat(j) = log(as/at);
    %     else
    %         cell_stat(j) = 20;
    %     end

    otherwise
        error('not recognized stat name') 

    end
end
