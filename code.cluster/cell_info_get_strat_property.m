function cell_stat = cell_info_get_strat_property(cell_info, property_type_arg, use_normalized_strat, varargin)
% Note the cell_info input here is expected to be the already filtered subset of interest, rather than the full array of all cells.
%   	# TODO: perhaps I should just add the celltype argument and do filtering here.
% varargin: used for certain properties, such as 'corr' against a certain type of known strat profile.

self = @cell_info_get_strat_property;

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

if ~iscell(property_type_arg)
    property_type_arg = {property_type_arg};
end
property_type = property_type_arg{1};
property_arg = property_type_arg(2:end);

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

    cell_info_elem = cell_info(k);

    onoff_boundary = 47;
    
    switch property_type
    case {'peak'}  % WARNING: potentially incorrect for cells saddling across the ON/OFF boundary
        if ~isempty(property_arg) && ~isempty(property_arg{1})
            onoff = property_arg{1};
        else
            onoff = '';
        end
        if ~isempty(property_arg) && ~isempty(property_arg{2})
            s = smooth(s); % smoothing
        end
        switch onoff
        case {'', 'both'}
            % nothing
        case {'on'}
            s = s(x>=onoff_boundary);
            x = x(x>=onoff_boundary);
        case {'off'}
            s = s(x<=onoff_boundary);
            x = x(x<=onoff_boundary);
        otherwise
            error('invalid argument: unrecognized onoff value')
        end
        [~, ind] = max(s);
        cell_stat(k) = x(ind);

	case {'sus_strat_vol', 'cum_sus', 'sus'}
        cell_stat(k) = cell_info_get_strat_property(cell_info(k), 'sus_on') + cell_info_get_strat_property(cell_info(k), 'sus_off');

	case {'trans_strat_vol', 'cum_trans', 'trans'}
        cell_stat(k) = cell_info_get_strat_property(cell_info(k), 'trans_on') + cell_info_get_strat_property(cell_info(k), 'trans_off');
        % cell_stat(k) = cell_stat(k) + sum(s(x==28 | x==62));

	case {'cum_on', 'on'}
        cell_stat(k) = sum(s(x>onoff_boundary));

	case {'cum_off', 'off'}
        cell_stat(k) = sum(s(x<onoff_boundary));

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

    case {'trans_on-trans_off'}
        cell_stat(k) = cell_info_get_strat_property(cell_info(k), 'trans_on') - cell_info_get_strat_property(cell_info(k), 'trans_off');

    case {'sus_off-trans_on'}
        cell_stat(k) = cell_info_get_strat_property(cell_info(k), 'sus_off') - cell_info_get_strat_property(cell_info(k), 'trans_on');

    case {'stdev'}
        cell_stat(k) = sqrt( (x.^2).'*s./sum(s) - (x.'*s./sum(s)).^2 );

    case {'stdev2' 'moment2'}  % second moment about a given reference point
        cell_stat(k) = sqrt( ((x-property_arg{1}).^2).'*s./sum(s) );

    case {'<'}  % amount of strat less than given value
        cell_stat(k) = sum(s(x<property_arg{1}));

    case {'><'}
        cell_stat(k) = sum(s(x>property_arg{1} & x<property_arg{2}));

	case {'corr', 'sac_corr'}
		sac = varargin{1};
		cell_stat(k) = s.'*sac / sqrt(sum(s.^2) * sum(sac.^2));

	case {'corr_unrml', 'corr_u', 'corr_uu'}
		sac = varargin{1};
		cell_stat(k) = s.'*sac;% / sqrt(sum(s.^2) * sum(sac.^2));

    case {'corr_sac_sac', 'sac_sac_2'}
        sac = varargin{1};
        cell_stat(k) = self(cell_info(k), 'sac_corr', use_normalized_strat, sac{1}).^2 ...
                     + self(cell_info(k), 'sac_corr', use_normalized_strat, sac{2}).^2;

    case {'prcntile', 'ptile'}
        cell_stat(k) = get_percentile([x s], property_arg{1});

    case {'asym_2an_p', 'asym2_p'}
        cell_stat(k) = log(norm(cell_info_elem.asymm_2an_prj));
        if (cell_info_elem.is_cutoff)
            cell_stat(k) = NaN;
        end

    case {'density'}
        cell_stat(k) = cell_info_elem.area_projection / cell_info_elem.area_hull;

    case {'size'}
        cell_stat(k) = cell_info_elem.area_hull;

    case {'maxdiameter'}
        cell_stat(k) = cell_info_elem.max_diameter;
        if ~isempty(property_arg)
            if ~strcmp(property_arg{1}, 'um')   % in um
                error('invalid arg');
            end
            cell_stat(k) = cell_stat(k) * 16.5 * 4 / 1000; % in um
        end


    % case {'sus/tran'}
    % 	cell_stat(k) = cell_info_get_strat_property(cell_info())
    %     if at
    %         cell_stat(k) = log(as/at);
    %     else
    %         cell_stat(k) = 20;
    %     end

    otherwise
        error('not recognized stat name') 

    end
end
