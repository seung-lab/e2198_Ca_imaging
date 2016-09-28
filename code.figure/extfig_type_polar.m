function cell_info_polarplot_pref_dir(cell_info, ca_dsos, tuning_ordered, onoff, normalize, vertalign, varargin)
% onoff: 0 = total, 1 = on, 2 = off, [1 2] = separate on and off

type_names={'37r','37d','37c','37v','7o','7id','7ic','7iv','2aw','63','73'};
type_names={'37c','37v','37r','37d','7o','7iv','7ir','7id','2aw','63','73'};
%idx_panel=[1 1 1 1 2 2 2 2  3 4 5];
%idx_color=[1 6 2 3 1 6 2 3  3 3 3];
%idx_color=[4 3 2 5 4 3 2 5  8 8 8];
type_onoff={0 0 0 0  0 0 0 0  0 1 0};

%titles = {'', '', '', '37', '', '', '', '7o / 7i','2aw','63 (On only)','72n'};
fontsize = 12;
if exist('onoff', 'var') && ~isempty(onoff)
    type_onoff = repmat({onoff}, size(type_names));
else
    onoff = -10;
end
%{
if ~exist('normalize', 'var')
    normalize = 0;
end
if ~exist('vertalign', 'var')
    vertalign = 0;
end
%}


fignumcol=4;
gridnumrow=2*length(onoff);
numgridrow = ceil(length(type_names) / fignumcol);
fignumrow = numgridrow * gridnumrow;


colors=distinguishable_colors(30);
%fig_h = figure(varargin{:});
fig_h = figure('Visible', 'off');
clf

fig_h.Visible = 'off';
fig_h.Position(3:4) = [1201 1855]

contact_file_path = 'gc_sac_contacts_20160610.mat';
contact_file_path = 'gc_sac_contacts.20160615.mat';
%contact_file_path = 'gc_sac_contacts.2um.20160621.mat';
%contact_file_path = 'gc_sac_contacts.0.5um.20160622.mat';
%contact_file_path = 'gc_sac_contacts.20160719.noGCOnOffSeparation.mat';
load(contact_file_path);    % gc_denom, gc_num


for j=1:numel(type_names)

    rowoffset = 0;
    col = mod(j-1, fignumcol)+1;
    gridrow = ceil(j / fignumcol);

	cells = get_cell_info(cell_info, type_names{j});
    cell_ids = [cells.cell_id];
    
    for layer = type_onoff{j}(:).'
        if layer==0
            layer = [1 2];
        end
        if isequal(layer, 2)
            luminance = 0.6;
        else
            luminance = 1;
        end

    no_ca = [];
    lines = [];

    for i=1:numel(cell_ids)
        idx=find(cell_ids(i)==ca_dsos.omni_id);
        if isempty(idx)
            no_ca(end+1) = i;
            ca_id = [];
        else
        	ca_id = ca_dsos.ca_id(idx);
        end


	    % sac input
	    rowoffset = [0];
	    subplot(fignumrow,fignumcol, ((gridrow-1)*gridnumrow+rowoffset)*fignumcol+col);

        idx=find(cell_ids(i)==gc_denom_keys);
        angle_denom = gc_denom_vals{idx};
        idx=find(cell_ids(i)==gc_num_keys);
        angle_num = gc_num_vals{idx};

        theta = pi/4 * [0:7].';
        bin_sum_num = rebin(angle_num);
        bin_sum_denom = rebin(angle_denom);
        bin_sum_num = sum(bin_sum_num(:,layer), 2);
        bin_sum_denom = sum(bin_sum_denom(:,layer), 2);

        binned_sac_input_rho=bin_sum_num./bin_sum_denom;
        %{
        [x,y]=pol2cart(theta,binned_sac_input_rho);
        theta_all = theta;
        [theta,rho]=cart2pol(sum(x),sum(y));
        if normalize
            rho = rho / mean(binned_sac_input_rho);
        end
        %}

        theta=pi*3/2-theta + pi/2;  % to final "standard" coord, SAC dir (not "predicted dir")

        lines(i, 1) = polarplot(theta([end 1:end]), binned_sac_input_rho([end 1:end]),'LineWidth',1,'Color',luminance*colors(i,:));
        hold on;

        % physiology
        if isempty(ca_id)
        	continue;
        end
	    rowoffset = [1];
	    subplot(fignumrow,fignumcol, ((gridrow-1)*gridnumrow+rowoffset)*fignumcol+col);
	    theta = pi/4 * [1:8].';
	    theta = pi/2+theta;  % to final "standard" coord
	    rho = tuning_ordered(layer,:,ca_id);
	    %rho = tuning_ordered(3,:,ca_id);
	    rho = sum(rho,1);
	    polarplot(theta([end 1:end]), rho([end 1:end]),'LineWidth',1,'Color',luminance*colors(i,:));
        hold on;
    end

    if ~isempty(lines(no_ca, :))
        set(lines(no_ca, :), 'LineStyle', '--');
    end

	end % for layer

	% sac input
	rowoffset = 0;
	subplot(fignumrow,fignumcol, ((gridrow-1)*gridnumrow+rowoffset)*fignumcol+col);
    title(type_names{j})
    ax = gca;
	%applystyle(gca());
	lim = rlim();
        if lim(2)>0.03
            rlim([0 0.07])
            ax.RTick = [0.05 0.1];
        else
            rlim([0 0.025])
            ax.RTick = [0.025];
        end

  	% physiology
	rowoffset = 1;
	subplot(fignumrow,fignumcol, ((gridrow-1)*gridnumrow+rowoffset)*fignumcol+col);
    ax = gca;
	%applystyle(gca());
	lim = rlim();
	if j<=8
        if lim(2)>=30
            rlim([0 60])
            ax.RTick(1) = []; % remove 0 tick
        else
            rlim([0 20])
            ax.RTick(1) = []; % remove 0 tick
        end
    else
            ax.RTick(1) = []; % remove 0 tick
    end
end

for fig = 1:fignumrow*fignumcol
    subplot(fignumrow,fignumcol, fig);

    ax=gca;
	if ~strcmp(class(ax), 'matlab.graphics.axis.PolarAxes')
		continue;
	end
    ax.FontSize = fontsize;
    ax.ThetaTick=0:45:315;
    ax.GridAlpha = 1;
    ax.GridColor = 0.8*[1 1 1];
    %ax.ThetaTickLabel(2:2:8)={''};
    ax.ThetaTickLabel(3:8)={''};
    ax.GridLineStyle = '--';
    %{
    lim = rlim();
    if normalize
        if lim(2)>0.03
            rlim([0 0.1])
        else
            rlim([0 0.02])
        end
        rlim([0 3])
    else
        if lim(2)>=30  % physiology
            rlim([0 60])
            ax.RTick(1) = []; % remove 0 tick
        elseif lim(2)>1  % physiology
            rlim([0 20])
            %ax.RTick = [0.05 0.1];
            ax.RTick(1) = []; % remove 0 tick
        elseif lim(2)>0.03
            rlim([0 0.07])
            ax.RTick = [0.05 0.1];
        else
            rlim([0 0.025])
            ax.RTick = [0.025];
        end
    end
    %}
end

fig_h.Position


    file = sprintf('../fig/extfig6');
    print(fig_h, '-r300', file, '-dpng');
    print(fig_h, '-r300', file, '-dsvg');
    print(fig_h, file, '-depsc');
    close(fig_h);

end %function


function x0to7 = rebin(x1to360)
    if any(isnan(x1to360))
        error('NaN values found')
    end
    x0to7 = squeeze(sum(reshape(circshift(x1to360, 23, 1), 45, 8, 2)));
end


function x0to7 = applystyle(ax)
	fontsize = 12;

	ax.FontSize = fontsize;
    ax.ThetaTick=0:45:315;
    ax.GridAlpha = 1;
    ax.GridColor = 0.8*[1 1 1];
    %ax.ThetaTickLabel(2:2:8)={''};
    ax.ThetaTickLabel(3:8)={''};
    ax.GridLineStyle = '--';
end