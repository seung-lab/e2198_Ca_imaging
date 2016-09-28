function cell_info_polarplot_pref_dir(cell_info,ca_dsos, onoff, normalize, vertalign, varargin)
% onoff: 0 = total, 1 = on, 2 = off, [1 2] = separate on and off

dir_sac2gc='sac2gc/';

type_names={'37r','37d','37c','37v','7o','7id','7ic','7iv','2aw','63','73'};
type_names={'37c','37v','37r','37d','7o','7iv','7ir','7id','2aw','63','73'};
idx_panel=[1 1 1 1 2 2 2 2  3 4 5];
%idx_color=[1 6 2 3 1 6 2 3  3 3 3];
idx_color=[4 3 2 5 4 3 2 5  8 8 8];
type_onoff={0 0 0 0  0 0 0 0  0 1 0};
%r_lim = []
titles = {'', '', '', '37', '', '', '', '7o / 7i','2aw','63 (SAC = On only)','73'};
fontsize = 12;
if exist('onoff', 'var') && ~isempty(onoff)
    type_onoff = repmat({onoff}, size(type_names));
else
    onoff = -10;
end
if ~exist('normalize', 'var')
    normalize = 0;
end
if ~exist('vertalign', 'var')
    vertalign = 0;
end

if vertalign
fignumrow=max(idx_panel);
fignumcol=1+2*length(onoff)
else
fignumcol=max(idx_panel);
fignumrow=1+2*length(onoff);
fignumrow=1+4*length(onoff);
end
colors=distinguishable_colors(10);
figure(varargin{:});
clf
ax = gca();
colors = ax.ColorOrder;
colors(end+1,:) = [0 0 0];

clf;

strat=cell_info_bin_strat(cell_info,0);
% rebin
for k = 1:length(strat)
    xx = strat{k}(1:end-7,:);  % discarding last few bins
    strat{k}=squeeze(mean(reshape(xx, 4, [], 2)));
end

cid=cell_info(find(strcmp({cell_info.type},type_names{1}),1)).cell_id;
strat_x=strat{cid}(:,1);


contact_file_path = 'gc_sac_contacts_20160610.mat';
contact_file_path = 'gc_sac_contacts.20160615.mat';
%contact_file_path = 'gc_sac_contacts.2um.20160621.mat';
%contact_file_path = 'gc_sac_contacts.0.5um.20160622.mat';
%contact_file_path = 'gc_sac_contacts.20160719.noGCOnOffSeparation.mat';
load(contact_file_path);    % gc_denom, gc_num


for j=1:numel(type_names)

    rowoffset = 0;
    % stratification 
    if vertalign
        subplot(fignumrow,fignumcol, (idx_panel(j)-1)*fignumcol+1 );
    else
    subplot(fignumrow,fignumcol,idx_panel(j));
    end
    cell_ids=[cell_info(strcmp({cell_info.type},type_names{j})).cell_id];
    
    strat_y=[];
    for i=1:numel(cell_ids)
        strat_y(:,i)=strat{cell_ids(i)}(:,2);
    end
    plot(strat_x/100,mean(strat_y,2)*100,'LineWidth',1,'Color',colors(idx_color(j),:));
    if mod(j,4)==0
        legend(type_names(j-3:j))
    end
    hold on;
    ax=gca;
    ax.FontSize = fontsize;
    %%{
    ax.GridLineStyle = '--';
    ax.GridAlpha = 0.3;
    %}
    ax.XLim=[0 1];
    ax.YTick = [];
    ax.XTick = [0 0.28 0.62 1];
    ax.XTickLabel(2:3) = {'', ''};
    %ax.XTickLabel(2:3) = {'Off SAC', 'On SAC'};
    %ax.XTickLabelRotation = 45;
    grid on
    if j>8
        title(type_names{j})
    end
    title(titles{j});
    
    for layer = type_onoff{j}(:).' %onoff(:).'
        if layer==0
            layer = [1 2];
        end
        if isequal(layer, 2)
            luminance = 0.6;
        else
            luminance = 1;
        end

    % sac input
    if vertalign
        coloffset = [1];
        subplot(fignumrow,fignumcol, (idx_panel(j)-1)*fignumcol+coloffset+1 );
    else
    rowoffset = [1 2];
    subplot(fignumrow,fignumcol,rowoffset*fignumcol+idx_panel(j));
    end
    
    for i=1:numel(cell_ids)
        %{
        %Old data by Matt
        search_pattern=sprintf('%s/sac2%d_gc*data.mat',dir_sac2gc,cell_ids(i));
        search_result=dir(search_pattern);
        if isempty(search_result)
            continue;
        end
        file_path=sprintf('%s/%s',dir_sac2gc,search_result.name);
        load(file_path);
        %}
        %%{
        % new data
        %display('new')
        % idx=find(cell_ids(i)==gc_denom(1,:));
        % angle_denom = gc_denom(2:end, idx);
        % idx=find(cell_ids(i)==gc_num(1,:));
        % angle_num = gc_num(2:end, idx);
        idx=find(cell_ids(i)==gc_denom_keys);
        angle_denom = gc_denom_vals{idx};
        idx=find(cell_ids(i)==gc_num_keys);
        angle_num = gc_num_vals{idx};
        %}

        theta=zeros(8,1);
        bin_sum_num=zeros(8,1);
        bin_sum_denom=zeros(8,1);

        for tt=0:8
            idx=mod(tt,8)+1;
            idx_angles=max(22+45*(tt-1)+1,1):min(22+45*tt,360);
            theta(idx)=(idx-1)*45*pi/180;
            bin_sum_num(idx)=bin_sum_num(idx)+nansum(nansum(angle_num(idx_angles,layer)));
            bin_sum_denom(idx)=bin_sum_denom(idx)+nansum(nansum(angle_denom(idx_angles,layer)));
        end
        %%{
        theta = pi/4 * [0:7].';
        bin_sum_num2 = bin_sum_num;
        bin_sum_denom2 = bin_sum_denom;
        bin_sum_num = rebin(angle_num);
        bin_sum_denom = rebin(angle_denom);
        bin_sum_num = sum(bin_sum_num(:,layer), 2);
        bin_sum_denom = sum(bin_sum_denom(:,layer), 2);
        %in_sum_denom
        %bin_sum_denom2
        assert(isequal(bin_sum_num2, bin_sum_num))
        assert(isequal(bin_sum_denom2, bin_sum_denom))
        %}
        %{
        % trying out 24 bins instead of 8 bins
        theta = pi/12 * [0:23].';
        bin_sum_num = rebinTo24(angle_num);
        bin_sum_denom = rebinTo24(angle_denom);
        bin_sum_num = sum(bin_sum_num(:,layer), 2);
        bin_sum_denom = sum(bin_sum_denom(:,layer), 2);
        %}
        binned_sac_input_rho=bin_sum_num./bin_sum_denom;
        [x,y]=pol2cart(theta,binned_sac_input_rho);
        theta_all = theta;
        [theta,rho]=cart2pol(sum(x),sum(y));  % !! this is doing a vec sum and not vec mean
        if normalize
            rho = rho / mean(binned_sac_input_rho);  % wrong. sum, not mean
        end

        theta=pi*3/2-theta + pi/2;  % to final "standard" coord, SAC dir (not "predicted dir")

        polarplot([0 theta],[0 rho],'LineWidth',1,'Color',luminance*colors(idx_color(j),:));

        hold on;
    end
    ax=gca;
    ax.FontSize = fontsize;
    ax.ThetaTick=0:45:315;
    ax.GridAlpha = 1;
    ax.GridColor = 0.8*[1 1 1];
    %ax.ThetaTickLabel(2:2:8)={''};
    ax.ThetaTickLabel(3:8)={''};
    ax.GridLineStyle = '--';
    lim = rlim();
    if normalize
        if lim(2)>0.03
            rlim([0 0.1])
        else
            rlim([0 0.02])
        end
        rlim([0 3])
    else
        if lim(2)>0.03
            rlim([0 0.1])
            ax.RTick = [0.05 0.1];
        else
            rlim([0 0.025])
            ax.RTick = [0.01 0.02];
            ax.RTick = [0.025];
        end
    end


    % physiology    
        %%{
        if strcmp(type_names{j}, '63')
            layer = [1 2];
        end
        %}
    if vertalign
        coloffset = [2];
        subplot(fignumrow,fignumcol, (idx_panel(j)-1)*fignumcol+coloffset+1 );
    else
    rowoffset = [3 4];
    subplot(fignumrow,fignumcol,rowoffset*fignumcol+idx_panel(j));
    end
    cell_ids=cell_ids(ismember(cell_ids,ca_dsos.omni_id));
    
    for i=1:numel(cell_ids)
        idx=find(cell_ids(i)==ca_dsos.omni_id);
        [xx, yy] = pol2cart(ca_dsos.ds_theta(idx,layer),ca_dsos.ds_r(idx,layer));
        [theta,rho]=cart2pol(sum(xx),sum(yy));
        if normalize
            rho = rho / sum(ca_dsos.r_mean(idx,layer));
        end
        
        theta=pi/2+theta;  % to final "standard" coord
        polarplot([0 theta],[0 rho],'LineWidth',1,'Color',luminance*colors(idx_color(j),:));    
        hold on;
    end
    ax=gca;
    ax.FontSize = fontsize;
    ax.ThetaTick=0:45:315;
    ax.GridAlpha = 1;
    ax.GridColor = 0.8*[1 1 1];
    %ax.ThetaTickLabel(2:2:8)={''};
    ax.ThetaTickLabel(3:8)={''};
    ax.GridLineStyle = '--';
    lim = rlim();
    if normalize
        if lim(2)>0.5
            rlim([0 1])
            ax.RTick = [0.5 1];
        else
            rlim([0 0.5])
            ax.RTick = [0.5];
        end
    else
        if lim(2)>6
            rlim([0 15])
            ax.RTick = [5 10 15];
        else
            rlim([0 6])
            ax.RTick = [6];
        end
    end

    end % for layer
    
end


end


function x0to7 = rebin(x1to360)
    if any(isnan(x1to360))
        error('NaN values found')
    end
    x0to7 = squeeze(sum(reshape(circshift(x1to360, 23, 1), 45, 8, 2)));
end


function x0to7 = rebinTo24(x1to360)
    if any(isnan(x1to360))
        error('NaN values found')
    end
    x0to7 = squeeze(sum(reshape(circshift(x1to360, 8, 1), 15, 24, 2)));
end
