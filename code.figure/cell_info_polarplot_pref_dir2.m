
function cell_info_polarplot_pref_dir(cell_info,ca_dsos, onoff, type_list)
% onoff: 0 = total, 1 = on, 2 = off, [1 2] = separate on and off
% [0 1 2]

if exist('type_list', 'var')
    print_vals = 1;
else
    print_vals = 0;
    gc_types = list_types(cell_info);
    type_list = [gc_types(:).' {'37' '7i'}];
end


fignumcol=6;
fignumrow=length(onoff);


colors=distinguishable_colors(40);


contact_file_path = 'gc_sac_contacts_20160610.mat';
contact_file_path = 'gc_sac_contacts.20160615.mat';
%contact_file_path = 'gc_sac_contacts.2um.20160621.mat';
%contact_file_path = 'gc_sac_contacts.0.5um.20160622.mat';
%contact_file_path = 'gc_sac_contacts.20160719.noGCOnOffSeparation.mat';
load(contact_file_path);    % gc_denom, gc_num


%for celltype = {'2aw'}
%for celltype = {'weirdos'}
for celltype = type_list(:).'
%for celltype = {'37'}
%for celltype = {'7i'}

    celltype = celltype{1};
    disp(celltype)

	cells = get_cell_info(cell_info, celltype);
    cell_ids = [cells.cell_id];
    %cell_ids_w_ca =cell_ids(ismember(cell_ids,ca_dsos.omni_id));
    %cell_ids = [cell_ids_w_ca setdiff([cells.cell_id], cell_ids_w_ca)];

summary_fig_h = figure(30);
clf;

    rowoffset = -1;

    for layer = onoff
        if layer==0
            layer = [1 2];
            titletext = 'on+off combined';
            luminance = 1;
        elseif isequal(layer, 2)
            luminance = 0.6;
            luminance = 1;
            titletext = 'off';
        else
            luminance = 1;
            titletext = 'on';
        end
        rowoffset = rowoffset+1;

    lines = {[], []};

    col = 0;

    % physiology
    col = col+1;
    subplot(fignumrow,fignumcol,rowoffset*fignumcol+col);

    ca_dirs = [];
    no_ca = [];
    %for k=1:length(cells)
    %   cell_id = cells(k).cell_id;
    for i=1:numel(cell_ids)
        idx=find(cell_ids(i)==ca_dsos.omni_id);
        if isempty(idx)
            %no_ca = true; linestyle = 'none';
            %display('empty')
            no_ca(end+1) = i;
            continue;
        else
            %no_ca = false; linestyle = ':';
        end
        [xx, yy] = pol2cart(ca_dsos.ds_theta(idx,layer),ca_dsos.ds_r(idx,layer));
        [theta,rho]=cart2pol(sum(xx),sum(yy));
        rho_mean = sum(ca_dsos.r_mean(idx,layer));
        rho=min(rho,20);  % to limit one due-to-noise outlier in 7iv 
        
        if ~isempty(idx)
            ca_dirs(end+1,:) = [cell_ids(i) rho/rho_mean rho];
        end

        theta=3/2*pi-theta;  % to match with the angle of sac input

        subplot(fignumrow,fignumcol,rowoffset*fignumcol+col);
        polarplot([0 theta],[0 rho],'LineWidth',1,'Color',luminance*colors(i,:));
        title('Calcium')
        hold on;

        subplot(fignumrow,fignumcol,rowoffset*fignumcol+col+1);
        polarplot([0 theta],[0 rho/rho_mean],'LineWidth',1,'Color',luminance*colors(i,:));
        title('Calcium Normalized')
        hold on;
    end
    ax=gca;
    ax.ThetaTick=0:45:315;
    
    %title(titletext);
    ax = gca();
    ax = axes('Position', ax.Position);
    %ax.Visible = 'Off';
    %axis off
    %ax.XAxis.Visible = 'Off';
    %ax.YTick = [];
    ax.Color = 'none';
    ylabel(titletext);

    col = col+1;


    if print_vals
        disp('Ca');
        num2str(sortrows(ca_dirs, 2))
    end

    % sac input
    %rowoffset = rowoffset+1;
    col = col+1;
    subplot(fignumrow,fignumcol,rowoffset*fignumcol+col);

    sacdirs = [];
    lines = [];

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
        if isempty(idx)
            warning(sprintf('cell not found %d in %s', cell_ids(i), celltype)); % only weirdos
            continue;
        end
        angle_denom = gc_denom_vals{idx};
        idx=find(cell_ids(i)==gc_num_keys);
        angle_num = gc_num_vals{idx};
        %}

        theta = pi/4 * [0:7].';
        bin_sum_num = rebin(angle_num);
        bin_sum_denom = rebin(angle_denom);
        bin_sum_num = sum(bin_sum_num(:,layer), 2);
        bin_sum_denom = sum(bin_sum_denom(:,layer), 2);

        binned_sac_input_rho=bin_sum_num./bin_sum_denom;
        [x,y]=pol2cart(theta,binned_sac_input_rho);
        [theta_meanvec,rho_meanvec]=cart2pol(mean(x),mean(y));

        theta_meanvec = pi+ theta_meanvec;
        theta=pi+theta;  % match Matt's angles

        sacdirs(end+1,:) = [cell_ids(i) rad2deg(theta_meanvec) rho_meanvec/mean(binned_sac_input_rho)];

        subplot(fignumrow,fignumcol,rowoffset*fignumcol+col);
        lines(i, 1) = polarplot([0 theta_meanvec],[0 rho_meanvec], 'LineWidth',1,'Color',luminance*colors(i,:));
        title('SAC Unnormalized')
        hold on;

        subplot(fignumrow,fignumcol,rowoffset*fignumcol+col+1);
        lines(i, 2) = polarplot([0 theta_meanvec],[0 rho_meanvec/mean(binned_sac_input_rho)],'LineWidth',1,'Color',luminance*colors(i,:));
        title('SAC Normalized')
        hold on;

        subplot(fignumrow,fignumcol,rowoffset*fignumcol+col+2);
        lines(i, 3) = polarplot(theta([end 1:end]), binned_sac_input_rho([end 1:end]),'LineWidth',1,'Color',luminance*colors(i,:));
        title('SAC')
        hold on;

        subplot(fignumrow,fignumcol,rowoffset*fignumcol+col+3);
        lines(i, 4) = polarplot(theta([end 1:end]), binned_sac_input_rho([end 1:end])./max(binned_sac_input_rho),'LineWidth',1,'Color',luminance*colors(i,:));
        title('SAC w max=1')
        hold on;
    end
    ax=gca;
    ax.ThetaTick=0:45:315;
    
    if ~isempty(lines(no_ca, :))
        set(lines(no_ca, :), 'LineStyle', '--');
    end

        if print_vals
            disp('SAC');
            num2str(sortrows(sacdirs, 2))
        end
    end % for layer
title(celltype)
%refresh
%drawnow
%close(gcf());


    folder = './tmp_summary';
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
    filepath = sprintf('%s\\%s.png',folder,['type_dirs_' celltype]);
    print(summary_fig_h, '-r300', filepath, '-dpng');

end


end


function x0to7 = rebin(x1to360)
    if any(isnan(x1to360))
        error('NaN values found')
    end
    x0to7 = squeeze(sum(reshape(circshift(x1to360, 23, 1), 45, 8, 2)));
end

