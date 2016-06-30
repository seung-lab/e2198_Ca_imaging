function cell_info_polarplot_pref_dir(cell_info,ca_dsos, onoff)
% onoff: 0 = total, 1 = on, 2 = off, [1 2] = separate on and off

dir_sac2gc='sac2gc/';

type_names={'37r','37v','37d','37c','7o','7iv','7id','7ic','2an','2aw','2o','1no', '2i'};
idx_panel=[1 1 1 1 2 2 2 2 3 4 5 6 7];
idx_color=[1 3 6 2 1 3 6 2 3 3 3 3 1];

fignumcol=max(idx_panel);
fignumrow=1+2*length(onoff);
colors=distinguishable_colors(10);
figure();
clf;

strat=cell_info_bin_strat(cell_info,0);
cid=cell_info(find(strcmp({cell_info.type},type_names{1}),1)).cell_id;
strat_x=strat{cid}(:,1);


contact_file_path = 'gc_sac_contacts_20160610.mat';
contact_file_path = 'gc_sac_contacts.20160615.mat';
%contact_file_path = 'gc_sac_contacts.2um.20160621.mat';
%contact_file_path = 'gc_sac_contacts.0.5um.20160622.mat';
load(contact_file_path);    % gc_denom, gc_num


for j=1:numel(type_names)

    rowoffset = 0;
    % stratification 
    subplot(fignumrow,fignumcol,idx_panel(j));
    cell_ids=[cell_info(strcmp({cell_info.type},type_names{j})).cell_id];
    
    strat_y=[];
    for i=1:numel(cell_ids)
        strat_y(:,i)=strat{cell_ids(i)}(:,2);
    end
    plot(strat_x,mean(strat_y,2),'LineWidth',1,'Color',colors(idx_color(j),:));
    hold on;
    ax=gca;
    ax.XLim=[0 100]; 
    if j>8
        title(type_names{j})
    end
    
    for layer = onoff(:).'
        if layer==0
            layer = [1 2];
        end
        if isequal(layer, 2)
            luminance = 0.6;
        else
            luminance = 1;
        end

    % physiology    
    rowoffset = rowoffset+1;
    subplot(fignumrow,fignumcol,rowoffset*fignumcol+idx_panel(j));
    cell_ids=cell_ids(ismember(cell_ids,ca_dsos.omni_id));
   
    for i=1:numel(cell_ids)
        idx=find(cell_ids(i)==ca_dsos.omni_id);
        [xx, yy] = pol2cart(ca_dsos.ds_theta(idx,layer),ca_dsos.ds_r(idx,layer));
        [theta,rho]=cart2pol(sum(xx),sum(yy));
        rho=min(rho,20);  % to limit one due-to-noise outlier in 7iv 
        
        theta=3/2*pi-theta;  % to match with the angle of sac input
        polarplot([0 theta],[0 rho],'LineWidth',1,'Color',luminance*colors(idx_color(j),:));    
        hold on;
    end
    ax=gca;
    ax.ThetaTick=0:45:315;
    lim = rlim();
    if lim(2)>5
        rlim([0 5])
    end
    %{
    % tmp scaling for presentation
    if idx_panel(j)==2
        rlim([0 3])
    end
    %}
    
    % sac input
    rowoffset = rowoffset+1;
    subplot(fignumrow,fignumcol,rowoffset*fignumcol+idx_panel(j));
    
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
        %%{
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
        [theta,rho]=cart2pol(sum(x),sum(y));

        theta=pi+theta;  % match Matt's angles

        polarplot([0 theta],[0 rho],'LineWidth',1,'Color',luminance*colors(idx_color(j),:));

        if rho > 0.2 %|| strcmp(type_names{j}, '37v')
            warning('too big')
            type_names(j)
            cell_ids(i)
            [theta*180/pi,rho]
            %%{
            %[theta_all binned_sac_input_rho]
            %[x y]
            num2str([binned_sac_input_rho bin_sum_num bin_sum_denom])
            %}
        end

        hold on;
    end
    ax=gca;
    ax.ThetaTick=0:45:315;
    lim = rlim();
    if lim(2)>0.2
        rlim([0 0.2])
        % tmp scaling for presentation
        rlim([0 0.1])
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
