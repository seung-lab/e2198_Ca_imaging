function cell_info_polarplot_pref_dir(cell_info,ca_dsos)

dir_sac2gc='/data/research/jk/e2198/e2198_Ca_imaging/sac2gc/gc_all/';

type_names={'37r','37v','37d','37c','7o','7iv','7id','7ic','2an','2aw','2o','1no'};
idx_panel=[1 1 1 1 2 2 2 2 3 4 5 6];
idx_color=[1 3 6 2 1 3 6 2 3 3 3 3];

fignumcol=max(idx_panel);
fignumrow=3;
colors=distinguishable_colors(10);
figure(1);
clf;

strat=cell_info_bin_strat(cell_info,0);
cid=cell_info(find(strcmp({cell_info.type},type_names{1}),1)).cell_id;
strat_x=strat{cid}(:,1);

for j=1:numel(type_names)

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
    
    % physiology    
    subplot(fignumrow,fignumcol,fignumcol+idx_panel(j));
    cell_ids=cell_ids(ismember(cell_ids,ca_dsos{:,6}));
   
    for i=1:numel(cell_ids)
        idx=find(cell_ids(i)==ca_dsos{:,6});
        % Shang: check if this is correct to combine on/off responses
        [x_on,y_on]=pol2cart(ca_dsos{idx,2}{1}(1),ca_dsos{idx,1}{1}(1));
        [x_off,y_off]=pol2cart(ca_dsos{idx,2}{1}(2),ca_dsos{idx,1}{1}(2));
        [theta,rho]=cart2pol(x_on+x_off,y_on+y_off);
        rho=min(rho,20);  % to limit one due-to-noise outlier in 7iv 
        theta=3/2*pi-theta;  % to match with the angle of sac input
        polarplot([0 theta],[0 rho],'LineWidth',1,'Color',colors(idx_color(j),:));    
        hold on;
    end
    ax=gca;
    ax.ThetaTick=0:45:315;
    
    % sac input
    subplot(fignumrow,fignumcol,2*fignumcol+idx_panel(j));
    
    for i=1:numel(cell_ids)
        search_pattern=sprintf('%s/sac2%d_gc*data.mat',dir_sac2gc,cell_ids(i));
        search_result=dir(search_pattern);
        if isempty(search_result)
            continue;
        end
        file_path=sprintf('%s/%s',dir_sac2gc,search_result.name);
        load(file_path);

        theta=zeros(8,1);
        bin_sum_num=zeros(8,1);
        bin_sum_denom=zeros(8,1);

        for tt=0:8
            idx=mod(tt,8)+1;
            idx_angles=max(22+45*(tt-1)+1,1):min(22+45*tt,360);
            theta(idx)=(idx-1)*45*pi/180;
            bin_sum_num(idx)=bin_sum_num(idx)+nansum(nansum(angle_num(idx_angles,:)));
            bin_sum_denom(idx)=bin_sum_denom(idx)+nansum(nansum(angle_denom(idx_angles,:)));
        end
        binned_sac_input_rho=bin_sum_num./bin_sum_denom;
        [x,y]=pol2cart(theta,binned_sac_input_rho);
        [theta,rho]=cart2pol(sum(x),sum(y));
        polarplot([0 theta],[0 rho],'LineWidth',1,'Color',colors(idx_color(j),:));
        hold on;
    end
    ax=gca;
    ax.ThetaTick=0:45:315;
    
end


end

