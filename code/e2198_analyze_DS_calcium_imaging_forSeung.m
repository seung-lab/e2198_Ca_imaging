

OVimg_fn = 'AVG_043_p2_g7_ZoomedOut_center.tif';
% roi_fn = 'AVG_043_p2_g7_ZoomedOut_center_011.rpb';
roi_fn = 'AVG_4classes_043_p2_g7_ZoomedOut_center_011.rpb';
stim_fn = 'sDS 8x45deg, narrow, 4.0s.QDS';
data_fn = ['035_p2_g7_DSbars_200um_export.mat';
           '036_p2_g7_DSbars_200um_export.mat';
           '037_p2_g7_DSbars_200um_export.mat';
           '038_p2_g7_DSbars_200um_export.mat';
           '039_p2_g7_DSbars_200um_export.mat';
           '040_p2_g7_DSbars_200um_export.mat';
           '041_p2_g7_DSbars_200um_export.mat';
           '042_p2_g7_DSbars_200um_export.mat';
           '034_p2_g7_DSbars_200um_export.mat']

typeclrs = [[0 125 125];...
           [0 255 255];...
           [0 0 255];...
           [255 0 255];...
           [255 0 100];...
           [255 255 0];...
           [255 150 0];...
           [0 255 0]];
       
% Build up data array from ROIs
load(data_fn(1,:));

n_rois = size(roi_sums,2)

stim_struct = load_stim(stim_fn);
[n_samples,t_min,t_max,T,dt,fs] = time_info(t_frame);
t_pre_frames = round(stim_struct(1).t_pre/dt);
t_stim_frames = round(stim_struct(1).t_stim/dt);
t_post_frames = round(stim_struct(1).t_post/dt);
nconds = 8; t_trial=31; n_reps = 5;
roi_sums_all = zeros(nconds*t_trial*n_reps ,n_rois);

for count = 1:size(data_fn,1)
       
    load(data_fn(count,:));
    
    xdata = t_scan'/1000;
    ydata = double(ch1_t_stim);
    xtimes = crossing_times(xdata,ydata,200);
    xtimes = xtimes(1:2:end);
    temp3 = abs(repmat(t_frame',[1 length(xtimes)])-repmat(xtimes',[length(t_frame) 1]));
    [x,allstimframes] = min(temp3);
    %%% Clip out each trial from the raw data

    [q,w] = find(sum(roi_sums,1)~=0);
    roi_sums_time_avg=mean(roi_sums,1);
    roi_sums=100*(roi_sums./repmat(roi_sums_time_avg,[size(roi_sums,1) 1])-1);
    
    for count2=1:length(w)
        thisroi = roi_sums(:,w(count2));
        thisroi_all = roi_sums_all(:,w(count2));
        thisroi_clip = [];
        %%% Clip out each trial from the raw data using stim times
        for repidx = 1:n_reps
            for condidx = 1:nconds
                     thisroi_clip = [thisroi_clip ; thisroi(allstimframes((repidx-1)*nconds+condidx):allstimframes((repidx-1)*nconds+condidx)+t_trial-1,:)];
            end
        end
%         if(sum(thisroi_all)==0)
            roi_sums_all(:,w(count2)) = thisroi_clip;
%         elseif max(thisroi_clip(200:end))>max(thisroi_all(200:end))
%             roi_sums_all(:,w(count2)) = thisroi_clip;
%         end
    end

end
allstimframes = [1:t_trial:size(roi_sums_all,1)];

figure; imagesc(roi_sums_all); axis square; colorbar

stim_struct = load_stim(stim_fn);
nconds = double(stim_struct(1).nconds);
if nconds>0
   t_pre_frames = round(stim_struct(1).t_pre/dt);
   t_stim_frames = round(stim_struct(1).t_stim/dt);
   t_post_frames = round(stim_struct(1).t_post/dt);
   t_trial = t_pre_frames+t_stim_frames+t_post_frames;
end

% nreps = length(allstimframes)./nconds;
% cond_order = [4 6 8 1 3 5 7 2];
%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%% These are correctly ordered for the EM images in which you're looking
%%% down through the IPL at the GCL
cond_order = [5 3 1 8 6 4 2 7]; 

[dx,dy,channels,frames,ch1,ch2] = read_tif(OVimg_fn);

figure; imagesc(ch1'); colormap('gray'); axis square; set(gca,'CLim',[0 40])
roi_struct = load_rois_from_rpb(roi_fn,typeclrs./255);
for count=1:length(roi_struct.roi_ids)
    roi_borders{count} = [get(roi_struct.border_h(count),'XData')' get(roi_struct.border_h(count),'YData')'];
end


%%
% Make summary figure for each cell
figdir = 'C:\Users\Kevin\Desktop\CalciumImaging\090708_epor\090708_epor_cell_summary';
for count=1:634
    thislabel = show_id_to_alpha(count)
    summary_fig_h = ...
        figure('Tag','summary_fig_h',...
        'Name',sprintf('cell_%s',thislabel),...
        'NumberTitle','off',...
        'PaperPositionMode','auto',...
        'InvertHardcopy','off',...
        'DoubleBuffer','on',...
        'Position',[341 79 465 619]);
    
    ymax = 0;
    ymin = 0;
    % subplot(double(nconds),2,[1 3 5],'Parent',summary_fig_h);
    colormap('gray(256)');
    temp = ch1;
    n_cols = size(temp,1); n_rows = size(temp,2);
    image_axes_h = axes('Parent',summary_fig_h,...
        'Tag','image_axes_h',...
        'YDir','reverse',...
        'DrawMode','normal',...
        'Visible','off',...
        'Units','Normalized',...
        'Position',[0.1 0.50 0.35 0.35],...
        'XLim',[0.5,n_cols+0.5],...
        'YLim',[0.5,n_rows+0.5],...
        'DataAspectRatio',[1 1 1],...
        'PlotBoxAspectRatio',[n_cols n_rows 1]);
    set(image_axes_h,'CLim',[0 30])
    temp = ch1;
    image_h = image('Parent',image_axes_h,...
        'Tag','image_h',...
        'CData',temp,...
        'SelectionHighlight','off',...
        'EraseMode','none',...
        'CDataMapping','scaled');
    
    
    ellipse_h=...
        line('Parent',image_axes_h,...
        'Color','r',...
        'Tag','border_h',...
        'XData',roi_borders{count}(:,1),...
        'YData',roi_borders{count}(:,2),...
        'ZData',repmat(2,size(roi_borders{count}(:,1))));
    com=border_com(roi_borders{count}');
    %     this_label_h(count)=...
    %         text('Parent',image_axes_h,...
    %         'Position',[com(1) com(2) 2],...
    %         'String',thislabel,...
    %         'HorizontalAlignment','center',...
    %         'VerticalAlignment','middle',...
    %         'Color','r',...
    %         'Tag','label_h',...
    %         'Clipping','on');
    
    [q,w] = sort(str2num(char(stim_struct.condnames)));
    if isempty(w)
        q = condnames;
        w = [1:length(condnames)];
    end
    
    subplot(double(nconds)+1,2,[1 2],'Parent',summary_fig_h);
    raw_optical_lines_h=line('Parent',gca,...
        'XData',dt*[0:size(roi_sums_all,1)-1]',...
        'YData',roi_sums_all(:,count),...
        'Visible','on',...
        'Color',[0.5 0.5 0.5],...
        'Visible','on');
    
    for y = 1:nconds
        thisy = w(y);
        subplot(double(nconds)+1,2,double(y)*2+2,'Parent',summary_fig_h);
        this_stimframes = allstimframes(thisy:nconds:end);
        for temp = 2:length(this_stimframes)
            avg_roi_sums(temp-1,:) = roi_sums_all(this_stimframes(temp):this_stimframes(temp)+t_pre_frames+t_stim_frames+t_post_frames-1,count);
            avg_optical_lines_h(temp)=line('Parent',gca,...
                'XData',dt*[0:(t_pre_frames+t_stim_frames+t_post_frames-1)],...
                'YData',avg_roi_sums(temp-1,:),...
                'Visible','on',...
                'Color',[0.5 0.5 0.5],...
                'Visible','on');
        end
        avg_optical_lines_h(1)=line('Parent',gca,...
            'XData',dt*[0:(t_pre_frames+t_stim_frames+t_post_frames-1)],...
            'YData',mean(avg_roi_sums),...
            'Visible','on',...
            'Color','r',...
            'LineWidth',2,...
            'Visible','on');
        thismean = mean(avg_roi_sums);
        roi_sum_area(y) = sum(thismean-repmat(mean(thismean(1:5)),[length(thismean) 1])');
        
        ymax = max([ymax max(mean(avg_roi_sums))]);
        ymin = min([ymin min(mean(avg_roi_sums))]);
        ylabel(q(y));
    end
    ymax = ymax+2;
    ymin = ymin-2;
    subplot(double(nconds)+1,2,double(y)*2+2,'Parent',summary_fig_h);
    xlabel('sec.');
    
    %     xlim = get(avg_optical_axes_h,'XLim');
    subplot(double(nconds)+1,2,[1 2],'Parent',summary_fig_h);
    set(gca,'YLim',[ymin ymax]);
    title(sprintf('cell %s, #%d',thislabel,count))
    for y = 1:nconds
        thisy = w(y);
        subplot(double(nconds)+1,2,double(y)*2+2,'Parent',summary_fig_h);
        set(gca,'YLim',[ymin ymax]);
        %         set(gca,'XLim',xlim);
    end
    %
    roi_sum_area = roi_sum_area/(max(roi_sum_area))
    subplot(double(nconds)+1,2,[11 13 15],'Parent',summary_fig_h);
    thisscale = [0 0.5 1];
    tout = [pi/4 pi/2 3*pi/4 pi 5*pi/4 3*pi/2 7*pi/4 2*pi pi/4]';
    rout = [roi_sum_area' ;roi_sum_area(1)];
    polar2(tout,rout*thisscale(3),thisscale,[1 0 0]);
    hold on;
    hline = findobj(gca,'Type','line');
    set(hline,'LineWidth',2)
    
    thisname = sprintf('e2198_#%0.3d_cell_%s',count,thislabel);
    print(gcf, '-r300', sprintf('%s\\%s.png',figdir,thisname), '-dpng');
    
    close all
    
end
   

%%
% Example of reshaping matrix 
for count = 1:n_reps
    for count2 = 1:nconds
        for count3 = 1:n_rois
                count4 = cond_order(count2);
                thiscond = roi_sums_all(allstimframes((count-1)*nconds+count4):allstimframes((count-1)*nconds+count4)+t_trial-1,count3)';
                allrois(count3,(count2-1)*t_trial+1:(count2-1)*t_trial+t_trial,count)=thiscond;
            
        end
    end
end

% Due to adaptation seen in Off RGCs for first trial, only use last 4
% trials
allrois_avgreps_selectedreps = mean(allrois(:,:,4:5),3);
figure; imagesc(allrois_avgreps_selectedreps)

%%% Type statistics
% This was a very preliminary analysis...
cells_nr = find(roi_struct.type==1);
cells_off = find(roi_struct.type==2);
cells_on = find(roi_struct.type==4);
cells_onoff = [find(roi_struct.type==6) ; find(roi_struct.type==7)];
cells_ds = find(roi_struct.type==6);
cells_other =  find(roi_struct.type==8);

ds1 = [1 5 17 25]; % peak at 2 use 'r'
ds2 = [9 11 12 20 24 28];  % peak at 3 4 use 'g'
ds3 = [3 6 8 13 16 19 27 30]; % peak 5 6 use 'm'
ds4 = [2 7 10 15 18 23 32]; % peak at  7 8 1 use [1 0.7 0]

other = [14 26 4 29 31]

figure; imagesc(allrois_avgreps_selectedreps(cells_nr,:))
figure; imagesc(allrois_avgreps_selectedreps(cells_ds,:))
figure; imagesc(allrois_avgreps_selectedreps(cells_on,:))
figure; imagesc(allrois_avgreps_selectedreps(cells_off,:))
figure; imagesc(allrois_avgreps_selectedreps(cells_onoff,:))
figure; imagesc(allrois_avgreps_selectedreps(cells_other,:))

