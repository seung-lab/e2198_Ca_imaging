
%%
% Make summary figure for each cell


% old single exponential fit
%coeffs16_reshape = reshape(coeffs16{1,2}(:,3:end-16).', 2, 8, 634);
%[~, tuning_onoff] = tuning_from_fit(coeffs16{1,2});
%tau_offset_list = coeffs16{1,2}(:,1:2);
% exp-exp fit
[tuning, tuning_onoff] = tuning_from_fit(coeffs16{3,2});
tau_offset_list = coeffs16{3,2}(:,1:3);

%coeffs16_reshape = reshape(coeffs16{1}(:,3:end).', 2, 8, 634);
%coeffs_ordered = coeffs_reshape(:, order, :);
%vert_offsets = coeffs16{1}(:,1).';

figdir = '~/dev/e2198_Ca_imaging/cell_summary/';
%figdir = '/omniData/e2198_reconstruction/e2198_Ca_imaging/new_cell_summary/'
mkdir(figdir);
%{
addpath('/usr/people/smu/seungmount/research/jinseopk/e2198/bin/analysis/')
gc = cell_info_typedef_gc();
cell_info=cell_info_set_type();
%}

OVimg_fn = 'AVG_043_p2_g7_ZoomedOut_center.tif';
ca_overview = imread(OVimg_fn);


%{
rawdata_fn = ['034_p2_g7_DSbars_200um.cfd';
           '035_p2_g7_DSbars_200um.cfd';
           '036_p2_g7_DSbars_200um.cfd';
           '037_p2_g7_DSbars_200um.cfd';
           '038_p2_g7_DSbars_200um.cfd';
           '039_p2_g7_DSbars_200um.cfd';
           '040_p2_g7_DSbars_200um.cfd';
           '041_p2_g7_DSbars_200um.cfd';
           '042_p2_g7_DSbars_200um.cfd'];
       
% Read raw CFD file
[dx,dy,channels,frames,ch1,ch2,t_scan,t_frame,ch1_t_stim] = read_cfd(['../rawdata/' rawdata_fn(1,:)]);
stim_fn = 'sDS 8x45deg, narrow, 4.0s.QDS';
stim_struct = load_stim(stim_fn);
%}
data_fn = ['035_p2_g7_DSbars_200um_export.mat';
           '036_p2_g7_DSbars_200um_export.mat';
           '037_p2_g7_DSbars_200um_export.mat';
           '038_p2_g7_DSbars_200um_export.mat';
           '039_p2_g7_DSbars_200um_export.mat';
           '040_p2_g7_DSbars_200um_export.mat';
           '041_p2_g7_DSbars_200um_export.mat';
           '042_p2_g7_DSbars_200um_export.mat';
           '034_p2_g7_DSbars_200um_export.mat'];
load(data_fn(1,:)); %for t_frame

[n_samples,t_min,t_max,T,dt,fs] = time_info(t_frame);
t_pre_frames = round(stim_struct(1).t_pre/dt);
t_stim_frames = round(stim_struct(1).t_stim/dt);
t_post_frames = round(stim_struct(1).t_post/dt);
% nconds: the 8 direction conditions
%nconds = 8; t_trial=31; n_reps = 5;


dt = 0.128;
t_trial=31; % 31 frames per trial
allstimframes = [1:t_trial:size(roi_sums_all,1)];

aa = [53 133 161 363 519 ...
297 379 411 440 533 601 ...
057 241 596];


visible = 0


for count = 1:n_rois
%for count = 79
%for count = 53
    ind = count;
    
    %close all;

    %% group into folders
    omni_id = cell_dict(cell_dict(:,2)==ind, 1);
    omni_id = cell_dict_j(cell_dict_j(:,1)==ind, 2);
    if ~isempty(omni_id) %find(cell_dict(:,2)==ind)
        if length(omni_id)>1    %TODO hack
            warning(horzcat('dup ', num2str(ind), ' ', num2str(omni_id.')));
            omni_id = omni_id(1);
        end
        celltype = cell_info([cell_info.cell_id]==omni_id).type;
        if isempty(celltype)
            classes = {'cantknow', 'GC', 'AC', 'BC', 'uncertain'};  % class: cantknow=0, GC=1, AC=2, BC=3, uncertain=4
            celltype = cell_info([cell_info.cell_id]==omni_id).class;
            celltype = classes{celltype+1};
        end
        subdir = celltype;
        folder = [figdir subdir];
    else
        folder = [figdir 'no_match'];
    end
    if ~exist(folder, 'dir')
        mkdir(folder);
    end

    %%-- figure()
    %thislabel = show_id_to_alpha(count)
    %thislabel = num2str([count omni_id]);
    thislabel = num2str(omni_id);
    summary_fig_h = ...
        figure('Tag','summary_fig_h',...
        'Name',sprintf('cell_%s',thislabel),...
        'NumberTitle','off',...
        'PaperPositionMode','auto',...
        'InvertHardcopy','off',...
        'DoubleBuffer','on',...
        'Position',[100 79 600 919]);
    if ~visible
        summary_fig_h.Visible = 'off';
    end
    n_rows_subplot = double(nconds)+3;
    
    ymax = 0;
    ymin = 0;

    %%-- Ca overview image
    subplot(n_rows_subplot,2,[3 5 7],'Parent',summary_fig_h);
    %imshow(ca_overview);
    imshow(ca_overview.');  % in final orientation for paper
    image_axes_h = gca();
    set(image_axes_h,'CLim',[0 30])
    
    ellipse_h=...
        line('Parent',image_axes_h,...
        'Color','r',...
        'LineWidth',0.8,...
        'Tag','border_h',...
        'YData',roi_borders{count}(:,1),...  % X/Y flipped for final orientation for paper
        'XData',roi_borders{count}(:,2),...
        'ZData',repmat(2,size(roi_borders{count}(:,1))));
    
    %[q,w] = sort(str2num(char(stim_struct.condnames)));
    angles_double = str2num(char(angles));
    angles_double = angles_double + 90;  % convert to final coord for paper
    angles_double(angles_double>360) = angles_double(angles_double>360) - 360;
    [ordered, order] = sort(angles_double);
    [q,w] = sort(angles_double);
    %q = q([3:8 1:2])
    if isempty(w)
        q = condnames;
        w = [1:length(condnames)];
    end
    
    %%-- full time course responses 
    subplot(n_rows_subplot,2,[1 2],'Parent',summary_fig_h);
    raw_optical_lines_h=line('Parent',gca,...
        'XData',dt*[0:size(roi_sums_all,1)-1]',...
        'YData',roi_sums_all(:,count),...
        'Visible','on',...
        'Color',[0.5 0.5 0.5],...
        'Visible','on');
    xlim([0 dt*(size(roi_sums_all,1)-1)]);
    
    %%-- response for each condition(direction)
    for y = 1:nconds    %from 45 deg to 360
        thisy = w(y);
        subplot(n_rows_subplot,2,double(y)*2+2,'Parent',summary_fig_h);
        this_stimframes = allstimframes(thisy:nconds:end);
        for temp = 2:length(this_stimframes)    % 2:5
            % avg_roi_sums (4x31 array) is not average yet, thismean is.
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
            'Color',[ 0    0.4470    0.7410],...
            'LineWidth',2,...
            'Visible','on');
        
        ymax = max([ymax max(mean(avg_roi_sums))]);
        ymin = min([ymin min(mean(avg_roi_sums))]);
        ylabel(q(y));
    end
    ymax = ymax+2;
    ymin = ymin-2;
    subplot(n_rows_subplot,2,double(y)*2+2,'Parent',summary_fig_h);
    xlabel('sec.');
    

    %%-- title and global ylim
    %     xlim = get(avg_optical_axes_h,'XLim');
    subplot(n_rows_subplot,2,[1 2],'Parent',summary_fig_h);
    set(gca,'YLim',[ymin ymax]);
    title(sprintf('cell %s, #%d',thislabel,count))
    for y = 1:nconds
        thisy = w(y);
        subplot(n_rows_subplot,2,double(y)*2+2,'Parent',summary_fig_h);
        set(gca,'YLim',[ymin ymax]);
        %         set(gca,'XLim',xlim);
    end
    %

    %%-- polar plot
    subplot(n_rows_subplot,2,[11 13 15],'Parent',summary_fig_h);
    %polar_tuning(tuning_onoff(:,:,ind), order);
    polar_tuning(tuning(:,:,ind), order);
    

    %%{
    %%-- plot curve fit
    subplot(n_rows_subplot,2,[19 20 21 22],'Parent',summary_fig_h);
    %plot([roi_sums_means_flatten(:,ind) fit16{1}(:,ind) fit16{3}(:,ind)])
    % Both single exp and double exp fits:
    %plot([roi_sums_means_flatten(:,ind) fit16{1,2}(:,ind)]);hold on;plot([fit16{3,2}(:,ind)], '-.')
    % Double exp fit only:
    plot([roi_sums_means_flatten(:,ind) fit16{3,2}(:,ind)]);
    ax = gca(); ax.XTick = t1t2_to_ti(8.8215,16.625); grid on;  %  [8.8215,16.625] = 1+[1 2]/0.128
    ax.XTickLabels = sort([1:4:32 2:4:32]);
    if diff(ylim())<8
        lim = ylim;
        if lim(2)<6
            lim(2) = 6;
        end
        ylim(lim);
    end
    title(sprintf('tau = %.1f  %.1f, offset = %.1f', deal(tau_offset_list(ind, :))))
    %}


    thisname = sprintf('e2198_#%0.3d_cell_%d',count,omni_id);
    %%{
    print(gcf, '-r200', sprintf('%s\\%s.png',folder,thisname), '-dpng');
    close(summary_fig_h);
    %}
    
end
%}

%{
tau s
done: tuning curve for type
stratification
%}
