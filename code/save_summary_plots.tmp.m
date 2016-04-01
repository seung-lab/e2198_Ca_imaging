
%%
% Make summary figure for each cell


coeffs16_reshape = reshape(coeffs16{1,2}(:,3:end-16).', 2, 8, 634);
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



%for count = 1:n_rois
for count = 1
    ind = count;
    
    %close all;

    %% group into folders
    omni_id = cell_dict(cell_dict(:,2)==ind, 1);
    if ~isempty(omni_id) %find(cell_dict(:,2)==ind)
        celltype = cell_info([cell_info.cell_id]==omni_id).type;
        subdir = celltype;
        folder = [figdir subdir];
        if ~exist(folder, 'dir')
            mkdir(folder);
        end
    else
        folder = figdir;
    end

    %%-- figure()
    %thislabel = show_id_to_alpha(count)
    thislabel = num2str(count)
    summary_fig_h = ...
        figure('Tag','summary_fig_h',...
        'Name',sprintf('cell_%s',thislabel),...
        'NumberTitle','off',...
        'PaperPositionMode','auto',...
        'InvertHardcopy','off',...
        'DoubleBuffer','on',...
        'Position',[341 79 465 1019]);
    n_rows_subplot = double(nconds)+2;
    
    ymax = 0;
    ymin = 0;
    %%-- Ca overview image
    subplot(n_rows_subplot,2,[3 5 7],'Parent',summary_fig_h);
    colormap('gray(256)');
    temp = ca_overview;
    n_cols = size(temp,1); n_rows = size(temp,2);
    %%{
    imshow(temp);
    %image(temp, ...
    %    'CDataMapping','scaled');
    image_axes_h = gca();
    set(image_axes_h,'CLim',[0 30])
    %}
    %image_axes_h = gca();
    %{
    %image_axes_h = axes(... %'Parent',summary_fig_h,...
    set(image_axes_h,...    'Tag','image_axes_h',...
        'YDir','reverse',...
        'DrawMode','normal',...
        'Visible','on',...
        'Units','Normalized',...
        ... %'Position',[0.1 0.50 0.35 0.35],...
        'XLim',[0.5,n_cols+0.5],...
        'YLim',[0.5,n_rows+0.5],...
        'DataAspectRatio',[1 1 1],...
        'PlotBoxAspectRatio',[n_cols n_rows 1]);
    %}
    %{
    set(image_axes_h,'CLim',[0 30])
    image_h = image('Parent',image_axes_h,...
        'Tag','image_h',...
        'CData',temp,...
        'SelectionHighlight','off',...
        'EraseMode','none',...
        'CDataMapping','scaled');
    %}
    
    %%{
    ellipse_h=...
        line('Parent',image_axes_h,...
        'Color','r',...
        'Tag','border_h',...
        'XData',roi_borders{count}(:,1),...
        'YData',roi_borders{count}(:,2),...
        'ZData',repmat(2,size(roi_borders{count}(:,1))));
    %}
    %com=border_com(roi_borders{count}');
    com=mean(roi_borders{count}', 2);
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
    %%{
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
            'Color','r',...
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
    %}
    
    %%{
    %%-- polar plot
    subplot(n_rows_subplot,2,[11 13 15],'Parent',summary_fig_h);
    polar_tuning(coeffs16_reshape(:,:,ind), order);
    

    %%{
    %%-- plot curve fit
    subplot(n_rows_subplot,2,[19 20],'Parent',summary_fig_h);
    %plot([roi_sums_means_flatten(:,ind) fit16{1}(:,ind) fit16{3}(:,ind)])
    plot([roi_sums_means_flatten(:,ind) fit16{1,2}(:,ind)]);hold on;plot([fit16{1,1}(:,ind) fit16{3,2}(:,ind)], '-.')
    ax = gca(); ax.XTick = t1t2_to_ti(9,16); grid on;
    if diff(ylim())<8
        lim = ylim;
        if lim(2)<6
            lim(2) = 6;
        end
        ylim(lim);
    end
    %}


    thisname = sprintf('e2198_#%0.3d_cell_%s',count,thislabel);
    %{
    print(gcf, '-r300', sprintf('%s\\%s.png',figdir,thisname), '-dpng');
    close;
    %}
    
end
%}
