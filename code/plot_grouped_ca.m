function plot_grouped_ca(cell_info, cell_dict_j, roi_sums_means_flatten, types, type_styles, subplots, normalization, ...
     y_max, draw_hat, legendtext, value_to_scale_to_1, shift, colororder) %plot_average, plot_individual)

% cell_dict_j: cell_dict_j, or the calcium ids as a vector.

draw_stimulus_shade = 1;

if ~exist('normalization', 'var') || isempty(normalization)
    normalization = 'full'; % on / off / full / none
end
if ~exist('y_max', 'var')
    y_max = 'max'; % value / 'max'
end
if ~exist('draw_hat', 'var')
    if strcmp(normalization, 'on') || strcmp(normalization, 'off')
        draw_hat = 1;
    else
        draw_hat = 0;
    end
end

if ~exist('shift', 'var')
    shift = zeros(size(types));
end

if ~iscell(types) %ischar(types)
    types = {types};
end
ntypes = length(types);

[nframes, n_rois] = size(roi_sums_means_flatten);
if exist('value_to_scale_to_1', 'var') && ~isempty(value_to_scale_to_1)
    %value_to_scale_to_1 = value_to_scale_to_1(ca_ids);

    % average across conditions
    min_of_mean = squeeze(min(mean(reshape(roi_sums_means_flatten, 31, 8, []), 2)));

    roi_sums_means_flatten = roi_sums_means_flatten - repmat(min_of_mean.', nframes, 1);
    roi_sums_means_flatten = roi_sums_means_flatten ./ repmat(value_to_scale_to_1, nframes, 1);
end

if subplots == 1
    subplots = 1;
    plot_individual = 1;
    plot_average = 1;
    plot_stderr = 0;
elseif subplots == 2
    subplots = 0;
    plot_individual = 1;
    plot_average = 1;
    plot_stderr = 0;
elseif subplots == 3
    subplots = 0;
    plot_individual = 1;
    plot_average = 0;
    plot_stderr = 0;
elseif subplots == 4
    subplots = 0;
    plot_individual = 0;
    plot_average = 1;
    plot_stderr = 1;
elseif subplots == 5    % same as 4 but plotting standard deviation instead of s.e.m.
    subplots = 0;
    plot_individual = 0;
    plot_average = 1;
    plot_stderr = 2;    % standard deviation instead of s.e.m.
elseif subplots == 6    % same as 4 but plotting individuals too
    subplots = 0;
    plot_individual = 1;
    plot_average = 1;
    plot_stderr = 1;
else
    subplots = 0;
    plot_individual = 0;
    plot_average = 1;
    plot_stderr = 0;
end
%plot_average = &&
%plot_individual
%subplots = plot_individual
%if ntypes
%subplots = 0;
if subplots
plot_individual = 1;
end

cla;
ax = gca();
if draw_stimulus_shade  % TODO: not implemented for subplots
    shadecolor = 0.93 * [1 1 1];
    %plot([], []);  % sets up the default without
    fill([0 1 1 0].'/0.128, 150*[1 1 0 0], shadecolor, 'LineStyle', 'none'); %'FaceAlpha', alphav, 
    hold on
    fill([2 4 4 2].'/0.128, 150*[1 1 0 0], shadecolor, 'LineStyle', 'none');
    ax.Box = 'off'; % somehow "fill" auto sets box to on
    ax.ColorOrderIndex = 1;
end
ax.ColorOrder(8,:) = ax.ColorOrder(5,:) * 0.6; %dark green. 0.3* 0.4660    0.6740    0.1880
if ntypes>8
    %colormap(distinguishable_colors(ntypes));
    ax.ColorOrder = distinguishable_colors(ntypes);
end
if exist('type_styles', 'var') && isnumeric(type_styles) && ~isempty(type_styles)  % treat as color order
    ax.ColorOrder = ax.ColorOrder(type_styles,:);
end
legends = {};
legendentries = [];
ymax = 0;
ncells = [];
for k = 1:ntypes
    ctype = types{k};
    if iscellstr(ctype)
        legends{end+1} = strjoin(ctype);
    elseif ischar(ctype)  % char
        legends{end+1} = ctype;
    elseif isnumeric(ctype)
        legends{end+1} = num2str(ctype);
    end

    if isvector(cell_dict_j)
        ca_ids = cell_dict_j{k};  % this argument is taken to be calcium ids already.
    else
    cells = get_cell_info(cell_info, ctype);
    idx = ismember(cell_dict_j(:,2), [cells.cell_id]);
    ca_ids = cell_dict_j(idx,1);
    end
    ca = roi_sums_means_flatten(:,ca_ids);
    if exist('value_to_scale_to_1', 'var')
        ca = ca(:, isfinite(ca(1,:)));
    end
    % average across trials and directions
    %if size(ca,1)~=31
    ca = squeeze(mean(reshape(ca, 31, 8, []), 2));
    %end

    % make first frame as 0s. Maybe I should make this -1 shift as the default "shift" argument?
    %ca = [ca; ca(1,:)];
    ca = circshift(ca, -1, 1);

    % average across cells
    mean_ca = mean(ca, 2);
    % "normalization"
    switch normalization
    case 'full'
        mean_ca_norm = (mean_ca - min(mean_ca)) / (max(mean_ca) - min(mean_ca));
        ca = (ca - repmat(min(ca),31,1)) ./ repmat((max(ca) - min(ca)) ,31,1);
    case 'on'
        mean_ca_norm = (mean_ca - min(mean_ca)) / (max(mean_ca(7:15)) - min(mean_ca));
        ca = (ca - repmat(min(ca),31,1)) ./ repmat((max(ca(7:15,:)) - min(ca)) ,31,1);
    case 'off'
        mean_ca_norm = (mean_ca - min(mean_ca)) / (max(mean_ca(16:end)) - min(mean_ca));
        ca = (ca - repmat(min(ca),31,1)) ./ repmat((max(ca(16:end,:)) - min(ca)) ,31,1);
        %mean_ca_norm = (mean_ca - min(mean_ca(12:20)));
        %mean_ca_norm = mean_ca_norm / max(mean_ca_norm(16:end));
    case 'none'
        mean_ca_norm = mean_ca;
    otherwise
        error 'invalid arg'
    end
    if exist('value_to_scale_to_1', 'var')
        mean_ca_norm = mean_ca;  % no double normalization
    end

    mean_ca_norm = circshift(mean_ca_norm, shift(k));
    ca = circshift(ca, shift(k));

    ncells = [ncells size(ca,2)];
    ca_std = std(ca, [], 2);
    ca_stderr = ca_std ./ sqrt(size(ca,2));
    if plot_stderr  % mean of normalized curves instead
        mean_ca_norm = mean(ca, 2);
        scaling = (max(mean_ca_norm) - min(mean_ca_norm));
        mean_ca_norm = (mean_ca_norm - min(mean_ca_norm)) / scaling;
        ca_std = ca_std / scaling;
        ca_stderr = ca_stderr / scaling;
    end

    ymax = max(max(mean_ca_norm), ymax);

        ax = gca();
        ColorOrderIndex = ax.ColorOrderIndex;
        kolor = ax.ColorOrder(mod(ax.ColorOrderIndex-1,length(ax.ColorOrder)) + 1, :);
    if subplots && plot_individual
        subplot(ntypes, 1, k);
        plot(ca);
    else
        if plot_stderr
        alphav = 0.5;
        alphav = 0.3;
        bounds = [mean_ca_norm+ca_stderr mean_ca_norm-ca_stderr];
        %chartlines = plot(bounds, 'Color', kolor);
        if plot_stderr==1 % standard err
            to_fill = [mean_ca_norm+ca_stderr; flip(mean_ca_norm-ca_stderr)];
        elseif plot_stderr==2 % standard deviation
            to_fill = [mean_ca_norm+ca_std; flip(mean_ca_norm-ca_std)];
        else
            error
        end
        fill([1:31 31:-1:1].', to_fill, kolor, 'FaceAlpha', alphav, 'LineStyle', 'none');%, 'EdgeAlpha', 0);
        ax.Box = 'off'; % somehow "fill" auto sets box to on
        hold on;
        %plot([mean_ca_norm+ca_std mean_ca_norm-ca_std], ':','Color', kolor);
        %plot([mean_ca_norm+2*ca_stderr mean_ca_norm-2*ca_stderr], ':','Color', kolor);
        end
        if plot_individual
            chartlines = plot(ca, 'Color', kolor);
        end
        ax.ColorOrderIndex = ColorOrderIndex;
        if ~plot_average
            ax.ColorOrderIndex = ColorOrderIndex+1;
        end
    end
    hold on;
    %plot(mean_ca, '.k', 'MarkerSize', 12);
    %plot(mean_ca, '.-', 'LineWidth', 1);
    %plot(mean_ca, '.-', 'LineWidth', 1, 'MarkerSize', 10);
    %plot(mean_ca, 'LineWidth', 3);
    if subplots
        linewidth = 3;
    elseif 1 && 0.2 > ax.Position(4) - ax.Position(2)
        linewidth = 1;
    else
        linewidth = 2;
    end

    if exist('type_styles', 'var') && ~isnumeric(type_styles) && k <= length(type_styles)
        style = type_styles{k};
    else
        style = '';
    end
    if plot_average
        chartlines = plot(mean_ca_norm, style, 'LineWidth', linewidth);
    end
    legendentries(end+1) = chartlines(end);
    %{
    if k > 2
        plot(mean_ca_norm, 'LineWidth', 3);
    else
        plot(circshift(mean_ca_norm, -7), 'LineWidth', 3);
    end
    %}
    ax = gca();
    %ax.Color = 'none';
    %ax.Visible = 'off';
    %xlim([0 100]);
    ax.XTick = t1t2_to_ti(8,16); grid on;
    ax.XTick = [0:4]/0.128;
    if ~strcmp(normalization, 'none')
        ax.YTick = [];
    end

    if subplots
    xticklabel = ax.XTickLabel;
    ax.XTickLabel = [];
    ylabel(ctype, 'Rotation', 0, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    end
end

if subplots
ax.XTickLabel = xticklabel;
else
end

%legend off;

ax.XTick = [0:4]/0.128;
xlim([0 4]/0.128);
%ylim([0 1.1])
ax.XTickLabel = num2str([0:4].');
%ax.XTickLabel = {'1.0', '2.0', '4.0'}; %num2str(ax.XTick.' * 0.128);
grid off
xlabel 'Time (s)'

if ~strcmp(normalization, 'none')
    %ymax = ylim();
    %ymax = ymax(2);
    %if ymax>1
    %    ylim([0 2.1])
    %end
    ylim([0 ymax]);
    ylim([0 2.2]);
    %ax.YGrid = 'on';
    %ax.YTick = [0 1];
    %ax.YTickLabel = [];
    if strcmp(y_max, 'max')
        y_max = ymax;
    end
    ylim([0 y_max]);
end

if draw_hat
    switch normalization
    case 'full'
        %line([0 31], [1 1], 'Color', 0.7*[1 1 1], 'LineWidth', 0.2);
        line([7 31], [1 1], 'Color', 0.7*[1 1 1], 'LineWidth', 0.4);
    case 'on'
        line([7 15]-1, [1 1], 'Color', 0.7*[1 1 1], 'LineWidth', 0.4);  % "-1": account for the one frame shift in making first frame as 0s.
    case 'off'
        line([15 31], [1 1], 'Color', 0.7*[1 1 1], 'LineWidth', 0.4);
    case 'none'
    otherwise
        error 'invalid arg'
    end
end

title(['n_cells = ', num2str(ncells)], 'Interpreter', 'none');

% ah for 2017a had to move this section to the end in order not to show an entry for the hat...
if subplots
else
    if exist('legendtext', 'var')
        legends = legendtext;
    end
    [legh,objh] = legend(legendentries, legends);
    set(objh,'linewidth',2);
end

%{
ylabel('$\Delta F \over F$', 'Interpreter', 'LaTex', ...
    'FontSize', 14, 'Rotation', 0, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
%}

%title 'Calcium'


%{ 
types = {'1wt'  '4ow' '6sw' '8w'};
types = {'1wt'  '1no' '1ni' '2an' '2aw' '2o' '2i' '3o' '3i' '4on' '4ow' '4i'};
'5so' '5si' '5to' '5ti'}
types = { '6sn' '6sw' '6t' '7o' '7i' '8w' '8n' '9w' '9n'};

figure;plot_grouped_ca(cell_info, cell_dict_j, roi_sums_means_flatten, types)
%}
