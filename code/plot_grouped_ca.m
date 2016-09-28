function plot_grouped_ca(cell_info, cell_dict_j, roi_sums_means_flatten, types, type_styles, subplots, normalization, y_max, draw_hat, legendtext) %plot_average, plot_individual)

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

if ~iscell(types) %ischar(types)
    types = {types};
end
ntypes = length(types);

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
ax.ColorOrder(8,:) = ax.ColorOrder(5,:) * 0.6; %dark green. 0.3* 0.4660    0.6740    0.1880
if ntypes>8
    %colormap(distinguishable_colors(ntypes));
    ax.ColorOrder = distinguishable_colors(ntypes);
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

    cells = get_cell_info(cell_info, ctype);
    idx = ismember(cell_dict_j(:,2), [cells.cell_id]);
    ca_ids = cell_dict_j(idx,1);
    ca = roi_sums_means_flatten(:,ca_ids);
    % average across trials and directions
    ca = squeeze(mean(reshape(ca, 31, 8, []), 2));
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
    elseif plot_individual
        chartlines = plot(ca, 'Color', kolor);
        ax.ColorOrderIndex = ColorOrderIndex;
        if ~plot_average
            ax.ColorOrderIndex = ColorOrderIndex+1;
        end
    elseif plot_stderr
        alphav = 0.5;
        alphav = 0.3;
        bounds = [mean_ca_norm+ca_stderr mean_ca_norm-ca_stderr];
        %chartlines = plot(bounds, 'Color', kolor);
        fill([1:31 31:-1:1].', [mean_ca_norm+ca_stderr; flip(mean_ca_norm-ca_stderr)], kolor, 'FaceAlpha', alphav, 'LineStyle', 'none');%, 'EdgeAlpha', 0);
        ax.Box = 'off'; % somehow "fill" auto sets box to on
        hold on;
        %plot([mean_ca_norm+ca_std mean_ca_norm-ca_std], ':','Color', kolor);
        %plot([mean_ca_norm+2*ca_stderr mean_ca_norm-2*ca_stderr], ':','Color', kolor);
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
    else
        linewidth = 2;
    end

    if exist('type_styles', 'var') && k <= length(type_styles)
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
    if exist('legendtext', 'var')
        legends = legendtext;
    end
    legend(legendentries, legends);
end

%legend off;

ax.XTick = [0:4]/0.128;
xlim([0 4]/0.128);
%ylim([0 1.1])
ax.XTickLabel = num2str([0:4].');
%ax.XTickLabel = {'1.0', '2.0', '4.0'}; %num2str(ax.XTick.' * 0.128);
grid off
xlabel 'time (s)'

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
        line([7 31], [1 1], 'Color', 0.7*[1 1 1], 'LineWidth', 0.2);
    case 'on'
        line([7 15], [1 1], 'Color', 0.7*[1 1 1], 'LineWidth', 0.2);
    case 'off'
        line([15 31], [1 1], 'Color', 0.7*[1 1 1], 'LineWidth', 0.05);
    case 'none'
    otherwise
        error 'invalid arg'
    end
end

title(['n_cells = ', num2str(ncells)], 'Interpreter', 'none');


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
