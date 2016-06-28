function plot_grouped_ca(cell_info, cell_dict_j, roi_sums_means_flatten, types, type_styles, subplots) %plot_average, plot_individual)

if ~iscell(types) %ischar(types)
    types = {types};
end
ntypes = length(types);

%plot_average = &&
%plot_individual
%subplots = plot_individual
subplots = 0;

cla;
legends = {};
for k = 1:ntypes
    ctype = types{k};
    if iscell(ctype)
        legends{end+1} = strjoin(ctype);
    else  % char
        legends{end+1} = ctype;
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
    mean_ca_norm = (mean_ca - min(mean_ca)) / (max(mean_ca) - min(mean_ca));

    if subplots
        subplot(ntypes, 1, k);
        plot(ca);
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
    plot(mean_ca_norm, style, 'LineWidth', linewidth);
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
    ax.YTick = [];

    if subplots
    xticklabel = ax.XTickLabel;
    ax.XTickLabel = [];
    ylabel(ctype, 'Rotation', 0, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    end
end

if subplots
ax.XTickLabel = xticklabel;
else
    legend(legends);
end

%legend off;

ax.XTick = [0:4]/0.128;
%xlim([0 4]/0.128);
ylim([0 1.1])
ax.XTickLabel = num2str([0:4].');
%ax.XTickLabel = {'1.0', '2.0', '4.0'}; %num2str(ax.XTick.' * 0.128);
grid off
xlabel 'time (s)'
%{
ylabel('$\Delta F \over F$', 'Interpreter', 'LaTex', ...
    'FontSize', 14, 'Rotation', 0, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
%}

title 'Calcium'


%{ 
types = {'1wt'  '4ow' '6sw' '8w'};
types = {'1wt'  '1no' '1ni' '2an' '2aw' '2o' '2i' '3o' '3i' '4on' '4ow' '4i'};
'5so' '5si' '5to' '5ti'}
types = { '6sn' '6sw' '6t' '7o' '7i' '8w' '8n' '9w' '9n'};

figure;plot_grouped_ca(cell_info, cell_dict_j, roi_sums_means_flatten, types)
%}
