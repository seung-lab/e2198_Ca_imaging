
get_super_types
%figure;plot_grouped_ca(cell_info, cell_dict_j, roi_sums_means_flatten, ca_groups.types, [],[],[], 1.1);
%legend(ca_groups.name);
figure;plot_grouped_ca(cell_info, cell_dict_j, roi_sums_means_flatten, ca_groups.types, [],4,[], 1.1, [], ca_groups.name);
ax = gca();
ax.FontSize = 13;

figure;plot_grouped_ca(cell_info, cell_dict_j, roi_sums_means_flatten, ca_groups3.types, [],4,[], 1.1, [], ca_groups3.name, [], ca_groups3.shift);
ax = gca();
ax.FontSize = 13;
%return

%figure;plot_grouped_ca(cell_info, cell_dict_j, roi_sums_means_flatten, types_alpha, [],[],[], 1.1);
figure;plot_grouped_ca(cell_info, cell_dict_j, roi_sums_means_flatten, types_alpha, [],4,[], 1.1);
ax = gca();
ax.FontSize = 13;

figure;plot_grouped_ca(cell_info, cell_dict_j, roi_sums_means_flatten, types_oodsgc, [],4,[], 1.1);
ax = gca();
ax.FontSize = 13;

% For axis only
figure;plot_grouped_ca(cell_info, cell_dict_j, roi_sums_means_flatten, '37', [],[],[], 1.1);
ax = gca();
ax.FontSize = 13;
title('Calcium');
ylabel('$\Delta F \over F$', 'Interpreter', 'LaTex', ...
    'FontSize', 16, 'Rotation', 0, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
xlim([0,4]/0.128)
f = gcf();
f.Position(3:4) = [570,422];

figure;plot_grouped_ca(cell_info, cell_dict_j, roi_sums_means_flatten, {'1ni' '1no'}, [],4,[], 1.1, []);
ax = gca();
ax.FontSize = 13;