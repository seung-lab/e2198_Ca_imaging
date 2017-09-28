%generate_print.m

%cell_info_polarplot_pref_dir2(cell_info,ca_dsos, [0 1 2])


cell_info_polarplot_pref_dir2(cell_info,ca_dsos, [1 2 0])

cell_info_polarplot_pref_dir(cell_info,ca_dsos, [], 0, 0, 4)

% fig 6
fig4_ca
typefit = fit_type_means(roi_sums_xcond_typemeans);
typefit_responsive = fit_type_means(roi_sums_xcond_responsive_typemeans(:,'trace')); 
% really poor fits: 7iv, (1wt) hmmm seems like I should have used a min() for baseline initial guess rather than a fixed value...
typefit = fit_type_means(roi_sums_xcond_renormalized_typemeans(:,'trace')); 
plot_type_onoff


% 1ni vs 1no
ca_ids = ca_quality_control_by_type(cell_dict_j, cell_info, roi_sums_all_reshape, {'1ni', '1no'}, 1:31, 'qi', 0.5);
figure; plot_grouped_ca(cell_info, ca_ids, roi_sums_means_flatten, {'1ni', '1no'}, [],4,[], 1.1, []);
figure_size_x2([1, 0.7]); set(gca, 'FontSize', 13)

% 2o vs other outer marginal
types_sus_off_wo_2o = {'1ni' '1no'   '1wt' '2an' '2aw' '2i' '3o' '1ws' '3i'   '25' '27' '28'};
ca_ids = ca_quality_control_by_type(cell_dict_j, cell_info, roi_sums_all_reshape, {types_sus_off_wo_2o, '2o'}, 1:31, 'qi', 0.5);
%figure; ax=gca();
%ax.ColorOrder = ax.ColorOrder([1 4 2:7],:); hold on
figure; plot_grouped_ca(cell_info, ca_ids, roi_sums_means_flatten, {types_sus_off_wo_2o, '2o'}, [3 5 2:7],4,[], 1.2, []);
figure_size_x2([1, 0.7]); set(gca, 'FontSize', 13)


% extfig6 -> extfig7
>> tuning_ordered = [tuning_onoff; tuning];
>> tuning_ordered = tuning_ordered(:,order,:);
>> extfig_type_polar(cell_info, ca_dsos, tuning_ordered, [], [], [], 30)


%extfig10: need color
stats_by_type.m

% supplementary
fig_petal

% extfig6
fig_ca_fit


TODOs:
axis tick intervals