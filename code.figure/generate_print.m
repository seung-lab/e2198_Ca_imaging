%generate_print.m

%cell_info_polarplot_pref_dir2(cell_info,ca_dsos, [0 1 2])


cell_info_polarplot_pref_dir2(cell_info,ca_dsos, [1 2 0])

cell_info_polarplot_pref_dir(cell_info,ca_dsos, [], 0, 0, 4)

% fig 6
fig4_ca
typefit = fit_type_means(roi_sums_xcond_typemeans);
plot_type_onoff


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
