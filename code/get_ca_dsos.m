function ca_dsos = get_ca_dsos(tuning, order, cell_dict_j)
% tuning: 2 x 8 x ca_id

ca_dsos = [];
for elem = cell_dict_j.'
    [ca_id, omni_id] = deal(elem(1), elem(2));
    cell_dsos = polar_tuning2(tuning(:,:,ca_id), order, false);
    [cell_dsos.omni_id, cell_dsos.ca_id] = deal(omni_id, ca_id);
    ca_dsos = [ca_dsos; cell_dsos];
end

%{
[ordered, order] = sort(str2num(char(angles)));
ca_dsos_area = get_ca_dsos(roi_sums_area_fakeonoff, order, cell_dict_j);
ca_dsos_area(1:5, :)
cell_info_polarplot_pref_dir(cell_info,ca_dsos_area)

[ordered, order] = sort(str2num(char(angles)));
ca_dsos = get_ca_dsos(tuning_onoff, order, cell_dict_j);
ca_dsos(1:5, :)
cell_info_polarplot_pref_dir(cell_info,ca_dsos)
%}
