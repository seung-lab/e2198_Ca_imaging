global roi_centers
global roi_borders


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
           '034_p2_g7_DSbars_200um_export.mat'];

rawdata_fn = ['034_p2_g7_DSbars_200um.cfd';
           '035_p2_g7_DSbars_200um.cfd';
           '036_p2_g7_DSbars_200um.cfd';
           '037_p2_g7_DSbars_200um.cfd';
           '038_p2_g7_DSbars_200um.cfd';
           '039_p2_g7_DSbars_200um.cfd';
           '040_p2_g7_DSbars_200um.cfd';
           '041_p2_g7_DSbars_200um.cfd';
           '042_p2_g7_DSbars_200um.cfd'];
       
typeclrs = [[0 125 125];...
           [0 255 255];...
           [0 0 255];...
           [255 0 255];...
           [255 0 100];...
           [255 255 0];...
           [255 150 0];...
           [0 255 0]];


roi_fn = 'AVG_4classes_043_p2_g7_ZoomedOut_center_011.rpb';
roi_struct = load_rois_from_rpb(roi_fn,typeclrs./255);
for count=1:length(roi_struct.roi_ids)
    roi_borders{count} = [get(roi_struct.border_h(count),'XData')' get(roi_struct.border_h(count),'YData')'];
    roi_centers(count,:) = mean(roi_borders{count}, 1);
end


% most of cells_ds
em_coords = [ % Ca id, em id, x, y, z
23 26103 1331 1942 1691
36 26115 1610 1076 4434
44 26101 1330 2760 2375
83 25005 1613 4682 5778
91 20254 1666 5043 7807
159 20137 1460 5773 10539
178 26130 1673 1573 10290
188 26094 1609 3774 8303
208 26158 1080 9406 1397
209 26162 949 9586 909
265 90002 1139 6836 2414
278 20179 1026 12051 4561
298 20233 1154 10213 6007
305 20210 1304 8896 5775
307 20213 1221 8174 4620
316 17161 1381 6979 6140
331 20239 1209 9850 7240
389 26029 1171 10579 11388
434 20016 752 14698 2675
481 20014 972 15769 3333
508 90001 787 15294 4401
515 26047 520 16745 5155
553 20125 953 14255 7474
590 20002 763 16027 10179
602 20096 817 14847 9118

%185 10015 1751 1072 11605
];


indices = em_coords(:, 1); 
ca_coords = [roi_centers(indices, :), ones(length(indices),1)];
%ca_coords = [roi_centers(indices, :), zeros(length(indices),1), ones(length(indices),1)];
%ca_coords = [roi_centers(indices, :), 1000*ones(length(indices),1), ones(length(indices),1)];
em_coords_only = [em_coords(:, 3:5), ones(length(indices),1)];

ca_to_em = ca_coords \ em_coords_only;

em_to_ca = em_coords_only \ ca_coords;

%{
% 505 is likely focused on x=~600 plane
index = 505
[roi_centers(index, :), 1] * ca_to_em
%}


addpath('/usr/people/smu/seungmount/research/jinseopk/e2198/bin/analysis/');
cell_info=cell_info_set_type();


cell_info=cell_info_get_soma_coord_omni(cell_info);
em_all = vali(cell_info.soma_coord);
ca_all = [vertcat(em_all{:}) ones(1576,1)] * em_to_ca;
ca_matches_all = get_ca_matches(ca_all);
n = size(ca_all,1);
ca_matches_all(:,2) = round(ca_matches_all(:,2)*10);
ca_matches_all = [ca_matches_all, zeros(n,4)];
for ii = 1:n
  backmap = [roi_centers(ca_matches_all(ii, 1), :), 1] * ca_to_em;
  backmap = backmap(1:3);
  em_dist = norm(em_all{ii} - backmap);
  ca_matches_all(ii, 3) = em_dist;
  ca_matches_all(ii, 4:6) = round(backmap);
end
ca_matches_all = round(ca_matches_all);

