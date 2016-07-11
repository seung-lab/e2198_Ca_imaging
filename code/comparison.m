%{
struct('name','37','annotation','ON-OFF DS','cells',[90002 90001 25005 20254 20245 20239 20233 20220 20213 20210 ...
                                            20179 20137 20125 20096 17161 17080 20016 20014 20002 26101 26056 26084 ...
                                            26032 26036 26137 26138 26178 26103 26094 26029 26162 26047 26115 ...
                                            26158 26165]);
%}

ooDSGCdiff = [	% ca id (TBD), em id, x, y, z, ca coords x y(TBD), dist(TBD)
    0   17080	692 18069 2340
    0   20220	1333 8189 7831
    0   20245	1363 6316 3372
    0   26032	996 13273 10788
    0   26036	1384 5398 12129
    0   26056	691 19129 3472
    0   26084	418 19447 8105
    0   26137	1697 1023 11143
    0   26138	1587 970 12055
    %0   26138	1756 936 12065
    0   26165	899 12847 1861
    0   26178	582 17100 1762
];

n = size(ooDSGCdiff,1);
ooDSGCdiff_in_ca = [ooDSGCdiff(:, end-2:end) ones(n,1)] * em_to_ca;

n_ca = size(roi_centers, 1);
for ii = 1:n
	dist = repmat(ooDSGCdiff_in_ca(ii, 1:2), n_ca, 1) - roi_centers;
	dist = sum(dist.^2, 2);
	[val, ind] = min(dist);
	ooDSGCdiff(ii,1) = ind;
	ooDSGCdiff(ii,6) = sqrt(val);
	%[ind sqrt(val)]
	
	[tmp, ind] = sort(dist);
	for ca_ind = ind(1:8).'
		if inpolygon(ooDSGCdiff_in_ca(ii, 1), ooDSGCdiff_in_ca(ii, 2), roi_borders{ca_ind}(:,1), roi_borders{ca_ind}(:,2))
			ooDSGCdiff(ii,1) = ca_ind;
			ooDSGCdiff(ii,6) *= -1;
			break
		end
	end
end

round(ooDSGCdiff(:,[1:2 6]))
find_closest_ca(ooDSGCdiff_in_ca(find(ooDSGCdiff(:,1)==574),:), 3)
%ooDSGCdiff




%{
% assuming global var roi_centers exists.
function ret = find_closest_ca(ca_coord, n_candidates)
	global roi_centers
	n_ca = size(roi_centers, 1);

	dist = repmat(ca_coord(1:2), n_ca, 1) - roi_centers;
	dist = sqrt(sum(dist.^2, 2));
	
	[tmp, ind] = sort(dist);
	ret = [ind(1:n_candidates), tmp(1:n_candidates)];
end
%}



addpath('/usr/people/smu/seungmount/research/jinseopk/e2198/bin/analysis/');
gc=cell_info_typedef_gc();

%onDSGCs = [gc(strcmp({gc.name},'7o')).cells gc(strcmp({gc.name},'7i')).cells]

cell_info=cell_info_set_type();
cell_info=cell_info_get_soma_coord_omni(cell_info)

ca_7i = vertcat(cell_info(strcmp({cell_info.type}, '7i')).soma_coord) * em_to_ca;
get_ca_matches(ca_7i)
vertcat(cell_info(strcmp({cell_info.type}, '7i')).cell_id)

find_closest_ca(ca_7i(find(vertcat(cell_info(strcmp({cell_info.type}, '7i')).cell_id)==26075),:),8)


compute_all_matches.m

backmap = [roi_centers(574, :), 1] * ca_to_em

ctype = '1wt'
[vertcat(cell_info(strcmp({cell_info.type}, ctype)).cell_id) ca_matches_all(strcmp({cell_info.type}, ctype), :)]
%ca_matches_all(strcmp({cell_info.type}, ctype), :)
%vertcat(cell_info(strcmp({cell_info.type}, ctype)).cell_id)

ca_matches_all(strcmp({cell_info.type}, '7i'), :)
vertcat(cell_info(strcmp({cell_info.type}, '7i')).cell_id)


find_closest_ca(ca_all(find(vertcat(cell_info.cell_id)==20147),:),8)


%ca_matches_all(ca_matches_all(:,1)==235,:)

%cell_info(vertcat(cell_info.cell_id)==20147))

cell_info_table = struct2table(cell_info);


%cell_info( find(cellfun(@isempty, {cell_info.soma_coord}) & ~strcmp({cell_info_gc.type}, '')), :);
length(find(~cellfun(@isempty, {cell_info_gc.type})))		%gc cells w type: 310
soma_missing = [cell_info( find(cellfun(@isempty, {cell_info.soma_coord}) & ~cellfun(@isempty, {cell_info_gc.type})), :).cell_id];
ca_matches_all_sort(ismember(ca_matches_all_sort(:,1), soma_missing), :)
existingcells = [cell_info(cellfun(@isempty, {cell_info_gc.type}), :).cell_id];
existingmatches = ca_matches_all_sort(ismember(ca_matches_all_sort(:,1), existingcells), :);

tocheck = ca_matches_all_sort(abs(ca_matches_all_sort(:,2))==abs(ca_matches_all_sort([2:end 1],2)) | abs(ca_matches_all_sort(:,2))==abs(ca_matches_all_sort([end, 1:end-1],2)), :);



%need to modify cell_info_set_type
%cell_info_gc=cell_info_set_type();




%{
addpath('/usr/people/smu/dev/e2198_Ca_imaging/code')
cell_property_mat_file = 'm2_cells_property_160212.mat';
cell_property_mat_file = 'm2_cells_property_warped_nosoma_160212.mat';
load(cell_property_mat_file,'cell_hull','xy_projection');
%}
load('roi_data.mat')
load('roi_data_stimconsts.mat')
compute_all_matches
cell_mapping_verified
cell_info
cell_property_compute_centroids


clear cell_hull xy_projection xy_projection_rot
cell_info=cell_info_get_strat(cell_info,cell_property_mat_file);
types = cellstr(unique(char(cell_info.type), 'rows'));
[cell_ids,prcntile_20,ctype,idx]=cell_info_hist_prcntile(cell_info,{'1ws','1wt','1no','1ni','2an','2aw','2o','2i','3o','3i','4on','4ow','4i','5to','5ti','5so','5si'  '6sn'    '6sw'    '6t'    '7i'    '7o'    '8n'    '8w'    '9n'    '9w'},0.2);
[cell_ids,prcntile_20,ctype,idx]=cell_info_hist_prcntile(cell_info,types(25:43),0.2);

for i=1:length(ctype); fprintf('%d %s %d~%d\n',cell_ids(i),ctype{i},2*(idx(i)-1),2*idx(i)); end


 '63'    '72n'    '72w'    '81'    '82i'    '82o'    '82w'    '83'    '91n'    '91w'
 '5to','5ti','5so','5si' '6sn'    '6sw'    '6t'    '7i'    '7o'    '8n'    '8w'    '9n'    '9w'


division = 38;
div = 17;
div = 21;
div = 23;

div = 22;
C = intersect(ctype(idx > div), ctype(idx <=div))
for c = C.'
	below = idx <= div & strcmp(ctype, c);
	above = idx > div & strcmp(ctype, c);
	fprintf('%s \t %d \t %d \n', c{1}, sum(below), sum(above))
end
	cell_ids(below)
	cell_ids(above)
	%num2str(cell_ids(below), cell_ids(above))
	sum(below), cell_ids(below)
	sum(above), cell_ids(above)


for 
ind = inx < 38; ctype(ind)

 for i=1:length(ctype); fprintf('%d %s %d~%d\n',cell_ids(i),ctype{i},2*(idx(i)-1),2*idx(i)); end



cell_info_multi_plot_strat(cell_info, {'5so'}, 0, 0, 1, '', 0.2)
[gc, ac, bc] = list_types(cell_info);
fignum = 9;
for ctype = gc(:).'
	cell_info_multi_plot_strat(cell_info, ctype, 0, 0, fignum, '', 0.2, Inf, true)
end


addpath('/usr/people/smu/dev/e2198_Ca_imaging/code')
addpath('/usr/people/smu/seungmount/research/jinseopk/e2198/bin/analysis/')
addpath('/usr/people/smu/dev/e2198_Ca_imaging/code.cluster/')
addpath('/usr/people/smu/dev/e2198_Ca_imaging/code.util/')
addpath('/usr/people/smu/dev/e2198_Ca_imaging/code.figure/')
addpath('/usr/people/smu/dev/e2198_Ca_imaging/code.analysis/')
load('cell_info_clustering.mat')
cell_info = cell_info_set_type(cell_info)

cell_info_skel = update_skeleton_strat(cell_info)


[tuning, tuning_onoff] = tuning_from_fit(coeffs16{3,2});
[ordered, order] = sort(str2num(char(angles)));
%ca_dsos = get_ca_dsos(tuning_onoff, order, cell_dict_j);
ca_dsos = get_ca_dsos([tuning_onoff; tuning], order, cell_dict_j);
cell_info_polarplot_pref_dir(cell_info,ca_dsos)


sac_soma_m2_warped = import_soma_centers_warp()


