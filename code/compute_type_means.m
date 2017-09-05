function typemeans = compute_type_means(cell_dict_j, cell_info, roi_sums_all_reshape, roi_sums_xcondmeans, renormalize, filtermethod, threshold)
% type temporal response function

if renormalize
	if strcmp(filtermethod, 'none')
		error('garbage in, garbage out')
	end
	% 31*n_rois
	roi_sums_xcondmeans = roi_sums_xcondmeans - repmat(min(roi_sums_xcondmeans,[],1), 31, 1);
	roi_sums_xcondmeans = roi_sums_xcondmeans ./ repmat(max(roi_sums_xcondmeans,[],1), 31, 1);
end

gc_types = list_types(cell_info);

typemeans = table();
for celltype = gc_types(:).'
	if strncmp(celltype{1}, '37', 2) || strncmp(celltype{1}, '7i', 2) || strncmp(celltype{1}, '7o', 2)
		%continue;
	end
	%include = ca_quality_control(roi_sums_all_reshape, framerange, filtermethod, threshold)
	ca_ids = get_ca_ids(cell_dict_j, cell_info, celltype, false);
	include = ca_quality_control(roi_sums_all_reshape(:,:,:,ca_ids), 1:31, filtermethod, threshold); % 'qi', 0.5);
	ca_ids = ca_ids(include);
	if isempty(include)
		typemeans(celltype, :) = {nan(31,1).', 0, []};
	else
		typemeans(celltype, :) = {mean(roi_sums_xcondmeans(:, ca_ids), 2).', length(ca_ids), {ca_ids}};
	end
end
typemeans.Properties.VariableNames = {'trace', 'n', 'ca_ids'};

for celltype = {'weirdos'    'cutoffs'} % {'orphans'    'weirdos'    'cutoffs'}
	typemeans(celltype, :) = [];
end

%{
gc_types = list_types(cell_info);
tmp = ca_quality_control(roi_sums_all_reshape(:,:,:,get_ca_ids(cell_dict_j, cell_info, '72', false)), 1:31, 'qi', 0.7)


cell_dict_j_qi_filtered = ca_quality_control(roi_sums_all_reshape, 1:31, 'qi', 0.7);
cell_dict_j_qi_filtered = cell_dict_j(ismember(cell_dict_j(:,1), cell_dict_j_qi_filtered), :);
%}