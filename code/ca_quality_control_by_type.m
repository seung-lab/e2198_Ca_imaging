function [ca_ids, cell_ids, filtered_cell_info] = ca_quality_control_by_type( ...
	cell_dict_j, cell_info, roi_sums_all_reshape, types, framerange, filtermethod, threshold, suppress_recursion)

self = @ca_quality_control_by_type;

if isempty(framerange)
	framerange = 1:31;
end
if ~exist('suppress_recursion', 'var')
	suppress_recursion = 0;
end

if iscell(types) && ~suppress_recursion
	% just calling self multiple types

	agg = cell(0,3);
	for t = types(:).'
		[agg{end+1, :}] = self(cell_dict_j, cell_info, roi_sums_all_reshape, t{1}, framerange, filtermethod, threshold, 1);
		%[agg{end+1, :}] = deal(a, b, c);
	end
	[ca_ids, cell_ids, filtered_cell_info] = deal(agg(:,1), agg(:,2), agg(:,3));
else

	[ca_ids, cell_info] = get_ca_ids(cell_dict_j, cell_info, types, false);
	include = ca_quality_control(roi_sums_all_reshape(:,:,:,ca_ids), framerange, filtermethod, threshold);
	ca_ids = ca_ids(include);
	filtered_cell_info = cell_info(include);
	cell_ids = [filtered_cell_info.cell_id];

end
