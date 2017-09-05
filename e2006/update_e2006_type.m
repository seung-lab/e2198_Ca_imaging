function cell_info = update_e2006_type(cell_info, assigntype, annotations)
% CAUTION: assuming cell_info.annotations has all values.

%strncmp({e2006.annotation}, 'gao', 3)
idx = ismember({cell_info.annotation}, annotations);
if ischar(assigntype)
	assigntype = {assigntype};
end
[cell_info(idx).type] = deal(assigntype{:});

