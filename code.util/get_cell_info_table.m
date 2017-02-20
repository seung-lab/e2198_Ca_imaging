function cell_info_table = get_cell_info_table(cell_info_struct, varargin)

	if length(varargin)>0
		cell_info_struct = get_cell_info(cell_info_struct, varargin{:});
	end

	cell_info_table = struct2table(cell_info_struct, 'AsArray',true);

	% convert [] to '' for the 'type' field
	%cell_info_table.type = cellstr(char(cell_info_struct.type));  % assuming no trailing spaces...
	cell_info_table.type(cellfun(@isempty,cell_info_table.type)) = {''};
