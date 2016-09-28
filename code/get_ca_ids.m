function ca_ids = get_ca_ids(cell_dict_j, cell_info, query, varargin)
	% warning: not order preserved
    cells = get_cell_info(cell_info, query, varargin{:});
    idx = ismember(cell_dict_j(:,2), [cells.cell_id]);
    ca_ids = cell_dict_j(idx,1);
end
