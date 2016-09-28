function [cell_info_struct, idx] = get_ca_cell_info(cell_dict_j, cell_info, query, varargin)
    cells = get_cell_info(cell_info, query, varargin{:});
    idx = ismember([cells.cell_id], cell_dict_j(:,2));
    cell_info_struct = cells(idx);
end
