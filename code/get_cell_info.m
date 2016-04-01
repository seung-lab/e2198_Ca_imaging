function [cell_info_struct, idx] = get_cell_info(cell_info, query, fuzzy)
	% query: cell id(s), or cell type(s)

	if isempty(query)
		error('empty query');
		return;
	end

	if isnumeric(query) && length(query)==1
		idx = find(vertcat(cell_info.cell_id)==query);
	elseif ischar(query)
		celltype = query;
		idx = find(strcmp({cell_info.type}, celltype));
	elseif isnumeric(query) || ( iscell(query) && ~isempty(query) && isnumeric(query{1}) )
		% multiple cell IDs
		if iscell(query)
			query = cell2mat(query);
		end
		idx = [];
		% return in the queried order
		for cell_id = query(:).'
			idx = [idx find([cell_info.cell_id]==cell_id)];
		end
	elseif iscell(query) && ~isempty(query) && ischar(query{1})
		% multiple cell types
		idx = [];
		% ordered by type
		for ii = 1:length(query)
			celltype = query{ii};
			new = find(strncmp({cell_info.type}, celltype, length(celltype)));
			if isempty(new)
				warning('alsdflsadf')
			end
			idx = [idx new];
		end
	else
		error('invalid argument');
	end
	cell_info_struct = cell_info(idx);
end
	%cell_info_struct = cell_info( vertcat(cell_info.cell_id)==cell_id );