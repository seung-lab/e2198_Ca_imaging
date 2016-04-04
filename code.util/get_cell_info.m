function [cell_info_struct, idx] = get_cell_info(cell_info, query, use_partial_match)
	% query: cell id(s), or cell type(s)
	% use_partial_match: for types, a type matches the query as long as it starts with the query string. default = true. 

	if isempty(query)
		error('empty query');
		return;
	end
	
	if ~exist('use_partial_match', 'var')
		use_partial_match = true;
	end

	if isnumeric(query) && length(query)==1
		idx = find(vertcat(cell_info.cell_id)==query);
	elseif ischar(query)
		celltype = query;
		if use_partial_match
			idx = find(strncmp({cell_info.type}, celltype, length(celltype)));
		else
			idx = find(strcmp({cell_info.type}, celltype));
		end
		if isempty(idx)
			warning(sprintf('No cell found for type "%s"', celltype));
		end
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
			if use_partial_match
				new = find(strncmp({cell_info.type}, celltype, length(celltype)));
			else
				new = find(strcmp({cell_info.type}, celltype));
			end
			if isempty(new)
				warning(sprintf('No cell found for type "%s"', celltype));
			end
			idx = [idx new];
		end
	else
		error('invalid argument');
	end
	cell_info_struct = cell_info(idx);
end
	%cell_info_struct = cell_info( vertcat(cell_info.cell_id)==cell_id );