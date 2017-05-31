function [cell_info, idx] = get_cell_info(cell_info, query, use_partial_match, keep_unknown_cell_ids)
	% query: cell id(s), or cell type(s), or a cell class
	% use_partial_match: for types, a type matches the query as long as it starts with the query string. default = true. 

	if isempty(query)
		error('empty query');
		return;
	end

	switch class(cell_info)
		case 'struct'
			was_table = 0;
		case 'table'
			was_table = 1;
			cell_info = table2struct(cell_info);
		otherwise
			assert(0)
	end
	
	if ~exist('use_partial_match', 'var')
		use_partial_match = true;
	end
	if ~exist('keep_unknown_cell_ids', 'var')
		keep_unknown_cell_ids = false;
	end

	% null cell
	cell_info(end+1).cell_id = 0;  % null cell
	nullcellidx = length(cell_info);
	cell_info(end).class = -1;

	if isnumeric(query) && length(query)==1
		if query < 10	% assume to be class ID
			idx = find(vertcat(cell_info.class)==query);
		else 	% cell ID
			idx = find(vertcat(cell_info.cell_id)==query);
		end
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
	elseif isnumeric(query) || ( iscell(query) && ~isempty(query) && all(cellfun(@isnumeric, query)) )
		% multiple cell IDs
		if iscell(query)
			query = cell2mat(query);
		end
		%tic
		idx = [];
		% return in the queried order
		for cell_id = query(:).'
			new = find([cell_info.cell_id]==cell_id);
			if keep_unknown_cell_ids && isempty(new)
				idx = [idx nullcellidx];
			else
				idx = [idx new];
			end
		end
		%toc
		if length(idx) ~= length(query)
			display('get_cell_info: not all cell IDs are found');
		end
	elseif iscell(query) && ~isempty(query) && iscellstr(query)
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
	elseif iscell(query) && any(cellfun(@isnumeric, query)) && any(cellfun(@ischar, query))
		% types and cell IDs mixed
		idx = [];
		for subquery = query(:).'
			[~, new] = get_cell_info(cell_info, subquery, use_partial_match);
			idx = [idx new];
		end
	else
		error('invalid argument');
	end
	if isempty(idx)
		warning('no result found');
	end
	cell_info = cell_info(idx);

	if was_table
		cell_info = struct2table(cell_info, 'AsArray',true);
	end
end
