function filtered = table_sql_where_in(table_records, colname, matchvals, not_in)
% not_in: optional. Specify 'not' for inverse selection or empty for normal selection.

	if exist('not_in', 'var') && ~isempty(not_in)
		if ischar(not_in) && strcmpi(not_in, 'not')
			not_in = 1;
		else
			warning('invalid input for not_in')
			not_in = 0;
		end
	else
		not_in = 0;
	end
    
    filtercol = table_records{:, colname};
    matchvals = matchvals(:).';
    idx = repmat(filtercol, size(matchvals)) == repmat(matchvals, size(filtercol)); % ...forgot there's an ismember() function
    idx = any(idx, 2);

    if not_in
    	idx = ~idx;
    end
    filtered = table_records(idx,:);
end
