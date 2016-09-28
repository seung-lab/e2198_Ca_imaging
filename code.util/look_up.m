%function mapped = map_data(dict_like, query)
function [value, idx] = look_up(lookup_table, key)   % support struct field name? ...SQL

%{
	xxxxxxxxxxxxxxxxxxx order not preserved.....

idx = ismember(lookup_table(:,1), query);
mapped = lookup_table(idx,2:end);
%}

%dict = containers.Map(lookup_table(:,1), lookup_table(:,2));
%value = values(dict, num2cell(key));

	idx = [];
	for src = key(:).'
		src
		idx(end+1) = find(lookup_table(:,1)==src);
	end
	%idx = arrayfun(@(x) find(lookup_vec==x), keys);
	if isvector(lookup_table)
		value = idx;
	else
		value = lookup_table(idx, 2:end);
	end


end
