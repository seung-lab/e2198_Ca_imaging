%function mapped = map_data(dict_like, query)
function [value, idx] = look_up(lookup_table, key)   % support struct field name? ...SQL

%{
	xxxxxxxxxxxxxxxxxxx order not preserved.....

idx = ismember(lookup_table(:,1), query);
mapped = lookup_table(idx,2:end);
%}

%dict = containers.Map(lookup_table(:,1), lookup_table(:,2));
%value = values(dict, num2cell(key));

	[~, idx] = ismember(key, lookup_table(:,1));

	if isvector(lookup_table)
		value = idx;
	else
		value = lookup_table(idx, 2:end);
	end


end
