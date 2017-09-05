function names = typename2displayname(names)
% input: string or cell array of strings

waschar = 0;
if ischar(names)
	waschar = 1;
	names = cellstr(names);
end

names(strcmp('91-', names)) = {'91'};

if waschar
	names = char(names);
end
