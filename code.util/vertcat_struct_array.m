function s = vertcat_struct_array(s1, s2)
% abiding to what would happen for s1(3).non_exist_field = []


f1 = fieldnames(s1);
f2 = fieldnames(s2);

f = setdiff(f1, f2);
for n = f(:).'
	s2(1).(n{1}) = [];
end

f = setdiff(f2, f1);
for n = f(:).'
	s1(1).(n{1}) = [];
end

s = [s1; s2];
