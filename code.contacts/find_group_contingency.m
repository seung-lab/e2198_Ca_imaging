function [total, individual, rowsum, group1, group2] = find_group_contingency(contingency_mat, contingency_labels, group1, group2)

[found1, idx1] = ismember(group1, contingency_labels);
[found2, idx2] = ismember(group2, contingency_labels);
if any(~found1) && nargout>1 && nargout<4
	error('not all group1 members are found')
end
if any(~found2) && nargout>1 && nargout<5
	error('not all group2 members are found')
end
group1 = group1(found1);
group2 = group2(found2);
idx1 = idx1(found1);
idx2 = idx2(found2);
individual = contingency_mat(idx1, idx2, :);

total = sum(reshape(individual, [], size(individual,3)));	% keeping 3rd dimension intact
total = reshape(total, 1, 1, []);
rowsum = sum(individual, 2);

function idx = ismember2(varargin)
	[~, idx] = ismember(varargin);


