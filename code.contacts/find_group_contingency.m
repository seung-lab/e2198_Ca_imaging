function [total, individual, rowsum] = find_group_contingency(contingency_mat, contingency_labels, group1, group2)

[~, idx1] = ismember(group1, contingency_labels);
[~, idx2] = ismember(group2, contingency_labels);
individual = contingency_mat(idx1, idx2);

total = sum(individual(:));
rowsum = sum(individual, 2);

function idx = ismember2(varargin)
	[~, idx] = ismember(varargin);


