function type_info = circstats_by_type(type_info, cell_info, ca_dsos)

addpath('~/seungmount/research/jinseopk/e2198/bin/analysis/CircStat')

gc_types = list_types(cell_info);

if isempty(type_info)
	type_info = table();
end

orig_state = warning;
warning('off','MATLAB:table:RowsAddedExistingVars');

for celltype = gc_types(:).'
%for celltype = {'2an' '2aw' '37c'}

	cells = get_cell_info(cell_info, celltype);
	cell_ids = [cells.cell_id];

	idx = ismember(ca_dsos.omni_id, cell_ids);
	%idx = look_up(ca_dsos.omni_id, cell_ids);

	for name = {'ds_theta', 'os_theta'}	% Oops BUG!! os_theta was -pi/2 to pi/2, needed to x2
		name = name{1};

		%angles = ca_dsos.ds_theta(idx, :);
		angles = ca_dsos{idx, name};
		if strcmp(name, 'os_theta')
			angles = 2*angles;
		end

		p = [];
		for col = 1:size(angles,2)
			%type_info.ds_p = circ_rtest(angles(:,1));
			p(col) = circ_rtest(angles(:,col));
		end
		type_info{celltype, [name(1:2) '_p']} = p;
	end
end

warning(orig_state);

end %func
