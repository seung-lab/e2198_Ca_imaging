%function omni_id = get_omni_id(ca_id)


coeffs16_reshape = reshape(coeffs16{1,2}(:,3:end-16).', 2, 8, 634);
tau_offset_list = coeffs16{1,2}(:,1:2);


figdir = '~/dev/e2198_Ca_imaging/cell_summary/';
mkdir(figdir);

alltypes = unique(cell_dict_type(:,3)).';

visible = 0;

[ordered, order] = sort(str2num(char(angles)));

ca_dsos = [];
for elem = cell_dict_j.'
    [ca_id, omni_id] = deal(elem(1), elem(2));
    cell_dsos = polar_tuning2(coeffs16_reshape(:,:,ca_id), order, false);
    [cell_dsos.omni_id, cell_dsos.ca_id] = deal(omni_id, ca_id);
    ca_dsos = [ca_dsos; cell_dsos];
end
%ca_dsos = struct2table()
%return;

%for celltype = {'7i', 'AC'}
for celltype = alltypes
	celltype = celltype{:}	 % convert cell to normal string

	%ca_ids = [cell_dict_type{strcmp(cell_dict_type(:,3), celltype), 2}];
    cells = get_cell_info(cell_info, celltype);
    if isempty(cells)
        continue;  % the cell_dict_type array had 'AC', 'uncertain' classes as types
    end
    omni_ids = [cells.cell_id];
    ca_ids = cell_dict_j(ismember(cell_dict_j(:,2), omni_ids), 1);
	ca_ids
    ncells = length(ca_ids);
    taus = tau_offset_list(ca_ids,1);
    tau_mean = mean(taus);
    tau_std = std(taus);


    %summary_fig_h = figure();
    summary_fig_h = figure('Position',[0 0 1200 800]);
    if ~visible
    	summary_fig_h.Visible = 'off';
	end
    polar_tuning2(coeffs16_reshape(:,:,ca_ids), order);

	%subplot(2, 3, 2); title('');
	%subplot(2, 3, 3); title('');
	subplot(2, 3, 1)
    title(sprintf('%s :  n = %d, tau = %.1f (%.1f)', celltype, ncells, tau_mean, tau_std));

    folder = [figdir '/' celltype];
	if ~exist(folder, 'dir')
		mkdir(folder);
	end
	filepath = sprintf('%s\\%s.png',folder,['type_summary_' celltype]);
    print(summary_fig_h, '-r300', filepath, '-dpng');
    % also make a copy in the type's own folder
    %copyfile(filepath, );
    if ~visible
    	close(summary_fig_h);
    end
end

