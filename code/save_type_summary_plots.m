%function omni_id = get_omni_id(ca_id)

%{
coeffs16_reshape = reshape(coeffs16{1,2}(:,3:end-16).', 2, 8, 634);
[~, tuning_onoff] = tuning_from_fit(coeffs16{1,2});
tau_offset_list = coeffs16{1,2}(:,1:2);
%}

[~, tuning_onoff] = tuning_from_fit(coeffs16{3,2});
tau_offset_list = coeffs16{3,2}(:,1:3);     % this tau is likely more variable than the single exp tau


figdir = '~/dev/e2198_Ca_imaging/cell_summary/';
mkdir(figdir);

%alltypes = unique(cell_dict_type(:,3)).';
alltypes = cellstr(unique(char(cell_info.type), 'rows')).';
if isempty(alltypes{1})     % remove the empty type
    alltypes = alltypes(2:end);
end

visible = 0;

[ordered, order] = sort(str2num(char(angles)));

%for celltype = {'7i', 'AC'}
for celltype = alltypes
	celltype = celltype{:}	 % convert cell to normal string

	%ca_ids = [cell_dict_type{strcmp(cell_dict_type(:,3), celltype), 2}];
    cells = get_cell_info(cell_info, celltype, 0);
    omni_ids = [cells.cell_id];
    ca_ids = cell_dict_j(ismember(cell_dict_j(:,2), omni_ids), 1);
    if isempty(ca_ids)
        continue;
    end
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
    polar_tuning2(tuning_onoff(:,:,ca_ids), order);

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

