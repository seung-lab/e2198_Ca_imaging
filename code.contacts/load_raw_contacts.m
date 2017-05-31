function [contacts, cellids] = load_raw_contacts(cell_ids)

basepath = 'contacts/raw3d_445-1167';

contacts = {};

cellids = [];
for cell_id = cell_ids(:).'
	%{
	search_pattern=sprintf('%s/%d.mat',basepath,c.cell_id);
	search_result=dir(search_pattern);
	if isempty(search_result)
	    error('cell not found');
	end
	file_path=sprintf('%s/%s',dir_sac2gc,search_result.name);
	%}

	file_path = sprintf('%s/%d.mat', basepath, cell_id);
	if ~exist(file_path, 'file')
		display(sprintf('cell not found: %d \n', cell_id));
		continue
	end
	tmp = load(file_path);  % var: 'contacts'
	contacts{cell_id} = tmp.contacts;
	cellids(end+1) = cell_id;
end
%display(sprintf('cells:  \n  %s', num2str(cellids)))
