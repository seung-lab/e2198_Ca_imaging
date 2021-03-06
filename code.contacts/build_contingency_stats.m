%function [contact_table, stats, cell_stats] = all_contacts2(cell_info)
function [type_stats, cell_stats, cell_info_new] = build_contingency_stats(cell_info)

s = load_contact_vars();
contact_vars = s;

cells = get_cell_info_table(cell_info, s.cell_ids, [], false);
%cells = cells(:,{'cell_id', 'type'});
mask = cellfun(@isempty,cells.type); % cantknow=0, GC=1, AC=2, BC=3, uncertain=4, etc=5
cells.type(mask) = cellstr(num2str(9000 + 100 * cells.class(mask)));	% *100 to avoid confusion with normal GC types with partial matching in get_cell_info


% enum = unique(cells.type) %contains '  0' or ' 100'
cell_info_new = table2struct(cells);
cells.type = categorical(cells.type);  % arh! leading spaces are removed

enum = unique(cells.type);  % arh! contains '0' converted from '  0'
%type_stats = table('VariableNames',{'cell1_type' 'cell2_type' })

type_stats = table();
cell_stats = table();

partial = table(cells.type(1), cells.type(1), 'VariableNames',{'cell1_type' 'cell2_type' });
for type1 = enum(:).'
	display(type1)
	type1str = char(type1);
	cells1 = get_cell_info(cell_info_new, type1str);
	contactvoxels = load_raw_contacts([cells1.cell_id]);
tic
	for type2 = enum(:).'
		type2str = char(type2);
		[aggregate, individual] = contact_summary(cell_info_new, type1str, type2str, 0, contact_vars, contactvoxels);
		if isempty(aggregate)
			continue;  % one cell to self
		end

		partial{:,:} = [type1 type2];

		type_stats(end+1,:) = [partial aggregate];
		cell_stats = [cell_stats; ...
				individual repmat(partial,size(individual,1),1) ];
	end
	toc
end

type_statsReorder = sortrows(type_stats, {'cell2_type', 'cell1_type'});
if ~isequal(type_statsReorder.cell2_type, type_stats.cell1_type) || ~isequal(type_statsReorder.cell1_type, type_stats.cell2_type)
	% this trick below is relying on it being symmetric
	warning('ahhhhhhhhhhh')
end

type_stats.fasciculation2way = sqrt(type_stats.fasciculation .* type_statsReorder.fasciculation);
type_stats.surfInGridsWContacts_surfInMaskedGrids = sqrt(type_stats.surf1inGrids1wContacts ./ type_stats.maskedSurface1 ...
	.* type_statsReorder.surf1inGrids1wContacts ./ type_statsReorder.maskedSurface1);

