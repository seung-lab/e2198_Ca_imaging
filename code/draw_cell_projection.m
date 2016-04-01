
%%{
figure; hold on
for cell_info_elem = cells.'
	omni_id = cell_info_elem.cell_id;
	%vec = [cell_info_elem.soma_coords_warped([2 3])/4; cell_info_elem.soma_coords_warped([2 3])/4 + cell_info_elem.radii_hull * cell_info_elem.asymm_index];
	vec = [cell_info_elem.soma_coords_warped_mip2_zscaled(1:2); cell_info_elem.soma_coords_warped_mip2_zscaled(1:2) + cell_info_elem.radii_hull * cell_info_elem.asymm_index];
	plot(cell_hull{omni_id}(:,1), cell_hull{omni_id}(:,2), '.-', xy_projection{omni_id}(:,1), xy_projection{omni_id}(:,2), '.'); %, ...
		%vec(:,1), vec(:,2), '.-');
	plot(vec(:,1), vec(:,2), 'o-', 'LineWidth', 2);
end
hold off
%}

%{
figure; hold on
for cell_info_elem = cells.'
	omni_id = cell_info_elem.cell_id;
	%vec = [cell_info_elem.soma_coords_warped([2 3])/4; cell_info_elem.soma_coords_warped([2 3])/4 + cell_info_elem.radii_hull * cell_info_elem.asymm_index];
	vec = [cell_info_elem.soma_coords_warped_mip2_zscaled(1:2)];
	plot(cell_hull_rot{omni_id}(:,1), cell_hull_rot{omni_id}(:,2), '.-', xy_projection_rot{omni_id}(:,1), xy_projection_rot{omni_id}(:,2), '.'); %, ...
		%vec(:,1), vec(:,2), '.-');
	plot(vec(:,1), vec(:,2), 'o-', 'LineWidth', 2);
end
hold off
%}