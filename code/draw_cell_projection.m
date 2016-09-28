
%%{
figure; hold on

colors=distinguishable_colors(50);
ax = gca();
ax.ColorOrder = colors;
for cell_info_elem = cells.'
	omni_id = cell_info_elem.cell_id;
	if isempty(cell_info_elem.soma_coords_warped_mip2_zscaled) continue; end
	%vec = [cell_info_elem.soma_coords_warped([2 3])/4; cell_info_elem.soma_coords_warped([2 3])/4 + cell_info_elem.radii_hull * cell_info_elem.asymm_index];
	vec = [cell_info_elem.soma_coords_warped_mip2_zscaled(1:2); cell_info_elem.soma_coords_warped_mip2_zscaled(1:2) + cell_info_elem.radii_hull * cell_info_elem.asymm_index];
	%plot(cell_hull{omni_id}(:,1), cell_hull{omni_id}(:,2), '.-', xy_projection{omni_id}(:,1), xy_projection{omni_id}(:,2), '.'); %, ...
		%vec(:,1), vec(:,2), '.-');
		%plot(cell_hull{omni_id}(:,1), cell_hull{omni_id}(:,2), '.-', xy_projection{omni_id}(:,1), xy_projection{omni_id}(:,2), '.'); %, ...
		plot(xy_projection{omni_id}(:,1), xy_projection{omni_id}(:,2))%, '.', 'MarkerSize', 0.1);
	
	%plot(vec(:,1), vec(:,2), 'o-', 'LineWidth', 2);
	%vec = soma_coords_warped(soma_coords(:,1)==omni_id, 2:3);
	%plot(vec(:,1), vec(:,2), 'o-', 'LineWidth', 2);
end
axis equal
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