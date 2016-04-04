
%cell_property_mat_file = 'm2_cells_property_160212.mat';
cell_property_mat_file = 'm2_cells_property_nosoma_160212.mat';
%clear cell_hull xy_projection xy_projection_rot

if ~exist('cell_hull', 'var')
	load(cell_property_mat_file,'cell_hull','xy_projection');
	% cell arrays with indices being cell id
end

%centroids_hull = [];
%radii_hull = [];
centroids_hull = cell(length(cell_hull), 1);
radii_hull = cell(length(cell_hull), 1);
radii_hull2 = cell(length(cell_hull), 1);
area_hull = cell(length(cell_hull), 1);
area_projection = cell(length(cell_hull), 1);
cell_hull_simplified = cell(length(cell_hull), 1);

xy_projection_rot = cell(length(cell_hull), 1);
cell_hull_rot = cell(length(cell_hull), 1);

% https://en.wikipedia.org/wiki/Centroid#Centroid_of_polygon
%for hull = cell_hull
for ii = 1:length(cell_hull)
	hull = cell_hull{ii};
	if isempty(hull)
		continue;
	end
	% x_i * y_i+1 - x_i+1 * y_i
	xyxy = hull(:,1) .* hull([2:end 1], 2) - hull(:,2) .* hull([2:end 1], 1);

	% x_i + x_i+1
	xx = hull(:,1) + hull([2:end 1], 1);
	yy = hull(:,2) + hull([2:end 1], 2);

	area = sum(xyxy) / 2; %polyarea
	
	centroids_hull{ii} = [sum(xx.*xyxy) sum(yy.*xyxy)] / 6 / area;
	[K, A] = convhull(hull, 'simplify', true);
	hull = hull(K, :);
	cell_hull_simplified{ii} = hull;
	assert(isequal(area,A))
	%radii_hull{ii} = sqrt(area / pi);
	area_hull{ii} = area;
	area_projection{ii} = size(xy_projection{ii},1);
end


if ~exist('surfaceMapping_10', 'var')
	load m2_surface_mapping_10.mat 		% surfaceMapping_10	1x1 struct
end


cells = cell_info;
%cells = struct2cell(cells);
tmp = ~cellfun(@isempty, {cells.soma_coord});
cells = cells(tmp);
%soma_coords = [vertcat(cells.cell_id) vertcat(cells.soma_coord)];

if ~exist('soma_coords_warped', 'var')
	soma_coords = [];
	for cell_info_elem = cells.'
	soma = cell_info_elem.soma_coord([2 3 1]);
	%soma_coords(end+1, :) = 4 * recon_preproc_warp_voxels(soma/4,soma/4,surfaceMapping_10,10,1);
	%soma_coords(end+1, :) = 4 * recon_preproc_warp_voxels(soma/4,soma/4,surfaceMapping_10,10,0);
	soma = recon_preproc_warp_voxels(soma/4,soma/4,surfaceMapping_10,10,0);
	soma(2) = soma(2) * 23/16.5;
	soma_coords(end+1, :) = soma;
	end

	%soma_coords = vertcat(cells.soma_coord);
	%soma_coords = 4 * recon_preproc_warp_voxels(soma_coords/4,soma_coords/4,surfaceMapping_10,10,1);
	%soma_coords(:,1) = soma_coords(:,1) * 23/16.5;
	soma_coords = [vertcat(cells.cell_id) soma_coords];
	%soma_coords(:,4) = soma_coords(:,4) * 23/16.5;
end
soma_coords_warped = soma_coords;


for ii = 1:size(soma_coords, 1)
	cell_id = soma_coords(ii, 1);
	idx = find(vertcat(cell_info.cell_id)==cell_id);

	hull = cell_hull_simplified{cell_id};
	hmin = min(hull);
	hmax = max(hull);
	mincut = hull < repmat(hmin+10, size(hull, 1), 1);	% magic number 10
	maxcut = hull > repmat(hmax-10, size(hull, 1), 1);
	%diff(hull(mincut(:,1), :))
	if any([mincut; maxcut] > 1)
		%is_cutoff = true;
		is_cutoff = false;
	else
		is_cutoff = false;
	end
	exclude = [hmin; hmax];
	excluded = ismember(hull(:,1), exclude(:,1)) | ismember(hull(:,2), exclude(:,2));
	if ~isempty(find(excluded))
		%fprintf('%d excluded vertex %d \n', cell_id, find(excluded))
		%fprintf('%d excluded vertex %d \n', cell_id, sum(excluded))
	end
	count = sum(~excluded);
	vecs = hull(~excluded, :) - repmat(soma_coords(ii, 2:3), count, 1);
	radii_hull{cell_id} = mean( sqrt(sum(vecs.^2, 2)) );

	tmp = permute(hull, [1 3 2]);	% n*1*2
	tmp = repmat(tmp, 1, size(tmp,1));	% n*n*2
	pairwiseDiff = tmp - permute(tmp, [2 1 3]);
	pairwiseDiff = sum(pairwiseDiff.^2, 3);	% squared distance, n*n*1
	max_diameter = sqrt(max(pairwiseDiff(:)));
	
	cell_info(idx).is_cutoff = is_cutoff;
	cell_info(idx).asymm_index = (centroids_hull{cell_id} - soma_coords(ii, 2:3)) / radii_hull{cell_id};
	cell_info(idx).asymm_index = (centroids_hull{cell_id} - soma_coords(ii, 2:3)) / sqrt(area_hull{cell_id} / pi);

	cell_info(idx).radii_hull = radii_hull{cell_id};
	cell_info(idx).soma_coords_warped_mip2_zscaled = soma_coords(ii, 2:4);
	%cell_info(idx).soma_coords_warped = soma_coords(ii, 2:4);
	%todo, dir of warped 

	cell_info(idx).max_diameter = max_diameter;
	cell_info(idx).area_hull = area_hull{cell_id};
	cell_info(idx).area_projection = area_projection{cell_id};
end

% use direction of the 2an group
cells = get_cell_info(cell_info, {'2an'});
tmp =  vertcat(cells.asymm_index);
tmpnorm = sqrt(sum(tmp.^2, 2));
tmp = tmp./repmat(tmpnorm, 1, 2);
dir2an = mean(tmp);
%% ans =   -0.4009    0.6884
trans = [dir2an; -dir2an(2) dir2an(1)];
trans = trans.';

for k = 1:size(soma_coords, 1)
	cell_id = soma_coords(k, 1);
	idx = find(vertcat(cell_info.cell_id)==cell_id);

	hull = cell_hull_simplified{cell_id};

	hull = hull - repmat(soma_coords(k, 2:3), size(hull,1), 1);
	hull = hull * trans;
	[xi,yi,ii] = polyxpoly(hull(:,1), hull(:,2), [0 0], [-1e10, 1e10]);
	if size(ii,1)==0
		fprintf('no intersects: %d \n', cell_id);
		% Left / Right
		if hull(1,1) < 0
			asym = Inf;
		else
			asym = 0;
		end
	else
		assert(size(ii,1) == 2); % two intersect points
		ii = ii(:,1);
		hull1 = [hull(1:ii(1), :); [xi, yi];  hull(ii(2)+1:end, :)];
		hull2 = [xi([2 1]), yi([2 1]); hull(ii(1)+1:ii(2), :); xi(2), yi(2)];
		a1 = polyarea(hull1(:,1), hull1(:,2));
		a2 = polyarea(hull2(:,1), hull2(:,2));
		% Left / Right
		if hull1(1,1) < 0 || hull2(3,1) > 0
			asym = a1 / a2;
		else
			asym = a2 / a1;
		end
	end
 
 	proj = xy_projection{cell_id};
 	proj = proj - repmat(soma_coords(k, 2:3), size(proj,1), 1);
	proj = proj * trans;
	asym_p = sum(proj(:,1)<0) / sum(proj(:,1)>0);

	cell_info(idx).asymm_2an = asym;
	cell_info(idx).asymm_2an_prj = asym_p;
	% asym / asym_p 
	if asym / asym_p > 5 || asym / asym_p < 0.2
		warning(sprintf('%d  %g %g', cell_id, asym, asym_p))
	end

	xy_projection_rot{cell_id} = proj;
	cell_hull_rot{cell_id} = hull;
end



cells = get_cell_info(cell_info, '37');
%cells = get_cell_info(cell_info, 90001);
cells = get_cell_info(cell_info, 90002);
cells = get_cell_info(cell_info, 70244);
cells = get_cell_info(cell_info, 70241);
%cells = get_cell_info(cell_info, 70237);

run draw_cell_projection.m
%{
figure; hold on
for cell_info_elem = cells(1:6).'
	omni_id = cell_info_elem.cell_id;
	vec = [cell_info_elem.soma_coords_warped([2 3])/4; cell_info_elem.soma_coords_warped([2 3])/4 + cell_info_elem.radii_hull * cell_info_elem.asymm_index];
	plot(cell_hull{omni_id}(:,1), cell_hull{omni_id}(:,2), '.-', xy_projection{omni_id}(:,1), xy_projection{omni_id}(:,2), '.'); %, ...
		%vec(:,1), vec(:,2), '.-');
	plot(vec(:,1), vec(:,2), 'o-', 'LineWidth', 2);
end
hold off
%}




%{
omni_id = 70008;

figure;
plot(cell_hull{omni_id}(:,1), cell_hull{omni_id}(:,2), '.-')
plot(xy_projection{omni_id}(:,1), xy_projection{omni_id}(:,2), '.')


omni_id = 20164;
figure;plot(cell_hull{omni_id}(:,1), cell_hull{omni_id}(:,2), '.-', xy_projection{omni_id}(:,1), xy_projection{omni_id}(:,2), '.')
%}