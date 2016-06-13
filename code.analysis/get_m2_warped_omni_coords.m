function warped_voxel_coords = get_m2_warped_omni_coords(points)
% points: Nx3
% return: Nx3

load /omniData/e2198_reconstruction/m2_surface_mapping_10.mat 		% surfaceMapping_10	1x1 struct

warped_voxel_coords = [];
offset_warped_vol=[-235 -8 -10]-1;
for point = points.'
	point = point([2 3 1]).'; 
	point = round(recon_preproc_warp_voxels(point/4,point/4,surfaceMapping_10,10,0));  % point has to be row vec 
	point=point([3 1 2])-offset_warped_vol;

	warped_voxel_coords(end+1, :) = point;
end
