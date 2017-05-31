function warped_voxel_coords = get_m2_warped_omni_coords(points)
% points: Nx3
% return: Nx3

tic
display('loading warping map')
load /omniData/e2198_reconstruction/m2_surface_mapping_10.mat 		% surfaceMapping_10	1x1 struct
toc
display('done')

warped_voxel_coords = [];
offset_warped_vol=[-235 -8 -10]-1;
for point = points.'
	if any(isnan(point))
		point = nan(1,3);
	else
		point = point([2 3 1]).'; 
		% WARNING: not safe with NaN's (hangs), or possibly other invalid (far away) coords
		point = round(recon_preproc_warp_voxels(point/4,point/4,surfaceMapping_10,10,0));  % point has to be row vec 
		point=point([3 1 2])-offset_warped_vol;
	end

	warped_voxel_coords(end+1, :) = point;
end
