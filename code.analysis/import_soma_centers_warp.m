function soma_voxel_coords = import_soma_centers_warp(filename, outfile)

filename = 'amacrine_soma_centers.csv';

if ~exist('surfaceMapping_10', 'var')
	load /omniData/e2198_reconstruction/m2_surface_mapping_10.mat 		% surfaceMapping_10	1x1 struct
end

soma_import = import_soma_centers(filename);

soma_voxel_coords = [];
offset_warped_vol=[-235 -8 -10]-1;
for row = soma_import.'
	id = row(1);
	soma = row(6:8);
	soma = soma([2 3 1]).'; 
	soma = round(recon_preproc_warp_voxels(soma/4,soma/4,surfaceMapping_10,10,0));  % soma has to be row vec 
	soma=soma([3 1 2])-offset_warped_vol;

	soma_voxel_coords(end+1, :) = [id soma];
end

%outfile = 'sac_soma_centers_m2_warped.csv';
if exist('outfile', 'var') && ~isempty(outfile)
	csvwrite(outfile, soma_voxel_coords)
end

% filename = 'amacrine_soma_centers.csv';
% sac_soma = import_soma_centers_warp()
