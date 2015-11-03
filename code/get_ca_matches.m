% assuming global var roi_centers and roi_borders exists.
% ca_coords: n x 3 
% returns: [ca_id dist] as one matrix. Dist is to center of ROI. "Negative" dist means inside ROI region.
function ret = get_ca_matches(ca_coords)
	global roi_borders
	n = size(ca_coords,1);
	ret = zeros(n, 2);
	for ii = 1:n
		close_rois = find_closest_ca(ca_coords(ii, 1:2), 8); % 8 closest cells
		ret(ii,1:2) = close_rois(1, 1:2);	% closest
		
		% for ca_ind = close_rois(:,1).'		% for the 8 closest cells
		for k = 1:8		% for the 8 closest cells
			ca_ind = close_rois(k,1);
			if inpolygon(ca_coords(ii, 1), ca_coords(ii, 2), roi_borders{ca_ind}(:,1), roi_borders{ca_ind}(:,2))
				if ca_ind ~= ret(ii,1)
					fprintf('closest was %d %d , inside is %d %d \n', ret(ii,:), ca_ind, close_rois(k, 2));
				end
				ret(ii,1) = ca_ind;
				ret(ii,2) = -close_rois(k, 2);
				break
			end
		end
	end
end
