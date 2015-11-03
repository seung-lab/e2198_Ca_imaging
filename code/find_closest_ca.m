% assuming global var roi_centers exists.
function ret = find_closest_ca(ca_coord, n_candidates)
	global roi_centers
	n_ca = size(roi_centers, 1);

	dist = repmat(ca_coord(1:2), n_ca, 1) - roi_centers;
	dist = sqrt(sum(dist.^2, 2));
	
	[tmp, ind] = sort(dist);
	ret = [ind(1:n_candidates), tmp(1:n_candidates)];
end
