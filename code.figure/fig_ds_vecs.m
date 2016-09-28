function fig_ds_vecs(dsos, cell_info)

	if exist('cell_info', 'var')
		cell_info = get_cell_info(cell_info, 1);  % class=1: GC only
	    idx = find(ismember(dsos.omni_id, [cell_info.cell_id]) & dsos.ds_r(:,3)>0);
	else
		%idx = 1:size(dsos,1);
		idx = find(dsos.ds_r(:,3)>0);
	end
    n = length(idx)
    zz = zeros(1, n);
size(dsos.ds_r(:,3))
	figure;
	nplots = 3;
	figure_size_x2([nplots,1])
	subplot(1, nplots, 1)
	theta = dsos.ds_theta(idx,3);
	theta = pi/2+theta;  % to final "standard" coord
	r = dsos.ds_r(idx,3);
	polarplot([zz; theta.'], [zz; r.']);
	rlim([0 10])
	subplot(1, nplots, 2)
	h = rose(theta);
	subplot(1, nplots, 3)
	h = rose(theta, 60);
	%h.XData

return
	hold on
	figure;fill(h.XData(57:64), h.YData(57:64), h.Color)
	hold on
	h.XData(61:64), h.YData(61:64)
	plot(h.XData(61), h.YData(61), 'o')
	%figure;plot(h.XData(1:60), h.YData(1:60))
	%figure;fill(h.XData, h.YData, h.Color)
end
