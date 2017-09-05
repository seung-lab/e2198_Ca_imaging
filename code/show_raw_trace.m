%c=496;figure;plot(roi_sums_raw(:,c)); set(gca, 'XTick', stimframes(:, patchdict(c))); grid on
%c=496;figure;plot(roi_sums_raw(:,c)); tk=stimframes(:, patchdict(c)); set(gca, 'XTick', sort([tk+1/.128;tk+2/.128])); grid on

function investigate_raw_trace(roi_sums_raw, stimframes, patchdict, ca_cell_id)

c=ca_cell_id;figure;plot(roi_sums_raw(:,c)); tk=stimframes(:, patchdict(c)); set(gca, 'XTick', sort([tk(1:8:end);tk+1/.128;tk+2/.128])); grid on
figure_size_x2([3,1]);
title(ca_cell_id)

