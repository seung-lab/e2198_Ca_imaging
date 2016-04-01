
coeffs16_reshape = reshape(coeffs16{1,2}(:,3:end-16).', 2, 8, 634);
%coeffs16_reshape = reshape(coeffs16{1}(:,3:end).', 2, 8, 634);
%coeffs_ordered = coeffs_reshape(:, order, :);
%vert_offsets = coeffs16{1}(:,1).';

figdir = '~/dev/e2198_Ca_imaging/cell_summary/';
figdir = '/omniData/e2198_reconstruction/e2198_Ca_imaging/new_cell_summary/'
mkdir(figdir);
%{
addpath('/usr/people/smu/seungmount/research/jinseopk/e2198/bin/analysis/')
gc = cell_info_typedef_gc();
cell_info=cell_info_set_type();
%}

%for ind = 1:634
for ind = 41 %48
	%close all;

	omni_id = cell_dict(cell_dict(:,2)==ind, 1);
	if ~isempty(omni_id) %find(cell_dict(:,2)==ind)
		celltype = cell_info([cell_info.cell_id]==omni_id).type;
		subdir = celltype;
		folder = [figdir subdir];
		if ~exist(folder, 'dir')
			mkdir(folder);
		end
	else
		folder = figdir;
	end

	figure;
	subplot(2,1,1)
	polar_tuning(coeffs16_reshape(:,:,ind), order);


	subplot(2,1,2)
	%plot([roi_sums_means_flatten(:,ind) fit16{1}(:,ind) fit16{3}(:,ind)])
	plot([roi_sums_means_flatten(:,ind) fit16{1,2}(:,ind)]);hold on;plot([fit16{1,1}(:,ind) fit16{3,2}(:,ind)], '-.')
	ax = gca(); ax.XTick = t1t2_to_ti(9,16); grid on;
	if diff(ylim())<8
		lim = ylim;
		if lim(2)<6
			lim(2) = 6;
		end
		ylim(lim);
	end



	title(num2str(ind));

	%{
	print(gcf, '-r300', sprintf('%s\\%s.png',folder,num2str(ind)), '-dpng');
	close;
	%}
end
