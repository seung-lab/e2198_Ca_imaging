function draw_e2006_cells(cell_info, query)

cells = get_cell_info(cell_info, query);
n = length(cells);

col = ceil(sqrt(n));
row = ceil(n/col);

figure;
if n==2
	row = 2;
	col = 1;
	figure_size_x2([1 2]);
end
for ii = 1:n
	c = cells(ii);
	a = c.flatskel;

	subplot(row, col*3, ii*3-[2 1])
	%plot(1,1, '')
	%hold on
	plot(a(:,3), a(:,2), '.', 'MarkerSize', 1); grid on; axis equal;
	%plot3(a(:,1), a(:,2), a(:,3), '.'); grid on; axis equal;
	ylim([0 120]);
	ylimm = ylim();
	title([c.annotation '  ' c.type], 'Interpreter', 'none')

	subplot(row, col*3, ii*3)
	plot(a(:,1), a(:,2), '.', 'MarkerSize', 1); grid on;
	ylim(ylimm)
	xlim([-20 120])
end
