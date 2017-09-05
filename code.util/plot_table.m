function plot_table(tbl, field)  % dir / virt/horz
%e.g. plot_ds_os(ca_dsos_bytype, 'ds_p')

n = size(tbl,1);

figure;
plot(tbl.(field), 1:n);
ax = gca();

ax.YTick = 1:n;
ax.YTickLabels = tbl.Properties.RowNames(:);
ylim([0 n+1])
