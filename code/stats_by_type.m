
%{
function cell_stat = cell_info_get_type_property(cell_info, property_type_arg, use_normalized_strat, varargin)

if ~iscell(property_type_arg)
    property_type_arg = {property_type_arg};
end
property_type = property_type_arg{1};
property_arg = property_type_arg(2:end);
%}
%property_type = 'ca_sus_ind_on'; type_list = types_ca_has_on;
%property_type = 'ca_sus_ind_off'; type_list = types_ca_only_off;


%gc_types = list_types(cell_info);

get_super_types
ca_groups = ca_groups2;
n_groups = size(ca_groups,1);
%N = length(n_groups);

colors=distinguishable_colors(n_groups);

roi_sums_xcond_typemeans = table();
roi_sums_xcond_typemeans_rescale = table();
type_stat = table();

figure;
legends = {};
legendentries = [];

all_values = [];
%ca_groups.
%for group = ca_groups.name
for k = 1:n_groups
    %group = strat_groups.types.'
	ca_groups.name{k}
    typelist = ca_groups.types{k};
    if ca_groups.has_on(k)
    	normalization = 'on';
    	property_type = 'ca_sus_ind_on';
    else
    	normalization = 'off';
    	property_type = 'ca_sus_ind_off';
    end

    group_values = [];
	for celltype = typelist(:).'

		roi_sums_xcond_typemeans{celltype, :} = mean(roi_sums_xcondmeans(:, get_ca_ids(cell_dict_j, cell_info, celltype, false)), 2).';
		roi_sums_xcond_typemeans_rescale{celltype, :} = roi_sums_xcond_typemeans{celltype, :} - min(roi_sums_xcond_typemeans{celltype, :});
		roi_sums_xcond_typemeans_rescale{celltype, :} = roi_sums_xcond_typemeans_rescale{celltype, :} / max(roi_sums_xcond_typemeans_rescale{celltype, :});

		switch property_type
		case {'ca_sus_ind_on'}
			[maxi, maxind] = max(roi_sums_xcond_typemeans_rescale{celltype, 7:15});
			roi_sums_xcond_typemeans_rescale{celltype, :} = roi_sums_xcond_typemeans_rescale{celltype, :} / maxi;
			%stat_ca
			tmp = roi_sums_xcond_typemeans_rescale{celltype, :};
			%sus_idx = (1-tmp(13))*sign(tmp(13)-tmp(12));
			%sus_idx = (1-tmp(14))*sign(tmp(14)-tmp(13));
			%sus_idx = (1-tmp(15))*sign(tmp(15)-tmp(14));
			sus_idx = (1-tmp(15)) * (-1);
		

			cells = get_ca_cell_info(cell_dict_j, cell_info, celltype);
		stat = cell_info_get_strat_property(cells, 'sus_on-trans_on');
		stat1 = cell_info_get_strat_property(cells, 'sus_on');
		stat2 = cell_info_get_strat_property(cells, 'trans_on');
		stat = (stat1-stat2) ./ (stat1+stat2);
		%stat = mean(stat1-stat2) ./ mean(stat1+stat2); % only slight visible difference for 27,28

		case {'ca_sus_ind_off'}
			% off
			[maxi, maxind] = max(roi_sums_xcond_typemeans_rescale{celltype, 15:end});
			maxind = maxind+14;
			roi_sums_xcond_typemeans_rescale{celltype, :} = roi_sums_xcond_typemeans_rescale{celltype, :} / maxi;
			tmp = roi_sums_xcond_typemeans_rescale{celltype, :};
			%off
			frame = 24;
			frame = 23;
			sus_idx = (1-tmp(frame))*sign(maxind-frame);

			cells = get_ca_cell_info(cell_dict_j, cell_info, celltype);
			stat1 = cell_info_get_strat_property(cells, 'sus_off');
			stat2 = cell_info_get_strat_property(cells, 'trans_off');
			stat = (stat1-stat2) ./ (stat1+stat2);
			%stat = cell_info_get_strat_property(cells, 'sus');

		end  %switch


		%stat = median(stat);
		stat = mean(stat);

		xx = stat;
		yy = sus_idx;

		yy = yy + 1;
		group_values(end+1,:) = [xx yy];

		%type_stat{celltype, {'sus_idx', 'sus_on'}} = [sus_idx, stat];

		%k = size(type_stat, 1);
		%k = size(roi_sums_xcond_typemeans_rescale, 1);
		kolor = colors(mod(k,n_groups)+1,:);
		%h = scatter(xx, yy, 'filled','MarkerFaceColor',kolor,'MarkerEdgeColor',kolor);
		%h = scatter(xx, yy, 90, ca_groups.style1{k}{:}, 'filled');
		%h = scatter(xx, yy, 90, ca_groups.style2{k}{:}, 'filled');
		h = scatter(xx, yy, 90, ca_groups.style3{k}{:});
		
		text(xx+0.0105,yy, celltype);
		%scatter(repmat(sus_idx,size(stat)),stat, 'filled','MarkerFaceColor',colors(k,:),'MarkerEdgeColor',colors(k,:));
		hold on;
	end

	all_values = [all_values; group_values];
	[R, P] = corrcoef(group_values(:,1), group_values(:,2));
	disp([R(2), P(2)])
	legendentries(end+1) = h;
	legends{end+1} = [ca_groups.name{k} ' (' normalization ')'];
end %for k = 1:n_groups
%legend(legendentries, ca_groups.name);
legend(legendentries, legends);
plot([0 0], [-1 1], '--k');
figure_size_x2();

[R, P] = corrcoef(all_values(:,1), all_values(:,2))

xlabel 'Ca'
ylabel 'strat '
%xlim([-0.8 0])
ylim([-0.8 0]+1)

switch property_type
case {'ca_sus_ind_on'}
	xlabel 'Ca @ 15'
	ylabel(sprintf('strat\n (sus_on-trans_on)/(sus_on+trans_on)'), 'Interpreter', 'none')
	title('on sustainedness')
case {'ca_sus_ind_off'}
	xlabel 'Ca @ 23'
	ylabel(sprintf('strat\n (sus_off-trans_off)/(sus_off+trans_off)'), 'Interpreter', 'none')
	title('off sustainedness')
end %switch

ylabel 'Ca @ 15/23'
xlabel(sprintf('strat\n (sus_on-trans_on)/(sus_on+trans_on)\n or \n (sus_off-trans_off)/(sus_off+trans_off)'), 'Interpreter', 'none')
title('sustainedness')

%asdfsadfsa
%legend(type_stat.Properties.RowNames)
%legend(roi_sums_xcond_typemeans_rescale.Properties.RowNames)
hold off;
%scatter(type_stat{:,1},type_stat{:,2}, 'filled','MarkerFaceColor',colors(k,:),'MarkerEdgeColor',colors(k,:));


n_groups = size(ca_groups,1);
figure;
for k = 1:n_groups
    %group = strat_groups.types.'
    subplot(n_groups, 1, k)
    typelist = ca_groups.types{k};
    if ca_groups.has_on(k)
    	normalization = 'on';
    else
    	normalization = 'off';
    end
    plot_grouped_ca(cell_info, cell_dict_j, roi_sums_means_flatten, typelist, [], 0, normalization, 2.2);
    ylabel(ca_groups.name{k}, 'Rotation', 0, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    l = legend();
    %legend('show', 'Location', 'eastoutside');
    %l.Position(3) = 0.15;
    %l.Location = 'eastoutside';
    ax = gca();
    xtext = ax.XLabel;
    ax.XLabel = [];
    xticklabel = ax.XTickLabel;
    ax.XTickLabel = [];
end
ax.XTickLabel = xticklabel;
ax.XLabel = xtext;

ca_classical_groups = {
	types_oodsgc 'on'
	types_ondsgc 'on'
	types_alpha 'full'	%
	types_alpha_all 'full'	%
};
for k = 1:size(ca_classical_groups,1)
    [typelist, normalization] = ca_classical_groups{k,:};
    figure; plot_grouped_ca(cell_info, cell_dict_j, roi_sums_means_flatten, typelist, [], 0, normalization);
end
