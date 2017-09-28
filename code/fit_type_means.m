% fit_type_means
function [typefit, cellmeanfit] = fit_type_means(roi_sums_xcond_typemeans, type_list, cell_list)

method_list = {'exp', 'alpha', 'exp-exp'}; %, 'exp_w_ti'};

if exist('type_list', 'var')
	debug_fit = 1;
else
	debug_fit = 0;
	type_list = roi_sums_xcond_typemeans.Properties.RowNames;  % gc_types
end

typefit = repmat({table()}, 1, 3);

%for celltype = {'37v'}
%for celltype = {'72n' '1wt'}
for celltype = type_list(:).'

	celltype

	for method_id = 3 %[1 3] %1:length(method_list) % {'exp', 'alpha', 'exp-exp'}

		method = method_list{method_id};

		yactual = roi_sums_xcond_typemeans{celltype,:}(:);
		if isempty(yactual) || sum(isnan(yactual))
			continue;
		end
		[params, cost, ispoor] = fit_single_onoff(yactual, method, debug_fit);

		%typefit{method_id}{celltype, 'params'} = guess;
		%typefit{method_id}{celltype, 'cost'} = cost;
		typefit{method_id}(celltype, {'params' 'cost' 'poor'}) = {params, cost, ispoor};
	end

end


%figure;plot([roi_sums_means_flatten(:,ind) fit16{1}(:,ind) fit16{3}(:,ind)])

