function [dsos, ds_r, ds_theta, os_r, os_theta, r_mean] = polar_tuning2(x, order, varargin)
% x: 2*8(*N) or n*8(*N)
% return: 2x1 (or 2xN) or n*N (i.e. removed the 2nd dimension from x)

nvarargin = length(varargin);
optargs = {true};
optargs(1:nvarargin) = varargin;
[make_polar_plot] = optargs{:};


	%polar = @polar2;
	[n_onoff, ~, N] = size(x);
	if size(x,3)>1
		x = permute(x, [1 3 2]);
		x = reshape(x, [], 8);
	end
	x = x(:, order);
	n = size(x,1);	% = n_onoff*N
	xlenmean = mean( x , 2 );
	xmean = mean( x .* repmat(exp(i * (1:8)*pi/4), n, 1) , 2 );
	orientationmean = mean( x .* repmat(exp(i * (1:8)*pi/2), n, 1) , 2 );
	dirs = angle(xmean);
	dir_leng = abs(xmean);
	dir_leng_normalized = dir_leng ./ xlenmean;
	orientations = angle(orientationmean) / 2;
	orientation_leng = abs(orientationmean);
	orientation_leng_normalized = orientation_leng ./ xlenmean;

	[ds_r, ds_theta, os_r, os_theta, r_mean] = deal(dir_leng, dirs, orientation_leng, orientations, xlenmean);
	tmp = [ds_r, ds_theta, os_r, os_theta, r_mean];
	tmp = reshape(tmp, n_onoff, N, 5);
	[ds_r, ds_theta, os_r, os_theta, r_mean] = deal(tmp(:,:,1), tmp(:,:,2), tmp(:,:,3), tmp(:,:,4), tmp(:,:,5));
	tmp = permute(tmp, [2 1 3]);  % N*n_onoff*5. Rows can be displayed more easily by the "table" data type
	% WTH: cell2table({[1 2] [3 4]}) has different element/variable/column types from cell2table({[1; 2] [3; 4]}). 
	dsos = permute(num2cell(tmp, [2]), [1 3 2]); % convert to Nx5 cell array
	%dsos = cell2struct(dsos, {'ds_r', 'ds_theta', 'os_r', 'os_theta', 'r_mean'}, 2);
	dsos = cell2table(dsos, 'VariableNames', {'ds_r', 'ds_theta', 'os_r', 'os_theta', 'r_mean'});

	
	%%--- plotting 
	if ~make_polar_plot
		return;
	end

	dirs = [zeros(n,1) dirs];
	dir_leng = [zeros(n,1) dir_leng] * 5;
	dir_leng_normalized = [zeros(n,1) dir_leng_normalized];
	orientations = [orientations-pi orientations];
	orientation_leng = [orientation_leng orientation_leng] * 5;
	orientation_leng_normalized = [orientation_leng_normalized orientation_leng_normalized];

	tmp = [x(:, end) x];
	range = (0:8)*pi/4;
	xmax = x(~any(x>1000, 2), :);	% filter out rows having rare bad fits
	if isempty(xmax)  % no rows left...
		xmax = x;
	end
	xmax = max(xmax(:));
	if xmax<4
		xmax = 3.9;
	end

	if make_polar_plot == 2  %
		for k = 1:size(tmp,1)
			polarplot(range, tmp(k, :), 'LineWidth', 1);
			hold on
			polarplot(range, tmp(k, :), 'LineWidth', 1);
			ax = gca;
			ax.ColorOrderIndex = 2;
			polarplot(dirs(k,:), dir_leng_normalized(k, :).*xlenmean(k), '.-', 'LineWidth', 1);

		    ax.ThetaTick=0:45:315;
		    ax.GridAlpha = 1;
		    ax.GridColor = 0.8*[1 1 1];
		    ax.ThetaTickLabel(2:2:8)={''};
		    %ax.ThetaTickLabel(3:8)={''};
		    ax.GridLineStyle = '--';

            ax.RTick(1) = [];	% remove 0 tick
		end
		return;
	end

	%plots = [3, 2];

	%{
	for k = 1:6
		subplot(3, 2, k)
		if k<=4
			P = polar(0, xmax);
		else
			P = polar(0, 0.5);
		end
		set(P, 'Visible', 'off')
		hold on;
	end

	for k = 1:size(tmp,1)
		onoff = mod(k-1,2)+1;
		subplot(3, 2, onoff)
		polar(range, tmp(k, :));

		subplot(3, 2, onoff+2)
		polar(dirs(k,:), dir_leng(k, :), '.-.');

		subplot(3, 2, onoff+4)
		polar(dirs(k,:), dir_leng_normalized(k, :), '.-.');
	end
	%polar(range, tmp(2, :), '');
	%title(titletext);
	%}

polar = @polarplot;
	%%{
	for k = 1:n_onoff*3
		subplot(n_onoff, 3, k)
		if mod(k-1, 3)==0
			P = polar(0, xmax);
		else
			P = polar(0, 0.5);
		end
		set(P, 'Visible', 'off')
		hold on;
	end

	for k = 1:size(tmp,1)
		onoff = mod(k-1,n_onoff);
		subplot(n_onoff, 3, onoff*3 + 1)
		polar(range, tmp(k, :), 'LineWidth', 1);
		rlim([0 xmax])

		subplot(n_onoff, 3, onoff*3 + 2)
		if onoff==0
			title('orientation')
		end
		polar(orientations(k,:), orientation_leng_normalized(k, :), '.-.', 'LineWidth', 1);
		rlim([0 0.5])

		subplot(n_onoff, 3, onoff*3 + 3)
		if onoff==0
			title('direction')
		end
		polar(dirs(k,:), dir_leng_normalized(k, :), '.-.', 'LineWidth', 1);
		rlim([0 0.5])
	end
	%polar(range, tmp(2, :), '');
	%title(titletext);
	%}

end