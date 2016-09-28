function [] = polar_tuning(x, order, varargin)
% x: 2x8(xN)
% order: 1st element is expected to be 45deg.

	%polar = @polar2;
	if size(x,3)>1
		x = permute(x, [1 3 2]);
		x = reshape(x, [], 8);
	end
	x = x(:, order);
	n = size(x,1);
	xmean = mean( x .* repmat(exp(i * (1:8)*pi/4), n, 1) , 2 );
	dirs = angle(xmean);
	dir_leng = abs(xmean);
	dirs = [zeros(n,1) dirs];
	dir_leng = [zeros(n,1) dir_leng] * 5;

	tmp = [x(:, end) x];
	range = (0:8)*pi/4;
	xmax = max(x(:));
	if xmax<4
		xmax = 3.9;
	end

	P = polar(0, xmax);
	%P = polar2(0, xmax);
	set(P, 'Visible', 'off')
	hold on;

	styles = {'r', 'k'};
	if 0 && n==1
		styles = {'k'};
	end
	for k = 1:size(tmp,1)
		polar(range, tmp(k, :), styles{mod(k-1,2)+1});
		%polar2(range, tmp(k, :));

		polar(dirs(k,:), dir_leng(k, :), ['.-.' styles{mod(k-1,2)+1}]);
	end
	%polar(range, tmp(2, :), '');
	%title(titletext);

end