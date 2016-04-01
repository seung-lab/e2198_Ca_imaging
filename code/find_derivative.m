%# Create the gaussian filter with hsize = [5 5] and sigma = 2
%G = fspecial('gaussian',[5 5],2);
%# Filter it
%Ig = imfilter(I,G,'same');

sigma = 2;
size_ = 10;
x = linspace(-size_ / 2, size_ / 2, size_);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize

alldata = {};

% cell id
ind = 90;
for ind = 90 %c5ti

	peaks = zeros(8,2);
	peaks_filt = zeros(8,2);
	alldir = {};
	for k = 1:8
		y = allmeans(:, k, ind);
		%yfilt = filter (gaussFilter,1, y);
		yfilt = conv (y, gaussFilter, 'same');
		alldir{k}.filtered.x = yfilt;
		%{
		%figure; findpeaks(yfilt);
		der = diff(y);
		[pks,locs] = findpeaks(der);
		alldir{k}.pks = pks;
		alldir{k}.locs = locs;
		[sorted, I] = sort(pks);
		figure; findpeaks(der);
		hold on; plot(1:31, y, 1:31, yfilt); xlim auto; ylim auto;
		peaks(k,:) = sorted(end:-1:end-1);
		if I(end) > I(end-1)
			peaks(k,:) = peaks(k,end:-1:1);
		end
		%}

		%%{
		[a, b] = find_peaks(alldir{k}, y);
		[alldir{k}, peaks(k,:)] = find_peaks(alldir{k}, y);
		[alldir{k}.filtered, peaks_filt(k,:)] = find_peaks(alldir{k}.filtered, yfilt);
		%}
		%figure; findpeaks(der);
		%hold on; plot(1:31, y, 1:31, yfilt); xlim auto; ylim auto;
	end
	alldata{ind}.alldir = alldir;
	alldata{ind}.peaks = peaks;
	alldata{ind}.peaks_filt = peaks_filt;

	%{
	range = (1:8)*45;
	figure; plot(range, peaks(:,1), range, peaks(:,2));
	xlim([45, 360]); ax = gca(); ax.XTick = range;
	%}

	range = (0:8)*pi/4;
	range = range.';
	figure;
	s = max(peaks(:)) / max(peaks_filt(:));
	s = 3;
	polar(range, [peaks(end,1); peaks(:,1)], '');
	hold on;
	polar(range, [peaks(end,2); peaks(:,2)], '-.');
	polar(range, [peaks_filt(end,1); peaks_filt(:,1)].*s, '');
	polar(range, [peaks_filt(end,2); peaks_filt(:,2)].*s, '-.');
	title(num2str(ind));
end