function [out, onoff]= find_peaks(out, x)

	%figure; findpeaks(yfilt);
	der = diff(x);
	[pks,locs] = findpeaks(der);
	out.diff = der;
	out.pks = pks;
	out.locs = locs;
	[sorted, I] = sort(pks);
	%figure; findpeaks(der);
	%hold on; plot(1:31, y, 1:31, yfilt); xlim auto; ylim auto;
	onoff = sorted(end:-1:end-1);
	if I(end) > I(end-1)
		onoff = onoff(end:-1:1);
	end
end