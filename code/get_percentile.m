function [percentiles]=get_percentile(strat,percents, alignment)

% alignment: 'center' (default), 'start' ( strat[1] is 0%, strat[end] is 100%-lastbin ) , 'finish' ( strat[end] is 100% )


if size(strat, 2)==2
	strat = sortrows(strat, 1);   % put into ascending order
	%{
	if strat(1, 1) > strat(end, 1)
		strat = strat(end:-1:1, :);
	end
	%}
	x=strat(:,1);
	s=strat(:,2);
else  % 1
	x = 1:length(strat);
	s = strat;
end

% pad the two sides
s = [0; s; 0];
x = [x(1)-(x(2)-x(1)); x; x(end)+x(end)-x(end-1)];


percents = percents * sum(s);
c=cumsum(s);

if ~exist('alignment', 'var')
	alignment = 'center';
end

for j=1:numel(percents)
	idx=find(c>=percents(j),1,'first');
	switch alignment
	case {'center'}
		xstep = x(idx)-x(idx-1);
		shift = xstep/2;
	case {'finish'}
		xstep = x(idx)-x(idx-1);
		shift = 0;
	case {'start'}
		xstep = x(idx+1)-x(idx);
		shift = xstep;
	end
	if c(idx)~=percents(j)
		%percentiles(j) = x(idx) + (percents(j)-c(idx)) * (x(idx)-x(idx-1)) / (c(idx)-c(idx-1));
		percentiles(j) = x(idx) + shift + (percents(j)-c(idx)) * xstep / s(idx);
	else
		percentiles(j) = x(idx) + shift;
	end
end
