function g = generate_t_ti(tau, varargin) %, n = 8, t1t2 = [], method='exp')

%TODO: func name generate_alphas no longer match reality (not alpha)

nvarargin = length(varargin);
optargs = {8, [], 'exp'};
optargs(1:nvarargin) = varargin;
[n, t1t2, method] = optargs{:};

if isempty(t1t2)
t1 = 8.9468; t2 = 16.5457;
	if strcmp(method,'exp-exp')
		t1 = 9.3993; t2 = 16.7412;
	end
elseif length(t1t2)==2
	t1 = t1t2(1);
	t2 = t1t2(2);
end

% g(t) = t * exp(-t/tau)
% y = c + A1 * g(t - t1) + A2 * g(t - t2) + A2 * g(t - t2 + 31)
%   = c + A1 * g1 + A2 * g2

%%{
% t*exp(t/tau) and exp(t/tau)
t = 1:31*n;
t = 1:31*n*2;
if length(t1t2)==16
	ti = t1t2;
else
t1s = (0:n-1)*31 + t1;
t2s = (0:n-1)*31 + t2;
ti = [t1s; t2s];
ti = ti(:);
end

% for the initial downslope. Note these are always >0.
t_t1_ = t + 31*n - ti(end-1);
t_t2_ = t + 31*n - ti(end);

g = zeros(n*2, length(t));
switch method
case {'exp', 'exp_w_ti'}
	for k = 1:n*2
		t_ti = t - ti(k);  %t_ti(t_ti<0) = 0;
		%g(k,:) = t_ti .* exp(-t_ti/tau);
		g(k,:) = exp(-t_ti/tau);
		%g(k,t_ti<=0) = 0;
		g(k,t_ti<0) = 0;
	end
	% assuming last two peaks in the last condition causes the initial downslope at the beginning
	g(end-1,:) = g(end-1,:) + exp(-t_t1_/tau);
	g(end,:) = g(end,:) + exp(-t_t2_/tau);
case 'alpha'
	for k = 1:n*2
		t_ti = t - ti(k);  t_ti(t_ti<0) = 0;
		g(k,:) = t_ti .* exp(-t_ti/tau);
	end
	% assuming last two peaks in the last condition causes the initial downslope at the beginning
	g(end-1,:) = g(end-1,:) + t_t1_ .* exp(-t_t1_/tau);
	g(end,:) = g(end,:) + t_t2_ .* exp(-t_t2_/tau);
case 'exp-exp'
	for k = 1:n*2
		t_ti = t - ti(k);  t_ti(t_ti<0) = 0;
		g(k,:) = exp(-t_ti/tau(1)) - exp(-t_ti/tau(2));
		g(k,t_ti<=0) = 0;
	end
	% assuming last two peaks in the last condition causes the initial downslope at the beginning
	if 0  % original method
	g(end-1,:) = g(end-1,:) + exp(-t_t1_/tau(1)) - exp(-t_t1_/tau(2));
	g(end,:) = g(end,:) + exp(-t_t2_/tau(1)) - exp(-t_t2_/tau(2));
		g = g(:,1:end/2);
	else  % alternative method using all peaks (not just last two)
		g = g(:,1:end/2) + g(:,1+end/2:end);
	end
otherwise
	error('unknown method')
end

%{
% assuming last peak causes the initial downslope at the beginning
t_t2_ = t - t2 + 31;
%g(end,:) = g(end,:) + t_t2_ .* exp(-t_t2_/tau);
g(end,:) = g(end,:) + exp(-t_t2_/tau);
%}

return
%}

%{
% t*exp(t/tau)
t = 1:31*n;
t = 1:31*(n+1);
t_t1 = t - t1;  t_t1(t_t1<0) = 0;
t_t2 = t - t2;  t_t2(t_t2<0) = 0;
t_t2_ = t - t2 + 31;

g1 = t_t1 .* exp(-t_t1/tau);
g2 = t_t2 .* exp(-t_t2/tau);
%g2_ = t_t2_ .* exp(-t_t2_/tau);
g2_ = g2(31+1:end);

g = zeros(n*2, 31*n);
for k = 1:n
	g(2*k-1,31*(k-1)+1:end) = g1(1:31*(n-k+1));
	g(2*k  ,31*(k-1)+1:end) = g2(1:31*(n-k+1));
end
% assuming last peak causes the initial downslope at the beginning
g(2*k, :) = g(2*k, :) + g2_;
%}

%%{
% exp - exp
t = 1:31*(n+1);
t_t1 = t - t1;  t_t1(t_t1<0) = 0;
t_t2 = t - t2;  t_t2(t_t2<0) = 0;
%t_t2_ = t - t2 + 31;

g1 = exp(-t_t1/tau(1)) - exp(-t_t1/tau(2)); g1(t_t1<=0) = 0;
g2 = exp(-t_t2/tau(1)) - exp(-t_t2/tau(2)); g2(t_t2<=0) = 0;
g1_ = g1(31+1:end);
g2_ = g2(31+1:end);

g = zeros(n*2, 31*n);
for k = 1:n
	g(2*k-1,31*(k-1)+1:end) = g1(1:31*(n-k+1));
	g(2*k  ,31*(k-1)+1:end) = g2(1:31*(n-k+1));
end
% assuming last 2 peaks causes the initial downslope at the beginning
g(2*k, :) = g(2*k, :) + g2_;
g(2*k-1, :) = g(2*k-1, :) + g1_;
%}


end
