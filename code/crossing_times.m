function [times,signs]=f(t,x,x_star)

% this function returns an array of times at which x crosses x_star.
% whether a crossing has occured is determined by looking at whether
% x-x_star changes sign.  in this context, zero is considered a positive
% number.
%
% t and x must be vectors of same size
% x_star can be scalar or vector of same size at t and x
% each element of signs is positive or negative depending on the
% sign of the crossing
%
% ALT 2/2001

[n,t_min,t_max,T,dt,fs]=time_info(t);
x_sub_star=x-x_star;
x_sub_star_sign=sign(x_sub_star);
x_sub_star_sign=(x_sub_star_sign>=0);
crossings=x_sub_star_sign(2:n)~=x_sub_star_sign(1:n-1);
crossing_indices=find(crossings);
n_crossings=length(crossing_indices);
times=zeros(n_crossings,1);
for i=1:n_crossings
  j=crossing_indices(i);
  times(i)=t(j)+dt*x_sub_star(j)/(x_sub_star(j)-x_sub_star(j+1));
  signs(i)=sign(x_sub_star(j+1)-x_sub_star(j));
end
