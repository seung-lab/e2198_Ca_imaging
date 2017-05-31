%https://www.mathworks.com/matlabcentral/fileexchange/49179-nanpolyarea-x-y-/content/nanpolyarea.m
%nanpolyarea(x,y)
%by Bas Altena
 
%27 Jan 2015 (Updated 27 Jan 2015)
%Calculate the area of polygons, separated by NaN's

function area = nanpolyarea(x,y)
%NANPOLYAREA Area of polygons delimited by NaN's.
%   NANPOLYAREA(X,Y) returns the area of the polygon specified by
%   the vertices in the vectors X and Y.  If X and Y are matrices
%   of the same size, then POLYAREA returns the area of
%   polygons defined by the columns X and Y.  If X and Y are
%   arrays, POLYAREA returns the area of the polygons in the
%   first non-singleton dimension of X and Y.  
%
%   The polygon edges must not intersect.  If they do, POLYAREA
%   returns the absolute value of the difference between the clockwise
%   encircled areas and the counterclockwise encircled areas.
%
%   NANPOLYAREA(X,Y,DIM) returns the area of the polygons specified
%   by the vertices in the dimension DIM.
%
%   Class support for inputs X,Y:
%      float: double, single

%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 1.12.4.3 $  $Date: 2010/08/23 23:11:57 $
%

%% initialization
if nargin==1 
  error(message('MATLAB:polyarea:NotEnoughInputs')); 
end

if ~isequal(size(x),size(y)) 
  error(message('MATLAB:polyarea:XYSizeMismatch')); 
end

if (any(x<0) || any(y<0))
    error('method can not be used with negative coordinates');
end
    

[x,nshifts] = shiftdim(x);
y = shiftdim(y);
siz = size(x);

%% estimation
% based on "Math for Map makers" by Allan, pp.119

if ~isempty(x),
    
    % find all delimiters
    f = find(isnan(x));
    xNew = x;
    yNew = y;
    if sum(f) == 0
        % one polygon
        area = abs(sum(y.*[x(2:end); x(1)]) - sum(x.*[y(2:end); y(1)]))/2;
    else
        % multiple polygons
        if size(f,1) == 1
            xNew(f) = x(1); % place first element to NaN place
            yNew(f) = y(1);
        else
            xNew(f) = [x(1); x(f(1:end-1)+1)];
            yNew(f) = [y(1); y(f(1:end-1)+1)];
        end
        xNew = [xNew; x(f((end))+1)];
        yNew = [yNew; y(f(end)+1)];
        
        xNew(f+1) = []; % remove first element of polygon
        xNew(1) = []; 
        
        yNew(f+1) = []; % remove first element of polygon
        yNew(1) = []; 
        
        area = abs(sum(y(~isnan(y)).*xNew) - sum(x(~isnan(x)).*yNew))/2;
    end
else
  area = sum(x); % SUM produces the right value for all empty cases
end

area = shiftdim(area,-nshifts);