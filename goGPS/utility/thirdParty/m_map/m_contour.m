function [cs,h]=m_contour(long,lat,data,varargin)
%  M_CONTOUR Draws contour lines on a map
%    M_CONTOUR(LONG,LAT,DATA,...) draw contours on a map. Behavior
%    is the same as for CONTOUR except that LONG and LAT vectors or
%    matrices must be specified.
%
%    [CS,H]=M_CONTOUR(...) returns a contour matrix C and a vector
%    H of handles to LINE or PATCH objects for use by CLABEL.
%
%    See also CONTOUR

% Rich Pawlowicz (rich@ocgy.ubc.ca) 17/Jan/1998
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

% 9/Dec/98 - made sure bad things don't happen if all your lat/long
%            points are out of the plot region.
% 6/Nov/00 - eliminate returned stuff if ';' neglected (thx to D Byrne)

global MAP_PROJECTION

% Have to have initialized a map first

if isempty(MAP_PROJECTION)
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end

if min(size(long))==1 && min(size(lat))==1
 [long,lat]=meshgrid(long,lat);
end

[X,Y]=m_ll2xy(long,lat,'clip','on');

i=isnan(X);      % For these we set the *data* to NaN...
data(i)=NaN;

                 % And then recompute positions without clipping. THis
                 % is necessary otherwise contouring fails (X/Y with NaN
                 % is a no-no. 
if any(i(:)), [X,Y]=m_ll2xy(long,lat,'clip','off'); end 

if any(~i(:))
 [cs,h]=contour(X,Y,data,varargin{:});
 set(h,'tag','m_contour');
else
  cs=[];h=[];
end

if nargout==0
 clear cs h
end

