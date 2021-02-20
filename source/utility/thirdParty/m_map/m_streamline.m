function h=m_streamline(long,lat,u,v,varargin)
% M_STREAMLINE Makes a quiverplot on a map (QUIVER-style)
%    M_STREAMLINE(LONG,LAT,U,V) draws well-spaced streamlines (with direction 
%    arrows) with components (U,V) at the points (LONG,LAT) on the currently 
%    defined map.The arrays LONG and LAT, which define the coordinates for U and 
%    V, must be monotonic, but do not need to be uniformly spaced. The 
%    matrices LONG,LAT,U,V must all be the same size. U and V contain the 
%    eastward and northward components of velocity (in m/s or equivalent, NOT
%    degrees lat/long  per sec or equivalent). Arrow scaling is automatic.
% 
%   Note - this is basically a call to STREAMSLICE, which has some limitations and
%   cannot be used for ALL projections. In particular, STREAMSLICE requires that
%
%            [X,Y]=m_ll2xy(LONG,LAT,'clip','point')
%
%   returns X/Y matrices that themselves must be monotonic and plaid (as if produced
%   by MESHGRID). This means that a) none of the points are from outside the map
%   boundaries (since the M_LL2XY call will turn them into NaN), and b) you are
%   probably limited to cylindrical projections (miller, mercator, equidistant 
%   cylindrical). 
%
%   M_STREAMLINE(...,density) modifies the automatic spacing of the streamlines. 
%   Density must be greater than 0. The default value is 1; higher values 
%   produce  more streamlines on each plane. For example, 2 produces 
%   approximately twice as many streamlines, while 0.5 produces approximately 
%   half as many.
% 
%   M_STREAMLINE(...,'arrowsmode') determines if direction arrows are present 
%   or not. arrowmode can be:
% 
%   arrows   - Draw direction arrows on the streamlines (default).
%   noarrows - Do not draw direction arrows.
% 
%   M_STREAMLINE(...,'method') specifies the interpolation method to use. 
%   method can be
% 
%   linear   - Linear interpolation (default)
%   cubic    - Cubic interpolation
%   nearest  - Nearest-neighbor interpolation
%   See interp3 for more information on interpolation methods.
% 
%   M_STREAMLINE(axes_handle,...) plots into the axes object with the handle 
%   axes_handle instead of into the current axes object (gca).
% 
%   h = M_STREAMLINE(...) returns a vector of handles to the line objects 
%   created.
%
%   [vertices arrowvertices] = M_STREAMLINE(...) returns two cell arrays of 
%   vertices  for drawing the streamlines and the arrows. You can pass these 
%   values to any of the streamline drawing functions (streamline, streamribbon, 
%   streamtube).
% 

% Shi Weiheng (tfoterye@gmail.com) 29/Jun/17
%
% Based on m_quiver written by Prof Rich Pawlowicz.
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.
%
%  Nov/2017 - cleaned up comments.

global MAP_PROJECTION MAP_VAR_LIST

% Have to have initialized a map first

if isempty(MAP_PROJECTION)
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end



[X,Y]=m_ll2xy(long,lat,'clip','point');

if any(isnan(X(:)))
  error(['M_Map : ' mfilename ' : InvalidInputs - input data includes points outside the map area']);
end
  

[XN ,YN ]=m_ll2xy([long(:) long(:)]',[lat(:) lat(:)+.001]','clip','off');
[XE ,YE ]=m_ll2xy([long(:) long(:)+(.001)./cos(lat(:)*pi/180)]',[lat(:) lat(:)]','clip','off');
mU=u.*reshape(diff(XE),size(lat))*1000 + v.*reshape(diff(XN),size(lat))*1000;
mV=u.*reshape(diff(YE),size(lat))*1000 + v.*reshape(diff(YN),size(lat))*1000;

h=streamslice(X,Y,mU,mV,varargin{:});
set(h,'tag','m_streamline');

if nargout==0
 clear h
end
