function dist = m_xydist(x,y)
% M_XYDIST Spherical earth distance between points in map projection coordinates.  
%   RANGE=M_XYDIST(X,Y) gives the distance in kilometers between
%   successive points in the vectors X and Y which are coordinates
%   in the current map projection. This function is useful for finding
%   distances between points "picked off" a map, e.g.:
%
%   Example:
%   G = GINPUT(2); km_dist = M_XYDIST(G(:,1),G(:,2));
%
%   See also M_LLDIST

% Rich Pawlowicz (rich@ocgy.ubc.ca) 6/Nov/00, based largely on an
% idea by Deirdre Byrne.

% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

global MAP_PROJECTION MAP_VAR_LIST

if isempty(MAP_PROJECTION),
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end;

% get lon, lat coords
[lon,lat] = m_xy2ll(x,y);

dist = m_lldist(lon,lat);



