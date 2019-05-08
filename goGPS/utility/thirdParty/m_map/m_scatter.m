function h = m_scatter(long,lat,varargin)
% M_SCATTER Scatter/bubble plot
%    M_SCATTER(LON,LAT,Y,S,C) displays colored circles at the locations specified
%    by the vectors LON and LAT (which must be the same size).  
% 
%    S determines the area of each marker (in points^2). S can be a
%    vector the same length a X and Y or a scalar. If S is a scalar, 
%    MATLAB draws all the markers the same size. If S is empty, the
%    default size is used.
%    
%    C determines the colors of the markers. When C is a vector the
%    same length as X and Y, the values in C are linearly mapped
%    to the colors in the current colormap. When C is a 
%    length(X)-by-3 matrix, it directly specifies the colors of the  
%    markers as RGB values. C can also be a color string. See ColorSpec.
% 
%    M_SCATTER(LON,LAT) draws the markers in the default size and color.
%    M_SCATTER(LON,LAT,S) draws the markers at the specified sizes (S)
%    with a single color. This type of graph is also known as
%    a bubble plot.
%    M_SCATTER(...,M) uses the marker M instead of 'o'.
%    M_SCATTER(...,'filled') fills the markers.
%  
%    H = M_SCATTER(...) returns handles to the scatter objects created.
% 
%    Use M_PLOT for single color, single marker size scatter plots.
%

% RP 2016

global MAP_PROJECTION MAP_VAR_LIST

if isempty(MAP_PROJECTION)
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end

if nargin < 2
  help m_scatter
  return
end

[x,y] = m_ll2xy(long,lat);

h=scatter(x,y,varargin{:});

set(h,'tag','m_scatter');

if nargout == 0
  clear h
end

return
