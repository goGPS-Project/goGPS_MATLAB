function [long,lat]=m_xy2ll(X,Y)
% M_XY2LL Converts X,Y to long,lat coordinates using the current projection
%         [LONGITUDE,LATITUDE]=m_ll2xy(X,Y)
%         This is useful for finding locations using ginput.

% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

% 6/Nov/00 - eliminate returned stuff if ';' neglected (thx to D Byrne)
% 4/DEc/11 - isstr to ischar

global MAP_PROJECTION MAP_COORDS

if nargin==0 || ischar(X)
  disp(' Usage:');
  disp(' [LONGITUDE,LATITUDE]=m_xy2ll(X,Y);');
else
  [long,lat]=feval(MAP_PROJECTION.routine,'xy2ll',X,Y);
  if ~strcmp(MAP_COORDS.name,MAP_PROJECTION.coordsystem.name)
    if strcmp(MAP_COORDS.name,'geographic')
      [long,lat]=mc_coords('mag2geo',long,lat);
    else
      [long,lat]=mc_coords('geo2mag',long,lat);      
    end 
  end
end

if nargout==0
  clear long lat
end


