function [longGEO,latGEO,phiVecGEO,thetaVecGEO]=m_mag2geo(longMAG,latMAG,phiVecMAG,thetaVecMAG)
% M_MAG2GEO  Converts magnetic to geographic coordinates.
%   [longGEO,latGEO]=M_MAG2GEO(longMAG,latMAG) converts geomagnetic 
%   (dipole) coordinates to geographic coordinates.  IGRF 2000 is used 
%   to determine the location of the geomagnetic dipole. All in units of 
%   degrees with + longitudes east. All variables can be scalar or matrix 
%   but must have the same size.
%
%    Vector rotations can be carried using
% 
%  [longGEO,latGEO,phiVecGEO,thetaVecGEO]=M_MAG2GEO(longMAG,latMAG,phiVecMAG,thetaVecMAG)
%
%   where
%
% phiVecMAG   - east component of the vector in geomagnetic coordinates
% thetaVecMAG - north component of the vector in geomagnetic coordinates
% phiVecGEO   - east component of the vector in geographic coordinates
% thetaVecGEO - north component of the vector in geographic coordinates
%
% See also M_GEO2MAG

% References:
%
% Hapgood, M.A., Space Physics Coordinate Transformations:
% A User Guide, Planet. Space Sci., Vol. 40, N0. 5, 1992.

% R. Pawlowicz (rich@ocgy.ubc.ca)
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

if nargin==2
   [longGEO,latGEO]=mc_coords('mag2geo',longMAG,latMAG);
elseif nargin==4
   [longGEO,latGEO,phiVecGEO,thetaVecGEO]=mc_coords('mag2geo',longMAG,latMAG,phiVecMAG,thetaVecMAG);
else
   error('Wrong number of input parameters');
end  

