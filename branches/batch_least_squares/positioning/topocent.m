function [Az, El, D] = topocent(Xr, Xs)

% SYNTAX:
%   [Az, El, D] = topocent(Xr, Xs);
%
% INPUT:
%   Xr = receiver coordinates (X,Y,Z)
%   Xs = satellite coordinates (X,Y,Z)
%
% OUTPUT:
%   D = rover-satellite distance
%   Az = satellite azimuth
%   El = satellite elevation
%
% DESCRIPTION:
%   Computation of satellite distance, azimuth and elevation with respect to
%   the receiver.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) Kai Borre
% Kai Borre 09-26-97
%
% Adapted by Mirko Reguzzoni, Eugenio Realini, 2009
%----------------------------------------------------------------------------------------------

%conversion from geocentric cartesian to geodetic coordinates
[phi, lam, h] = cart2geod(Xr(1), Xr(2), Xr(3)); %#ok<NASGU>

%new origin of the reference system
X0(:,1) = Xr(1) * ones(size(Xs,1),1);
X0(:,2) = Xr(2) * ones(size(Xs,1),1);
X0(:,3) = Xr(3) * ones(size(Xs,1),1);

%computation of topocentric coordinates
cl = cos(lam); sl = sin(lam);
cb = cos(phi); sb = sin(phi);
F = [-sl -sb*cl cb*cl;
      cl -sb*sl cb*sl;
       0    cb   sb];
local_vector = F' * (Xs-X0)';
E = local_vector(1,:)';
N = local_vector(2,:)';
U = local_vector(3,:)';
hor_dis = sqrt(E.^2 + N.^2);

if hor_dis < 1.e-20
   %azimuth computation
   Az = 0;
   %elevation computation
   El = 90;
else
   %azimuth computation
   Az = atan2(E,N)/pi*180;
   %elevation computation
   El = atan2(U,hor_dis)/pi*180;
end

i = find(Az < 0);
Az(i) = Az(i)+360;

%receiver-satellite distance
D = sqrt(sum((Xs-X0).^2 ,2));
