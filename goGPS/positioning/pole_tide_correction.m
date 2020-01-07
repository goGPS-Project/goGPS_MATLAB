function [poletidecorr] = pole_tide_correction(time, XR, XS, SP3, phiC, lam)

% SYNTAX:
%   [poletidecorr] = pole_tide_correction(time, XR, XS, SP3, phiC, lam);
%
% INPUT:
%   time = GPS time
%   XR   = receiver position  (X,Y,Z)
%   XS   = satellite position (X,Y,Z)
%   SP3  = structure containing precise ephemeris and ERP data
%   phiC = receiver geocentric latitude (rad)
%   lam  = receiver longitude (rad)
%
% OUTPUT:
%   poletidecorr = pole tide correction terms (along the satellite-receiver line-of-sight)
%
% DESCRIPTION:
%   Computation of the pole tide displacement terms.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:
%  Contributors:     Andrea Gatti...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------
if (nargin < 5)
    [~, lam, ~, phiC] = cart2geod(XR(1,1), XR(2,1), XR(3,1));
end

poletidecorr = zeros(size(XS,1),1);

%interpolate the pole displacements
if (~isempty(SP3.ERP))
    if (length(SP3.ERP.t) > 1)
        m1 = interp1(SP3.ERP.t, SP3.ERP.m1, time, 'linear', 'extrap');
        m2 = interp1(SP3.ERP.t, SP3.ERP.m2, time, 'linear', 'extrap');
    else
        m1 = SP3.ERP.m1;
        m2 = SP3.ERP.m2;
    end
    
    deltaR   = -33*sin(2*phiC)*(m1*cos(lam) + m2*sin(lam))*1e-3;
    deltaLam =  9* cos(  phiC)*(m1*sin(lam) - m2*cos(lam))*1e-3;
    deltaPhi = -9* cos(2*phiC)*(m1*cos(lam) + m2*sin(lam))*1e-3;

    corrENU(1,1) = deltaLam; %east
    corrENU(2,1) = deltaPhi; %north
    corrENU(3,1) = deltaR;   %up
    
    %displacement along the receiver-satellite line-of-sight
    XRcorr = local2globalPos(corrENU, XR);
    corrXYZ = XRcorr - XR;
    for s = 1 : size(XS,1)
        LOS  = XR - XS(s,:)';
        LOSu = LOS / norm(LOS);
        poletidecorr(s,1) = -dot(corrXYZ,LOSu);
    end
end


