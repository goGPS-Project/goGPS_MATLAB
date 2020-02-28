function [oceanloadcorr] = ocean_loading_correction(time, XR, XS)

% SYNTAX:
%   [oceanloadcorr] = ocean_loading_correction(time, XR, XS);
%
% INPUT:
%   time = GPS time
%   XR   = receiver position  (X,Y,Z)
%   XS   = satellite position (X,Y,Z)
%
% OUTPUT:
%   oceanloadcorr = ocean loading correction terms (along the satellite-receiver line-of-sight)
%
% DESCRIPTION:
%   Computation of the ocean loading displacement terms.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b7
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:
%  Contributors:     ...
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

%ocean loading displacements matrix, station-dependent (see http://holt.oso.chalmers.se/loading/)
global ol_disp zero_time

oceanloadcorr = zeros(size(XS,1),1);
if (isempty(ol_disp))
    return
end

%terms depending on the longitude of the lunar node (see Kouba and Heroux, 2001)
fj = 1; %(at 1-3 mm precision)
uj = 0; %(at 1-3 mm precision)

%ref: http://202.127.29.4/cddisa/data_base/IERS/Convensions/Convension_2003/SUBROUTINES/ARG.f
tidal_waves = [1.40519E-4, 2.0,-2.0, 0.0, 0.00; ... % M2  - semidiurnal
               1.45444E-4, 0.0, 0.0, 0.0, 0.00; ... % S2  - semidiurnal
               1.37880E-4, 2.0,-3.0, 1.0, 0.00; ... % N2  - semidiurnal
               1.45842E-4, 2.0, 0.0, 0.0, 0.00; ... % K2  - semidiurnal
               0.72921E-4, 1.0, 0.0, 0.0, 0.25; ... % K1  - diurnal
               0.67598E-4, 1.0,-2.0, 0.0,-0.25; ... % O1  - diurnal
               0.72523E-4,-1.0, 0.0, 0.0,-0.25; ... % P1  - diurnal
               0.64959E-4, 1.0,-3.0, 1.0,-0.25; ... % Q1  - diurnal
               0.53234E-5, 0.0, 2.0, 0.0, 0.00; ... % Mf  - long-period
               0.26392E-5, 0.0, 1.0,-1.0, 0.00; ... % Mm  - long-period
               0.03982E-5, 2.0, 0.0, 0.0, 0.00];    % Ssa - long-period

refdate = datenum([1975 1 1 0 0 0]);

[week, sow] = time2weektow(zero_time + time);
dateUTC = datevec(gps2utc(datenum(gps2date(week, sow))));

%separate the fractional part of day in seconds
fday = dateUTC(4)*3600 + dateUTC(5)*60 + dateUTC(6);
dateUTC(4:end) = 0;

%number of days since reference date (1 Jan 1975)
days = (datenum(dateUTC) - refdate);

capt = (27392.500528 + 1.000000035*days)/36525;

%mean longitude of the Sun at the beginning of day
H0 = (279.69668 + (36000.768930485 + 3.03e-4*capt)*capt)*pi/180;

%mean longitude of the Moon at the beginning of day
S0 = (((1.9e-6*capt - 0.001133)*capt + 481267.88314137)*capt + 270.434358)*pi/180;

%mean longitude of the lunar perigee at the beginning of day
P0 = (((-1.2e-5*capt - 0.010325)*capt + 4069.0340329577)*capt + 334.329653)*pi/180;

corr = zeros(3,1);
for k = 1 : 11
    angle = tidal_waves(k,1)*fday + tidal_waves(k,2)*H0 + tidal_waves(k,3)*S0 + tidal_waves(k,4)*P0 + tidal_waves(k,5)*2*pi;
    corr  = corr + fj*ol_disp(1).matrix(1:3,k).*cos(angle + uj - ol_disp(1).matrix(4:6,k)*pi/180);
end
corrENU(1,1) = -corr(2,1); %east
corrENU(2,1) = -corr(3,1); %north
corrENU(3,1) =  corr(1,1); %up

%displacement along the receiver-satellite line-of-sight
XRcorr = local2globalPos(corrENU, XR);
corrXYZ = XRcorr - XR;
for s = 1 : size(XS,1)
    LOS  = XR - XS(s,:)';
    LOSu = LOS / norm(LOS);
    % oceanloadcorr(s,1) = dot(corrXYZ,LOSu);
    oceanloadcorr(s,1) = sum(conj(corrXYZ).*LOSu);
end
