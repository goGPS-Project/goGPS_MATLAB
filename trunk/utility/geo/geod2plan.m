function [EAST, NORTH, utm_zone] = geod2plan(lat, lon)

% SYNTAX:
%   [EAST, NORTH, utm_zone] = geod2plan(lat, lon);
%
% INPUT:
%   lat = latitude [rad]
%   lon = longitude [rad]
%
% OUTPUT:
%   EAST = EAST coordinate [m]
%   NORTH = NORTH coordinate [m]
%   utm_zone = UTM zone
%
% DESCRIPTION:
%   Conversion from geodetic coordinates to planimetric coordinates (UTM WGS84).

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Laboratorio di Geomatica, Polo Regionale di Como,
%    Politecnico di Milano, Italy
% Portions of code based on "deg2utm" by Rafael Palacios, Universidad Pontificia Comillas
%    Madrid, Spain
%----------------------------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%----------------------------------------------------------------------------------------------

%number of input points
n = size(lat);
n = n(1,1);

if (n > 0)
    %pre-allocation
    M = zeros(n,1);
    north_south = zeros(n,1);
    utm_zone(n,:) = '60 X';
    NORTH = zeros(n,1);
    EAST = zeros(n,1);
    
    %conversion algorithm
    %UTM parameters
    EsMc = 500000; % false East
    
    %WGS84 ellipsoid parameters
    %semi-major equatorial axis [m]
    SemiEq = 6378137;
    %flattening f=(a-b)/a
    f = 1 / 298.25722356;
    
    %squared eccentricity (a^2-b^2)/a^2
    EccQ = 1 - (1 - f)^2;
    
    %contraction factor
    contr = 0.9996;
    
    SemiEq = SemiEq * contr;
    ecc4 = EccQ * EccQ;
    ecc6 = ecc4 * EccQ;
    ecc8 = ecc6 * EccQ;
    
    k0 = SemiEq * (EccQ / 4 + ecc4 * 3 / 64 + ecc6 * 5 / 256 + ecc8 * 175 / 16384);
    k = SemiEq - k0;
    k1 = SemiEq * (ecc4 * 13 / 96 + ecc6 * 59 / 384 + ecc8 * 1307 / 8192);
    k2 = SemiEq * (ecc6 * 61 / 484 + ecc8 * 609 / 2048);
    k3 = SemiEq * (ecc8 * 49561 / 322560);
    c1 = (EccQ * 5 - ecc4) / 6;
    c2 = (ecc4 * 104 - ecc6 * 45) / 120;
    c3 = ecc6 * 1237 / 1260;
    
    %Sines, cosines and latitude powers
    Elix = [lat lon];
    latsessadec = lat ./ pi .* 180;
    lonsessadec = lon ./ pi .* 180;
    
    fiSin(:,1) = sin(Elix(:,1));
    fiCos(:,1) = cos(Elix(:,1));
    fiSin2(:,1) = fiSin .* fiSin;
    fiSin4(:,1) = fiSin2 .* fiSin2;
    fiSin6(:,1) = fiSin4 .* fiSin2;
    
    %UTM zone finding
    for i = 1 : n

        M(i,1) = fix((180 + lonsessadec(i,1)) / 6) + 1;
        
        if latsessadec(i,1) >= 0
            north_south(i,1) = 1; %1 north, 0 south
        else
            north_south(i,1) = 0;
        end
        
        if     (latsessadec(i,1) < -72), letter='C';
        elseif (latsessadec(i,1) < -64), letter='D';
        elseif (latsessadec(i,1) < -56), letter='E';
        elseif (latsessadec(i,1) < -48), letter='F';
        elseif (latsessadec(i,1) < -40), letter='G';
        elseif (latsessadec(i,1) < -32), letter='H';
        elseif (latsessadec(i,1) < -24), letter='J';
        elseif (latsessadec(i,1) < -16), letter='K';
        elseif (latsessadec(i,1) <  -8), letter='L';
        elseif (latsessadec(i,1) <   0), letter='M';
        elseif (latsessadec(i,1) <   8), letter='N';
        elseif (latsessadec(i,1) <  16), letter='P';
        elseif (latsessadec(i,1) <  24), letter='Q';
        elseif (latsessadec(i,1) <  32), letter='R';
        elseif (latsessadec(i,1) <  40), letter='S';
        elseif (latsessadec(i,1) <  48), letter='T';
        elseif (latsessadec(i,1) <  56), letter='U';
        elseif (latsessadec(i,1) <  64), letter='V';
        elseif (latsessadec(i,1) <  72), letter='W';
        else                             letter='X';
        end
        
        utm_zone(i,:) = sprintf('%02d %c', M(i, 1), letter);
    end
    
    for i = 1 : n
        %Longitude of the central meridian
        LonMeridianoCentrale = -177+6 * (M(i,1) - 1);
        
        %Distance of the point from the central meridian
        %la_sd --> distance in decimal degrees
        %la    --> distance in radians
        la_sd = lonsessadec(i,1) - LonMeridianoCentrale;
        la = la_sd / 180 * pi;
        
        if la == 0
            laSin = 0;
            laCos = 1;
            laCot = 0;
        else
            laSin = sin(la);
            laCos = cos(la);
            laCot = laCos / laSin ;
        end
        
        %longitude with respect to central meridian
        laCot2 = laCot * laCot;
        
        %psi
        psi = Elix(i,1) - EccQ * fiSin(i,1) * fiCos(i,1) *(1 + c1 * fiSin2(i,1) + c2 * fiSin4(i,1) + c3 * fiSin6(i,1));
        psiSin = sin(psi);
        psiCos = cos(psi);
        psiTan = psiSin / psiCos ;
        psiSin2 = psiSin * psiSin ;
        
        %omega
        ome = atan(psiTan / laCos);
        
        %sigma
        if laSin ~= 0
            sigSin =laSin * psiCos;
            sig = asin(sigSin);
        else
            sigSin = 0;
            sig = 0;
        end
        
        sigSin2 = sigSin * sigSin;
        sigSin4 = sigSin2 * sigSin2;
        sigCos2 = 1 - sigSin2;
        sigCos4 = sigCos2 * sigCos2;
        
        %chi
        chi = sig / 2 + pi / 4;
        chiTan = tan(chi);
        chiLog = log(chiTan);
        
        %constants
        aa = psiSin * psiCos * laCos * (1 + sigSin2) / sigCos4;
        bb = sigSin * (sigCos2 - 2 * psiSin2) / sigCos4;
        
        if laCot ~= 0
            a1 = (psiSin2 - sigSin4 * laCot2)/ sigCos4;
            b1 = 2 * sigSin2 * psiSin * laCot / sigCos4;
        else
            a1 = psiSin2 / sigCos4 ;
            b1 = 0;
        end
        a2 = a1 * a1 - b1 * b1 ;
        b2 = 2 * a1 * b1 ;
        a3 = a1 * a2 - b1 * b2 ;
        b3 = a1 * b2 + b1 * a2 ;
        rr = k0 - a1 * k1 + a2 * k2 - a3 * k3;
        tt = b1 * k1 - b2 * k2 + b3 * k3;
        
        %X/Y coordinates
        xx = k * ome + aa * rr + bb * tt;
        yy = k * chiLog + bb * rr - aa * tt;
        
        %North and East
        if north_south(i,1) == 1
            NoEq = 0;
        else
            NoEq = 10000000;
        end
        
        NORTH(i,1) = NoEq + xx;
        EAST(i,1) = EsMc + yy;
    end
    
else
    utm_zone = [];
    NORTH = [];
    EAST = [];
end
