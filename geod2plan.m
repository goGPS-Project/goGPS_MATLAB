function [EST, NORD, M, nord_sud] = geod2plan(lat, lon)

% SYNTAX:
%   [EST, NORD, M, nord_sud] = geod2plan(lat, lon);
%
% INPUT:
%   lat = latitude [rad]
%   lon = longitude [rad]
%
% OUTPUT:
%   EST = EAST coordinate [m]
%   NORD = NORTH coordinate [m]
%   M = UTM zone
%   nord_sud = (north=1), (south=0)
%
% DESCRIPTION:
%   Conversion from geodetic coordinates to planimetric coordinates (UTM WGS84).

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 alpha
%
% Copyright (C) 2009 Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
%
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
n=size(lat);
n=n(1,1);

EST = [];
NORD = [];

%conversion algorithm
%UTM parameters
EsMc=500000; % false East

%WGS84 ellipsoid parameters
%semi-major equatorial axis [m]
SemiEq=6378137;
%flattening f=(a-b)/a 
f=1/298.25722356;

%squared eccentricity (a^2-b^2)/a^2
EccQ=1-(1-f)^2;

%contraction factor
contr=0.9996;

SemiEq=SemiEq*contr;
ecc4=EccQ*EccQ;
ecc6=ecc4*EccQ;
ecc8=ecc6*EccQ;

k0=SemiEq*(EccQ/4+ecc4*3/64+ecc6*5/256+ecc8*175/16384);
k=SemiEq-k0;
k1=SemiEq*(ecc4*13/96+ecc6*59/384+ecc8*1307/8192);
k2=SemiEq*(ecc6*61/484+ecc8*609/2048);
k3=SemiEq*(ecc8*49561/322560);
c1=(EccQ*5-ecc4)/6;
c2=(ecc4*104-ecc6*45)/120;
c3=ecc6*1237/1260;

%Sines, cosines and latitude powers
Elix=[lat lon];
latsessadec=lat./pi.*180;
lonsessadec=lon./pi.*180;

fiSin(:,1)=sin(Elix(:,1));
fiCos(:,1)=cos(Elix(:,1));
fiSin2(:,1)=fiSin.*fiSin;
fiSin4(:,1)=fiSin2.*fiSin2;
fiSin6(:,1)=fiSin4.*fiSin2;

%UTM zone finding
for i=1:n
    if lonsessadec(i,1)>=0 & lonsessadec(i,1)<=180
       M(i,1)=fix((180+lonsessadec(i,1))/6)+1;
    else
       M(i,1)=fix((180-lonsessadec(i,1))/6)+1;
    end

    if latsessadec(i,1)>=0
       nord_sud(i,1)=1; %1 nord, 0 sud
    else   
       nord_sud(i,1)=0;
    end
end

%Longitude with respect to the central meridian
for i=1:n
    LonMeridianoCentrale(i,1)=-177+6*(M(i,1)-1);
end

%Distance of the point from the central meridian
%la_sd --> distance in decimal degrees
%la    --> distance in radians
for i=1:n
    la_sd(i,1)=lonsessadec(i,1)-LonMeridianoCentrale(i,1);
    la(i,1)=la_sd(i,1)/180*pi;
end

for i=1:n
    if la(i,1)==0
       laSin(i,1)=0;
       laCos(i,1)=1;
       laCot(i,1)=0;
    else
       laSin(i,1)=sin(la(i,1));
       laCos(i,1)=cos(la(i,1));
       laCot(i,1)=laCos(i,1)/laSin(i,1);
    end
end

%longitude with respect to central meridian
for i=1:n
    laCot2(i,1)=laCot(i,1)*laCot(i,1);
end

%psi
for i=1:n
    psi(i,1)=Elix(i,1)-EccQ*fiSin(i,1)*fiCos(i,1)*(1+c1*fiSin2(i,1)+c2*fiSin4(i,1)+c3*fiSin6(i,1));
    psiSin(i,1)=sin(psi(i,1));
    psiCos(i,1)=cos(psi(i,1));
    psiTan(i,1)=psiSin(i,1)/psiCos(i,1);
    psiSin2(i,1)=psiSin(i,1)*psiSin(i,1);
end

%omega
for i=1:n
    ome(i,1)=atan(psiTan(i,1)/laCos(i,1));
end

%sigma
for i=1:n
    if laSin(i,1)~=0
       sigSin(i,1)=laSin(i,1)*psiCos(i,1);
       sig(i,1)=asin(sigSin(i,1));
    else
       sigSin(i,1)=0;
       sig(i,1)=0;
    end

    sigSin2(i,1)=sigSin(i,1)*sigSin(i,1);
    sigSin4(i,1)=sigSin2(i,1)*sigSin2(i,1);
    sigCos2(i,1)=1-sigSin2(i,1);
    sigCos4(i,1)=sigCos2(i,1)*sigCos2(i,1);
end

%chi
for i=1:n
    chi(i,1)=sig(i,1)/2+pi/4;
    chiTan(i,1)=tan(chi(i,1));
    chiLog(i,1)=log(chiTan(i,1));
end

%constants
for i=1:n
    aa(i,1)=psiSin(i,1)*psiCos(i,1)*laCos(i,1)*(1+sigSin2(i,1))/sigCos4(i,1);
    bb(i,1)=sigSin(i,1)*(sigCos2(i,1)-2*psiSin2(i,1))/sigCos4(i,1);

    if laCot(i,1)~=0
       a1(i,1)=(psiSin2(i,1)-sigSin4(i,1)*laCot2(i,1))/sigCos4(i,1);
       b1(i,1)=2*sigSin2(i,1)*psiSin(i,1)*laCot(i,1)/sigCos4(i,1);
    else
       a1(i,1)=psiSin2(i,1)/sigCos4(i,1);
       b1(i,1)=0;
    end
    a2(i,1)=a1(i,1)*a1(i,1)-b1(i,1)*b1(i,1);
    b2(i,1)=2*a1(i,1)*b1(i,1);
    a3(i,1)=a1(i,1)*a2(i,1)-b1(i,1)*b2(i,1);
    b3(i,1)=a1(i,1)*b2(i,1)+b1(i,1)*a2(i,1);
    rr(i,1)=k0-a1(i,1)*k1+a2(i,1)*k2-a3(i,1)*k3;
    tt(i,1)=b1(i,1)*k1-b2(i,1)*k2+b3(i,1)*k3;
end

%X/Y coordinates
for i=1:n
    xx(i,1)=k*ome(i,1)+aa(i,1)*rr(i,1)+bb(i,1)*tt(i,1);
    yy(i,1)=k*chiLog(i,1)+bb(i,1)*rr(i,1)-aa(i,1)*tt(i,1);
end

%North and East
for i=1:n
    if nord_sud(i,1)==1
        NoEq(i,1)=0;
    else
        NoEq(i,1)=10000000;
    end

    NORD(i,1)=NoEq(i,1)+xx(i,1);
    EST(i,1)=EsMc+yy(i,1);
end
