function [ionocorr] = err_iono(ionoparams, Lat, Lon, Az, El, T)

% SYNTAX:
%   [ionocorr] = err_iono(ionoparams, Lat, Lon, Az, El, T);
%
% INPUT:
%   ionoparams = ionospheric correction parameters
%   Lat = rover latitude
%   Lon = rover longitude
%   Az = satellite azimuth
%   El = satellite elevation
%   T = time
%
% OUTPUT:
%   ionocorr = ionospheric correction
%
% DESCRIPTION:
%   Computation of the pseudorange correction due to ionospheric refraction.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.2.0 beta
%
% Copyright (C) 2009-2011 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Laboratorio di Geomatica, Polo Regionale di Como,
%    Politecnico di Milano, Italy
%
% Algorithm taken from Leick, A. (2004) "GPS Satellite Surveying - 2nd Edition"
% John Wiley & Sons, Inc., New York, pp. 301-303
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

%variable initialization
global v_light

if (sum(abs(ionoparams)) == 0)

    ionocorr = 0;
else
    
    %iono parameters
    a0 = ionoparams(1);
    a1 = ionoparams(2);
    a2 = ionoparams(3);
    a3 = ionoparams(4);
    b0 = ionoparams(5);
    b1 = ionoparams(6);
    b2 = ionoparams(7);
    b3 = ionoparams(8);
    El=abs(El);
    
    %conversion to semicircles
    Lat = Lat / 180;
    Lon = Lon / 180;
    Az = Az / 180;
    El = El / 180;
    
    %Klobuchar algorithm
    f=1+16*(0.53-El)^3;
    
    psi=(0.0137/(El+0.11))-0.022;
    
    phi=Lat+psi*cos(Az*pi);
    
    if (phi>0.416)
        phi=0.416;
    end
    
    if (phi<-0.416)
        phi=-0.416;
    end

    lambda=Lon+((psi*sin(Az*pi))/cos(phi*pi));

    ro=phi+0.064*cos((lambda-1.617)*pi);
    
    t=lambda*43200+T;
    
    while (t>=86400)
        t=t-86400;
    end
    
    while (t<0)
        t=t+86400;
    end
    
    p=b0+b1*ro+b2*ro^2+b3*ro^3;
    
    if (p<72000)
        p=72000;
    end
    
    a=a0+a1*ro+a2*ro^2+a3*ro^3;
    
    if (a<0)
        a=0;
    end
    
    x=(2*pi*(t-50400))/p;
    
    %ionospheric error
    if (abs(x)<1.57)
        ionocorr = v_light * f * (5e-9+a*(1-(x^2)/2+(x^4)/24));
    else
        ionocorr = v_light * f * 5e-9;
    end
    
end