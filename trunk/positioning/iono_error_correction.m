function [corr] = iono_error_correction(lat, lon, az, el, T, ionoparams)

% SYNTAX:
%   [corr] = iono_error_correction(lat, lon, az, el, T, ionoparams);
%
% INPUT:
%   lat = receiver latitude    [degrees]
%   lon = receiver longitude   [degrees]
%   az  = satellite azimuth    [degrees]
%   el  = satellite elevation  [degrees]
%   T   = time
%   ionoparams = ionospheric correction parameters
%
% OUTPUT:
%   corr = ionospheric error correction
%
% DESCRIPTION:
%   Computation of the pseudorange correction due to ionospheric refraction.
%   Klobuchar algorithm.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.1 beta
%
% Copyright (C) 2009-2012 Mirko Reguzzoni, Eugenio Realini
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

%global variable
global v_light

%----------------------------------------------------------------------------------------------

if ((nargin > 5) & (sum(abs(ionoparams)) > 0))

    %ionospheric parameters
    a0 = ionoparams(1);
    a1 = ionoparams(2);
    a2 = ionoparams(3);
    a3 = ionoparams(4);
    b0 = ionoparams(5);
    b1 = ionoparams(6);
    b2 = ionoparams(7);
    b3 = ionoparams(8);

    %elevation from 0 to 90 degrees
    el = abs(el);

    %conversion to semicircles
    lat = lat / 180;
    lon = lon / 180;
    az = az / 180;
    el = el / 180;

    f = 1 + 16*(0.53-el).^3;

    psi = (0.0137 ./ (el+0.11)) - 0.022;

    phi = lat + psi .* cos(az*pi);
    phi(phi > 0.416)  =  0.416;
    phi(phi < -0.416) = -0.416;

    lambda = lon + ((psi.*sin(az*pi)) ./ cos(phi*pi));

    ro = phi + 0.064*cos((lambda-1.617)*pi);

    t = lambda*43200 + T;

    for i = 1 : length(T)
        while (t(i) >= 86400)
            t(i) = t(i)-86400;
        end
        while (t(i) < 0)
            t(i) = t(i)+86400;
        end
    end

    % index = find(t >= 86400);
    % while ~isempty(index)
    %     t(index) = t(index)-86400;
    %     index = find(t >= 86400);
    % end

    % index = find(t < 0);
    % while ~isempty(index)
    %     t(index) = t(index)+86400;
    %     index = find(t < 0);
    % end

    a = a0 + a1*ro + a2*ro.^2 + a3*ro.^3;
    a(a < 0) = 0;

    p = b0 + b1*ro + b2*ro.^2 + b3*ro.^3;
    p(p < 72000) = 72000;

    x = (2*pi*(t-50400)) ./ p;

    %ionospheric error
    index = find(abs(x) < 1.57);
    corr(index,1) = v_light * f(index) .* (5e-9 + a(index) .* (1 - (x(index).^2)/2 + (x(index).^4)/24));

    index = find(abs(x) >= 1.57);
    corr(index,1) = v_light * f(index) .* 5e-9;

else
    %a simplified model could be used
    corr = zeros(size(el));
end