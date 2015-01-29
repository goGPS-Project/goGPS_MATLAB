function [corr] = iono_error_correction(lat, lon, az, el, time_rx, ionoparams, sbas)

% SYNTAX:
%   [corr] = iono_error_correction(lat, lon, az, el, time_rx, ionoparams, sbas);
%
% INPUT:
%   lat = receiver latitude    [degrees]
%   lon = receiver longitude   [degrees]
%   az  = satellite azimuth    [degrees]
%   el  = satellite elevation  [degrees]
%   time_rx    = receiver reception time
%   ionoparams = ionospheric correction parameters
%   sbas = SBAS corrections
%
% OUTPUT:
%   corr = ionospheric error correction
%
% DESCRIPTION:
%   Computation of the pseudorange correction due to ionospheric delay.
%   Klobuchar model or SBAS ionosphere interpolation.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Laboratorio di Geomatica, Polo Regionale di Como,
%    Politecnico di Milano, Italy
% Portions of code contributed by Giuliano Sironi, 2011
% Portions of code contributed by Antonio Herrera Olmo, 2012
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

global iono_model

v_light = goGNSS.V_LIGHT;

% If in ATM_model section in config file iono is set to 0 it does not use ionospheric model
% [ATM_model]
% iono=0
if (iono_model == 0)
    corr = zeros(size(el));
else
    %initialization
    corr = zeros(size(el));

    %if ionosphere parameters are available and SBAS corrections are disabled/not available
    if ((nargin == 6) & (sum(abs(ionoparams)) > 0)) | ...
	    ((nargin >  6) & (sum(abs(ionoparams)) > 0)  & (isempty(sbas)))
	
	%apply Klobuchar ionosphere model
	corr = klobuchar_model(lat, lon, az, el, time_rx, ionoparams);
	
	%if SBAS corrections are available (and requested by the user)
    elseif ((nargin > 6) & (~isempty(sbas)))
	
	%apply SBAS interpolated ionospheric delay (where possible)
	corr = sbas_iono_interp(lat, lon, az, el, sbas);
	
	%detect if some satellites could not be corrected by SBAS
	not_corr = isnan(corr);
	
	if (any(not_corr))
	    %apply Klobuchar ionosphere model where it was not possible to apply
	    % SBAS corrections
	    corr(not_corr) = klobuchar_model(lat, lon, az(not_corr), el(not_corr), time_rx, ionoparams);
	end
	
    else
	%a simplified model could be used
    end
end

% -------------------------------------------------------------------------
% End of function - start nested function declaration
% -------------------------------------------------------------------------

    function [delay] = klobuchar_model(lat, lon, az, el, time_rx, ionoparams)
	
	%initialization
	delay = zeros(size(el));
	
	%-------------------------------------------------------------------------------
	% KLOBUCHAR MODEL
	%
	% Algorithm taken from Leick, A. (2004) "GPS Satellite Surveying - 2nd Edition"
	% John Wiley & Sons, Inc., New York, pp. 301-303)
	%-------------------------------------------------------------------------------
	
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
	
	t = lambda*43200 + time_rx;
	t = mod(t,86400);
	
	% for i = 1 : length(time_rx)
	%    while (t(i) >= 86400)
	%        t(i) = t(i)-86400;
	%    end
	%    while (t(i) < 0)
	%        t(i) = t(i)+86400;
	%    end
	%end
	
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
	
	%ionospheric delay
	index = find(abs(x) < 1.57);
	delay(index,1) = v_light * f(index) .* (5e-9 + a(index) .* (1 - (x(index).^2)/2 + (x(index).^4)/24));
	
	index = find(abs(x) >= 1.57);
	delay(index,1) = v_light * f(index) .* 5e-9;
    end

    function [delay] = sbas_iono_interp(lat, lon, az, el, sbas)
	
	%initialization
	delay = NaN(size(el));
	
	%-------------------------------------------------------------------------------
	% SBAS IONOSPHERE INTERPOLATION
	%-------------------------------------------------------------------------------
	
	lat = lat * pi/180; %rad
	lon = lon * pi/180; %rad
	az  = az  * pi/180; %rad
	el  = el  * pi/180; %rad
	
	for i = 1 : length(az)
	    
	    %ionosphere pierce point coordinates and slant factor
	    [latpp, lonpp, fpp] = iono_pierce_point(lat, lon, az(i), el(i));
	    
	    %find the nodes of the cell that contains the ionosphere pierce point
	    [igp4, tv] = sel_igp(latpp, lonpp, sbas.igp, sbas.lat_igp, sbas.lon_igp);
	    
	    %the ionospheric vertical delay can be interpolated only if the
	    %ionosphere pierce point falls within a valid SBAS ionosphere grid
	    %cell (rectangular or triangular).
	    if (sum(isfinite(igp4)) >= 3)
		
		%interpolate the ionospheric vertical delay at the piercing point
		[ivd_pp] = interp_ivd(igp4, sbas.igp, sbas.ivd, latpp, lonpp, tv); %m
		
		%ionospheric delay
		delay(i,1) = fpp * ivd_pp; % m
	    end
	end
    end

end