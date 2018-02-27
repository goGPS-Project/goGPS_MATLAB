%   CLASS Atmosphere
% =========================================================================
%
% DESCRIPTION
%   Class to store static method to compute atmospheric corrections using
%   different models
%
% EXAMPLE
%   ls = LS();
%


%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.2 beta 1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by: Giulio Tagliaferro
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
classdef Atmosphere < handle
    properties
    end
    methods (Static)
        %-----------------------------------------------------------
        % TROPO
        %-----------------------------------------------------------
        function [delay] = saastamoinen_model(lat, lon, h, el)
            % SYNTAX:
            %   [delay] = Atmosphere.tropo_error_correction(lat, lon, h, el);
            %
            % INPUT:
            %   time_rx = receiver reception time
            %   lat = receiver latitude          [degrees]
            %   lon = receiver longitude         [degrees]
            %   h  = receiver ellipsoidal height [meters]  [nx1]
            %   el = satellite elevation         [degrees] [nx1]
            %
            % OUTPUT:
            %   corr = tropospheric error correction
            %
            % DESCRIPTION:
            %   Computation of the pseudorange correction due to tropospheric refraction.
            %   Saastamoinen algorithm using standard atmosphere accounting for
            %   height gradient for temperature pressure and humidity.
            %   --> multi epoch for static receiver
            gs = Go_State.getInstance;
            geoid = gs.getRefGeoid();
            
            %Saastamoinen model requires (positive) orthometric height
            if geoid.ncols > 0
                %geoid ondulation interpolation
                undu = getOrthometricCorr(lon, lat);
                h(undu > -300) = h(undu > -300) - undu(undu > -300);
            end
            h(h < 0) = 0;
            
            if (h < 5000)
                
                %conversion to radians
                el = abs(el) * pi/180;
                
                %Standard atmosphere - Berg, 1948 (Bernese)
                %pressure [mbar]
                Pr = goGNSS.STD_PRES;
                %temperature [K]
                Tr = goGNSS.STD_TEMP;
                %humidity [%]
                Hr = goGNSS.STD_HUMI;
                
                P = Pr * (1-0.0000226*h).^5.225;
                T = Tr - 0.0065*h;
                H = Hr * exp(-0.0006396*h);
                
                %----------------------------------------------------------------------
                
                %linear interpolation
                h_a = [0; 500; 1000; 1500; 2000; 2500; 3000; 4000; 5000];
                B_a = [1.156; 1.079; 1.006; 0.938; 0.874; 0.813; 0.757; 0.654; 0.563];
                
                t = zeros(length(T),1);
                B = zeros(length(T),1);
                
                for i = 1 : length(T)
                    
                    d = h_a - h(i);
                    [~, j] = min(abs(d));
                    if (d(j) > 0)
                        index = [j-1; j];
                    else
                        index = [j; j+1];
                    end
                    
                    t(i) = (h(i) - h_a(index(1))) ./ (h_a(index(2)) - h_a(index(1)));
                    B(i) = (1-t(i))*B_a(index(1)) + t(i)*B_a(index(2));
                end
                
                %----------------------------------------------------------------------
                
                e = 0.01 * H .* exp(-37.2465 + 0.213166*T - 0.000256908*T.^2);
                
                %tropospheric delay
                delay = ((0.002277 ./ sin(el)) .* (P - (B ./ max(0.01,(tan(el)).^2))) + (0.002277 ./ sin(el)) .* (1255./T + 0.05) .* e); % max to eliminate numeric instability near 0
                
            else
                delay = zeros(size(el));
            end
        end
        
        function [delay] = saastamoinen_model_GPT(time_rx, lat, lon, h, el)
            % SYNTAX:
            %   [delay] = Atmosphere.saastamoinen_model_GPT(time_rx, lat, lon, h, el)
            %
            % INPUT:
            %   time_rx = receiver reception time
            %   lat = receiver latitude          [degrees]
            %   lon = receiver longitude         [degrees]
            %   h  = receiver ellipsoidal height [meters]
            %   el = satellite elevation         [degrees] [nx1]
            %
            % OUTPUT:
            %   corr = tropospheric error correction
            %
            % DESCRIPTION:
            %   Computation of the pseudorange correction due to tropospheric refraction.
            %   Saastamoinen algorithm using P T from Global Pressure and Temperature
            %   (GPT), and H from standard atmosphere accounting for humidity height gradient.
            %   --> single epoch
            global zero_time geoid
            
            delay = zeros(size(el));
            
            [week, sow] = time2weektow(time_rx);
            date = gps2date(week, sow);
            [~, mjd] = date2jd(date);
            
            [pres, temp, undu] = gpt(mjd, lat*pi/180, lon*pi/180, h);
            if (exist('geoid','var') && isfield(geoid,'ncols') && geoid.ncols ~= 0)
                %geoid ondulation interpolation
                undu = getOrthometricCorr(lon, lat);
            end
            t_h = h;
            t_h(undu > -300) = t_h(undu > -300) - undu(undu > -300);
            ZHD_R = saast_dry(pres, t_h, lat);
            ZWD_R = saast_wet(temp, goGNSS.STD_HUMI, t_h);
            
            for s = 1 : length(el)
                [gmfh_R, gmfw_R] = gmf_f_hu(mjd, lat*pi/180, lon*pi/180, h, (90-el(s,1))*pi/180);
                delay(s,1) = gmfh_R*ZHD_R + gmfw_R*ZWD_R;
            end
        end
        %-----------------------------------------------------------
        % IONO
        %-----------------------------------------------------------
        function [delay] = klobuchar_model(lat, lon, az, el, sow, ionoparams)
            % SYNTAX:
            %   [delay] = Atmosphere. klobuchar_model(lat, lon, az, el, sow, ionoparams)
            %
            % INPUT:
            %   lat = receiver latitude          [degrees] [nx1]
            %   lon = receiver longitude         [degrees] [nx1]
            %   az  = satellite azimuth          [degrees] [nx1]
            %   el = satellite elevation         [degrees] [nx1]
            %   sow = second of week                       [nx1]
            % OUTPUT:
            %   corr = tropospheric error correction
            %
            % DESCRIPTION:
            %   Computation of the pseudorange correction due to ionosphere.
            %   --> multiple epoch for both static and moving target
            %-------------------------------------------------------------------------------
            % KLOBUCHAR MODEL
            %
            % Algorithm from Leick, A. (2004) "GPS Satellite Surveying - 2nd Edition"
            % John Wiley & Sons, Inc., New York, pp. 301-303)
            %-------------------------------------------------------------------------------
            %initialization
            delay = zeros(size(el));
            
            
            
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
            
            t = lambda*43200 + sow;
            t = mod(t,86400);
            
            
            a = a0 + a1*ro + a2*ro.^2 + a3*ro.^3;
            a(a < 0) = 0;
            
            p = b0 + b1*ro + b2*ro.^2 + b3*ro.^3;
            p(p < 72000) = 72000;
            
            x = (2*pi*(t-50400)) ./ p;
            
            %ionospheric delay
            index = find(abs(x) < 1.57);
            delay(index,1) = goGNSS.V_LIGHT * f(index) .* (5e-9 + a(index) .* (1 - (x(index).^2)/2 + (x(index).^4)/24));
            
            index = find(abs(x) >= 1.57);
            delay(index,1) = goGNSS.V_LIGHT * f(index) .* 5e-9;
        end
        
        function [lat_pp, lon_pp, iono_mf, k] = getPiercePoint(lat_rad, lon_rad, h_ortho, az_rad, el_rad, thin_shell_height, rcm)
            % Get the pierce point
            % INPUT:
            %   lat_rad             latitude of the receiver           [rad]
            %   lon_rad             longitude of the receiver          [rad]
            %   h_ortho             orthometric height of the receiver [m]
            %   az_rad              azimuth of the satellites          [rad]
            %   el_rad              elevation of the satellites        [rad]
            %   thin_shell_height   height of the pierce point         [m]
            %   rcm                 meridian radius curvature <optional>
            %
            % OUTPUT
            %   latpp               latitude pierce point [rad]
            %   lonpp               longitude pierce point [rad]
            %   iono_mf             iono mapping function
            %   k                   iono k factor ????
            % 
            % SYNTAX: 
            %   [latpp, lonpp, mfpp, k] = getPiercePoint(lat_rad, lon_rad, h_ortho, az_rad, el_rad, thin_shell_height, <rcm>)
            
            % Get radius of curvature at lat
            if nargin < 7
                rcm = getMeridianRadiusCurvature(lat_rad);
            end
            
            input_size = (size(az_rad));
            az_rad = az_rad(:);
            el_rad = el_rad(:);
            
            k = ((rcm + h_ortho)/((rcm + h_ortho) + thin_shell_height)) * cos(el_rad);
            psi_pp = (pi/2) - el_rad - asin(k);
            
            %set azimuth from -180 to 180
            az_rad = mod((az_rad+pi),2*pi)-pi;
            
            %latitude of the ionosphere piercing point
            lat_pp = asin(sin(lat_rad) * cos(psi_pp) + cos(lat_rad) * sin(psi_pp) .* cos(az_rad));
            
            %longitude of the ionosphere piercing point
            id_hl = ((lat_pp >  70*pi/180) & (tan(psi_pp).*cos(az_rad)      > tan((pi/2) - lat_rad))) | ...
                ((lat_pp < -70*pi/180) & (tan(psi_pp).*cos(az_rad + pi) > tan((pi/2) + lat_rad)));
            
            lon_pp = zeros(size(az_rad));
            lon_pp(id_hl) = lon_rad + pi - asin(sin(psi_pp(id_hl)) .* sin(az_rad(id_hl)) ./ cos(lat_pp(id_hl)));

            lon_pp(~id_hl) = lon_rad + asin(sin(psi_pp(~id_hl)) .* sin(az_rad(~id_hl)) ./ cos(lat_pp((~id_hl))));
            
            % using thin shell layer mapping function (Handbook of Global
            % Navigation System pp 185)
            if nargout > 2
                iono_mf = (1-(k)^2)^(-1/2);
            end
            
            if nargout > 3
                k = [cos(az_rad).*cos(el_rad);
                    -sin(az_rad).*cos(el_rad);
                    sin(el_rad)];
                % go to global system
                R = [-sin(lat_rad) cos(lon_rad) 0;
                    -sin(lat_rad)*cos(lon_rad) -sin(lat_rad)*sin(lon_rad) cos(lat_rad);
                    +cos(lat_rad)*cos(lon_rad) +cos(lat_rad)*sin(lon_rad) sin(lat_rad)];
                [k] = R'*k;
            end
            
            lat_pp = reshape(lat_pp, input_size(1), input_size(2));
            lon_pp = reshape(lon_pp, input_size(1), input_size(2));
        end
    end
end