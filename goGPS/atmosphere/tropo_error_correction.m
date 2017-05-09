function [corr] = tropo_error_correction(time_rx, lat, lon, h, el)

% SYNTAX:
%   [corr] = tropo_error_correction(time_rx, lat, lon, h, el);
%
% INPUT:
%   time_rx = receiver reception time
%   lat = receiver latitude          [degrees]
%   lon = receiver longitude         [degrees]
%   h  = receiver ellipsoidal height [meters]
%   el = satellite elevation         [degrees]
%
% OUTPUT:
%   corr = tropospheric error correction
%
% DESCRIPTION:
%   Computation of the pseudorange correction due to tropospheric refraction.
%   Saastamoinen algorithm.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.1 beta 2
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
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

global tropo_model

switch tropo_model
    case 0 %no model
        corr = zeros(size(el));
    case 1 %Saastamoinen model (with standard atmosphere parameters)
        corr = saastamoinen_model(lat, lon, h, el);
    case 2 %Saastamoinen model (with Global Pressure Temperature model)
        corr = saastamoinen_model_GPT(time_rx, lat, lon, h, el);
end

% -------------------------------------------------------------------------
% End of function - start nested function declaration
% -------------------------------------------------------------------------

    function [delay] = saastamoinen_model(lat, lon, h, el)
        global geoid

        %Saastamoinen model requires (positive) orthometric height
        if (exist('geoid','var') && isfield(geoid,'ncols') && geoid.ncols ~= 0)
            %geoid ondulation interpolation
            undu = grid_bilin_interp(lon, lat, geoid.grid, geoid.ncols, geoid.nrows, geoid.cellsize, geoid.Xll, geoid.Yll, -9999);
            h = h - undu;
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
            delay = ((0.002277 ./ sin(el)) .* (P - (B ./ (tan(el)).^2)) + (0.002277 ./ sin(el)) .* (1255./T + 0.05) .* e);
        else
            delay = zeros(size(el));
        end
    end

    function [delay] = saastamoinen_model_GPT(time_rx, lat, lon, h, el)
        global zero_time geoid

        delay = zeros(size(el));

        [week, sow] = time2weektow(time_rx + zero_time);
        date = gps2date(week, sow);
        [~, mjd] = date2jd(date);

        [pres, temp, undu] = gpt(mjd, lat*pi/180, lon*pi/180, h);
        if (exist('geoid','var') && isfield(geoid,'ncols') && geoid.ncols ~= 0)
            %geoid ondulation interpolation
            undu = grid_bilin_interp(lon, lat, geoid.grid, geoid.ncols, geoid.nrows, geoid.cellsize, geoid.Xll, geoid.Yll, -9999);
        end
        ZHD_R = saast_dry(pres, h - undu, lat);
        ZWD_R = saast_wet(temp, goGNSS.STD_HUMI, h - undu);

        for s = 1 : length(el)
            [gmfh_R, gmfw_R] = gmf_f_hu(mjd, lat*pi/180, lon*pi/180, h, (90-el(s,1))*pi/180);
            delay(s,1) = gmfh_R*ZHD_R + gmfw_R*ZWD_R;
        end
    end
end
