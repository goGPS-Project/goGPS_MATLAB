function nmeastring = NMEA_RMC_gen(pos_R, time)

% SYNTAX:
%   nmeastring = NMEA_RMC_gen(pos_R, time);
%
% INPUT:
%   pos_R = estimated ROVER position (X,Y,Z)
%   time = GPS time of measurements
%
% OUTPUT:
%   nmeastring = $GPRMC sentence (NMEA)
%
% DESCRIPTION:
%   Returns a $GPRMC sentence in NMEA 0183 format.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
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

%date
if (nargin > 1)
    if (size(time) == [1,6]) %#ok<BDSCA>
        if (time(1) < 100)
            time(1) = four_digit_year(time(1));
        end
        date = datevec(gps2utc(datenum(time)));
    else
        [week, sow] = time2weektow(time);
        date = datevec(gps2utc(datenum(gps2date(week, sow))));
        %date = datevec((time / (3600*24)) - fix(time / (3600*24)));
    end
else
    %get current date
    date = clock;
end

%approximate coordinates X Y Z
X = pos_R(1);
Y = pos_R(2);
Z = pos_R(3);

%conversion from cartesian to geodetic coordinates
[phi, lam] = cart2geod(X, Y, Z);

%conversion from radians to degrees
lam = abs(lam*180/pi);
phi = abs(phi*180/pi);

%-----------------------------------------------------------------------------------------------
% FORMAT COORDINATES ACCORDING TO NMEA STANDARD
%-----------------------------------------------------------------------------------------------

%longitude
lam_deg = floor(lam);
lam_min = (lam - lam_deg) * 60;
lam_nmea = lam_deg * 100 + lam_min;

[lam_str] = sprintf('%d', lam_deg);
if (length(lam_str) == 1 )
    [lam_nmea] = sprintf('00%.8f', lam_nmea);
elseif (length(lam_str) == 2 )
    [lam_nmea] = sprintf('0%.8f', lam_nmea);
else
    [lam_nmea] = sprintf('%.8f', lam_nmea);
end

%latitude
phi_deg = floor(phi);
phi_min = (phi - phi_deg) * 60;
phi_nmea = phi_deg * 100 + phi_min;

[phi_str] = sprintf('%d', phi_deg);
if (length(phi_str) == 1 )
    [phi_nmea] = sprintf('0%.8f', phi_nmea);
else
    [phi_nmea] = sprintf('%.8f', phi_nmea);
end

%emisphere definition
if (lam >= 0)
    emi_EW = 'E';
else
    emi_EW = 'W';
end

if (phi >= 0)
    emi_NS = 'N';
else
    emi_NS = 'S';
end

%-----------------------------------------------------------------------------------------------
% OTHER NMEA PARAMETERS
%-----------------------------------------------------------------------------------------------

status = 'A';  %V = navigation receiver warning, A = Valid
speed = ''; %speed over ground (knots), to be changed

%-----------------------------------------------------------------------------------------------
% FORMAT DATA
%-----------------------------------------------------------------------------------------------

hour = sprintf('%02d', date(1,4));
minute = sprintf('%02d', date(1,5));
second = sprintf('%05.2f', date(1,6));

% [null, nchar] = size(hour); %#ok<ASGLU>
% if (nchar == 1)
%     [hour] = sprintf('0%s',hour);
% end
%
% [null, nchar] = size(minute); %#ok<ASGLU>
% if (nchar == 1)
%     [minute] = sprintf('0%s',minute);
% end
%
% [null, nchar] = size(second); %#ok<ASGLU>
% if (nchar == 1)
%     [second] = sprintf('0%s',second);
% end
%
% decsec = '.00';

day = sprintf('%02d', date(1,3));
month = sprintf('%02d', date(1,2));
year = sprintf('%02d', two_digit_year(date(1,1)));

% [null, nchar] = size(day); %#ok<ASGLU>
% if (nchar == 1)
%     [day] = sprintf('0%s',day);
% end
%
% [null, nchar] = size(month); %#ok<ASGLU>
% if (nchar == 1)
%     [month] = sprintf('0%s',month);
% end
%
% [null, nchar] = size(year); %#ok<ASGLU>
% if (nchar == 1)
%     [year] = sprintf('0%s',year);
% end

%-----------------------------------------------------------------------------------------------
% COMPOSITION OF THE NMEA SENTENCE
%-----------------------------------------------------------------------------------------------

nmeastring = sprintf('$GPRMC,%s%s%s,%c,%s,%c,%s,%c,%s,,%s%s%s,,,,',hour,minute,second,status,phi_nmea,emi_NS,lam_nmea,emi_EW,speed,day,month,year);

%-----------------------------------------------------------------------------------------------
% CHECKSUM COMPUTATION
%-----------------------------------------------------------------------------------------------

checksum = NMEA_checksum(nmeastring);
nmeastring = [nmeastring '*' checksum];
