function nmeastring = NMEA_GGA_gen(pos_R, nsat, time, HDOP)

% SYNTAX:
%   nmeastring = NMEA_GGA_gen(pos_R, nsat, time, HDOP);
%
% INPUT:
%   pos_R = estimated ROVER position (X,Y,Z)
%   nsat = number of visible satellites
%   time = GPS time of measurements
%   HDOP = horizontal dilution of precision
%
% OUTPUT:
%   nmeastring = $GPGGA sentence (NMEA)
%
% DESCRIPTION:
%   Returns a $GPGGA sentence in NMEA 0183 format.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 beta
%
% Copyright (C) 2009-2010 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Graduate School for Creative Cities, Osaka City University, Japan
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

if (nargin > 2)
    date = datevec((time / (3600*24)) - fix(time / (3600*24)));
else
    %get current date
    date = clock;
end

%approximate coordinates X Y Z
X = pos_R(1);
Y = pos_R(2);
Z = pos_R(3);

%conversion from cartesian to geodetic coordinates
[phi, lam, h] = cart2geod(X, Y, Z);

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

if (length(num2str(lam_deg)) == 1 )
    [lam_nmea] = sprintf('00%s', num2str(lam_nmea, '%.8f'));
elseif (length(num2str(lam_deg)) == 2 )
    [lam_nmea] = sprintf('0%s', num2str(lam_nmea, '%.8f'));
else
    lam_nmea = num2str(lam_nmea,'%.8f');
end

%latitude
phi_deg = floor(phi);
phi_min = (phi - phi_deg) * 60;
phi_nmea = phi_deg * 100 + phi_min;

if (length(num2str(phi_deg)) == 1 )
    [phi_nmea] = sprintf('0%s',num2str(phi_nmea,'%.8f'));
else
    phi_nmea = num2str(phi_nmea,'%.8f');
end

%height
if (length(num2str(floor(h))) == 1)
    [h] = sprintf('00%s',num2str(h,'%.3f'));
elseif (length(num2str(floor(h))) == 2)
    [h] = sprintf('0%s',num2str(h,'%.3f'));
else
    h = num2str(h,'%.3f');
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

surv_type = '1'; %0 = not valid, 1 = GPS, 2 = DGPS
if (nargin < 4)
    HDOP = 1;      %fake HDOP value
end
h_unit = 'M';

if (nargin > 1)
    nsat = num2str(nsat);
    if (length(nsat) == 1)
        [nsat] = sprintf('0%s',nsat);
    end
else
    nsat = '';
end

%-----------------------------------------------------------------------------------------------
% FORMAT DATA
%-----------------------------------------------------------------------------------------------

hour = num2str(date(1,4));
minute = num2str(date(1,5));
second = num2str(floor(date(1,6)));

[null, nchar] = size(hour); %#ok<ASGLU>
if (nchar == 1)
    [hour] = sprintf('0%s',hour);
end

[null, nchar] = size(minute); %#ok<ASGLU>
if (nchar == 1)
    [minute] = sprintf('0%s',minute);
end

[null, nchar] = size(second); %#ok<ASGLU>
if (nchar == 1)
    [second] = sprintf('0%s',second);
end

decsec = '.00';

%-----------------------------------------------------------------------------------------------
% COMPOSITION OF THE NMEA SENTENCE
%-----------------------------------------------------------------------------------------------

nmeastring = sprintf('$GPGGA,%s%s%s%s,%s,%c,%s,%c,%s,%s,%.2f,%s,%c,,,,',hour,minute,second,decsec,phi_nmea,emi_NS,lam_nmea,emi_EW,surv_type,nsat,HDOP,h,h_unit);

%-----------------------------------------------------------------------------------------------
% CHECKSUM COMPUTATION
%-----------------------------------------------------------------------------------------------

checksum = NMEA_checksum(nmeastring);
nmeastring = [nmeastring '*' checksum];
