function nmeastring = NMEA_string_generator(pos_R,nsat)

% SYNTAX:
%   nmeastring = NMEA_string_generator(pos_R,nsat);
%
% INPUT:
%   pos_R = estimated ROVER position (X,Y,Z)
%   nsat = number of visible satellites
%
% OUTPUT:
%   nmeastring = $GPGGA sentence (NMEA)
%
% DESCRIPTION:
%   Returns a $GPGGA sentence in NMEA 0183 format.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Media Center, Osaka City University, Japan
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

%get current date
current_date = clock;

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
HDOP = '1';      %fake HDOP value: to be changed
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

hour = num2str(current_date(1,4));
minute = num2str(current_date(1,5));
second = num2str(floor(current_date(1,6)));

[null, ncifre] = size(hour);
if (ncifre == 1)
    [hour] = sprintf('0%s',hour);
end

[null, ncifre] = size(minute);
if (ncifre == 1)
    [minute] = sprintf('0%s',minute);
end

[null, ncifre] = size(second);
if (ncifre == 1)
    [second] = sprintf('0%s',second);
end

decsec = '.00';

%-----------------------------------------------------------------------------------------------
% COMPOSITION OF THE NMEA SENTENCE
%-----------------------------------------------------------------------------------------------

nmeastring = sprintf('$GPGGA,%s%s%s%s,%s,%c,%s,%c,%s,%s,%s,%s,%c,,,,',hour,minute,second,decsec,phi_nmea,emi_NS,lam_nmea,emi_EW,surv_type,nsat,HDOP,h,h_unit);

%-----------------------------------------------------------------------------------------------
% CHECKSUM COMPUTATION (found at http://www.mathworks.com/matlabcentral/fileexchange/15080)
%-----------------------------------------------------------------------------------------------

checksum = 0;

% see if string contains the * which starts the checksum and keep string
% upto * for generating checksum
nmeastring = strtok(nmeastring,'*');

nmeastring_d = double(nmeastring);                    % convert characters in string to double values
for count = 2:length(nmeastring)                      % checksum computation ignores $ at start
    checksum = bitxor(checksum,nmeastring_d(count));  % checksum computation
    checksum = uint16(checksum);                      % make sure that checksum is unsigned int16
end

% convert checksum to hex value
checksum = double(checksum);
checksum = dec2hex(checksum);

% add leading zero to checksum if it is a single digit, e.g. 4 has a 0
% added so that the checksum is 04
if length(checksum) == 1
    checksum = strcat('0',checksum);
end

nmeastring = [nmeastring '*' checksum];
