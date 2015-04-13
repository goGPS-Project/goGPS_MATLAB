function nmeastring = NMEA_GGA_gen(pos_R, nsat, time, HDOP, mode)

% SYNTAX:
%   nmeastring = NMEA_GGA_gen(pos_R, nsat, time, HDOP, mode);
%
% INPUT:
%   pos_R = estimated ROVER position (X,Y,Z)
%   nsat = number of visible satellites
%   time = GPS time of measurements
%   HDOP = horizontal dilution of precision
%   mode = positioning mode
%
% OUTPUT:
%   nmeastring = $GPGGA sentence (NMEA)
%
% DESCRIPTION:
%   Returns a $GPGGA sentence in NMEA 0183 format.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
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

global geoid

%date
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

%initialize geoid ondulation
N = [];
N_unit = [];

if (exist('geoid','var') && isfield(geoid,'ncols') && geoid.ncols ~= 0)
    %geoid ondulation interpolation
    N = grid_bilin_interp(lam, phi, geoid.grid, geoid.ncols, geoid.nrows, geoid.cellsize, geoid.Xll, geoid.Yll, -9999);
    %orthometric height
    h = h - N;
    %ondulation units
    N_unit = 'M';
end

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

%height
% [h_str] = sprintf('%d', floor(h));
% if (length(h_str) == 1)
%     [h] = sprintf('00%.1f', h);
% elseif (length(h_str) == 2)
%     [h] = sprintf('0%.1f', h);
% else
%     [h] = sprintf('%.1f', h);
% end

%N
% [N_str] = sprintf('%d', floor(N));
% if (length(N_str) == 1)
%     [N] = sprintf('00%.1f', N);
% elseif (length(N_str) == 2)
%     [N] = sprintf('0%.1f', N);
% else
%     [N] = sprintf('%.1f', N);
% end

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

% survey type
% 0 = not valid
% 1 = code undifferenced
% 2 = code double differences
% 4 = code and phase double diff. (integer ambiguities fixed)
% 5 = code and phase double diff. (float ambiguities)
if (nargin > 4)
    if (mode == 1 | mode == 2)
        surv_type = '1';
    elseif (mode == 11 | mode == 12)
        surv_type = '2';
    elseif (mode == 13)
        surv_type = '2';
    elseif (mode == 14)
        surv_type = '5';
    else
        surv_type = '1';
    end
else
    surv_type = '1';
end

if (nargin < 4)
    HDOP = 1;      %fake HDOP value
end
h_unit = 'M';

if (nargin > 1)
    [nsat] = sprintf('%d', nsat);
    if (length(nsat) == 1)
        [nsat] = sprintf('0%s',nsat);
    end
else
    nsat = '';
end

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

% decsec = '.00';

%-----------------------------------------------------------------------------------------------
% COMPOSITION OF THE NMEA SENTENCE
%-----------------------------------------------------------------------------------------------

nmeastring = sprintf('$GPGGA,%s%s%s,%s,%c,%s,%c,%s,%s,%.2f,%.1f,%c,%.1f,%c,,',hour,minute,second,phi_nmea,emi_NS,lam_nmea,emi_EW,surv_type,nsat,HDOP,h,h_unit,N,N_unit);

%-----------------------------------------------------------------------------------------------
% CHECKSUM COMPUTATION
%-----------------------------------------------------------------------------------------------

checksum = NMEA_checksum(nmeastring);
nmeastring = [nmeastring '*' checksum];
