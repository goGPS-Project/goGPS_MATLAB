function GSVlines = NMEA_GSV_gen(sat, el, az, snr)

% SYNTAX:
%   GSVlines = NMEA_GSV_gen(sat, el, az, snr);
%
% INPUT:
%   sat = list of visible satellites
%   el = elevation [deg]
%   az = azimuth [deg]
%   snr = signal-to-noise ratio [dB]
%
% OUTPUT:
%   GSVlines = $GPGSV sentence(s) (NMEA)
%
% DESCRIPTION:
%   Returns (a) $GPGSV sentence(s) in NMEA 0183 format.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.3 alpha
%
% Copyright (C) 2009-2011 Mirko Reguzzoni, Eugenio Realini
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

%number of satellites
nsat = size(sat,1);

%-----------------------------------------------------------------------------------------------
% COMPOSITION OF THE NMEA SENTENCE
%-----------------------------------------------------------------------------------------------

%number of GSV sentences (max 4 satellites per sentence)
n = ceil(nsat/4);

%variable to store the n GSV sentences
GSVlines = [];

for i = 1 : n
    nmeastring = sprintf('$GPGSV,%d,%d,%d', n, i, nsat);
    for j = 1 : 4
        index = 4*(i-1) + j;
        if (index <= nsat)
            nmeastring = [nmeastring sprintf(',%d,%d,%d,%d', sat(index), round(el(index)), round(az(index)), round(snr(index)))];
        else
            nmeastring = [nmeastring ',,,,'];
        end
    end

    %checksum computation
    checksum = NMEA_checksum(nmeastring);
    nmeastring = [nmeastring '*' checksum];
    
    %add the new string to the GSV group
    GSVlines = [GSVlines nmeastring];
    
    %add new line character between GSV sentences
    %(but not at the end of the overall string)   
    if (i ~= n)
        GSVlines = [GSVlines '\n'];
    end    
end
