function [Eph] = RINEX_get_nav_GLO(file_nav_GLO)

% SYNTAX:
%   [Eph] = RINEX_get_nav_GLO(file_nav_GLO);
%
% INPUT:
%   file_nav_GLO = RINEX navigation file (GLONASS)
%
% OUTPUT:
%   Eph = matrix containing 17 ephemerides for each satellite
%
% DESCRIPTION:
%   Parse a RINEX GLONASS navigation file.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.0 beta
%
% Copyright (C) 2009-2012 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Sara Lucca
%
% Partially based on RINEXE.M (EASY suite) by Kai Borre
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

Eph = [];

%open navigation file
fid = fopen(file_nav_GLO,'rt');

answer = [];

%search for the end of the header
while (isempty(answer))
    %lettura della riga e cerco la scritta
    lin = fgetl(fid);
    answer = findstr(lin,'END OF HEADER');
end

i=0;

%parse the rest of the file and store ephemerides
while (~feof(fid))

    i = i+1;
    lin1 = fgetl(fid);
    lin2 = fgetl(fid);
    lin3 = fgetl(fid); %#ok<NASGU>
    lin4 = fgetl(fid); %#ok<NASGU>

    svprn  = str2num(lin1(1:2));
    year   = str2num(lin1(3:6)); %#ok<NASGU>
    month  = str2num(lin1(7:9)); %#ok<NASGU>
    day    = str2num(lin1(10:12)); %#ok<NASGU>
    hour   = str2num(lin1(13:15)); %#ok<NASGU>
    minute = str2num(lin1(16:18)); %#ok<NASGU>
    second = str2num(lin1(19:22)); %#ok<NASGU>
    TauN   = -str2num(lin1(23:41));
    GammaN = str2num(lin1(42:60));
    tk     = str2num(lin1(61:79));

    X      = str2num(lin2(4:22));
    vX     = str2num(lin2(23:41));
    aX     = str2num(lin2(42:60));
    Bn     = str2num(lin2(61:79));

    Y      = str2num(lin2(4:22));
    vY     = str2num(lin2(23:41));
    aY     = str2num(lin2(42:60));
    freq_num = str2num(lin2(61:79));

    %frequencies on L1 and L2
    freq_L1 = freq_num * 0.5625 + 1602.0;
    freq_L2 = freq_num * 0.4375 + 1246.0;

    Z      = str2num(lin2(4:22));
    vZ     = str2num(lin2(23:41));
    aZ     = str2num(lin2(42:60));
    E      = str2num(lin2(61:79));

    %save ephemerides
    Eph(1,i)  = svprn;
    Eph(2,i)  = TauN;
    Eph(3,i)  = GammaN;
    Eph(4,i)  = tk;
    Eph(5,i)  = X;
    Eph(6,i)  = vX;
    Eph(7,i)  = aX;
    Eph(8,i)  = Bn;
    Eph(9,i)  = Y;
    Eph(10,i) = vY;
    Eph(11,i) = aY;
    Eph(12,i) = freq_L1;
    Eph(13,i) = freq_L2;
    Eph(14,i) = Z;
    Eph(15,i) = vZ;
    Eph(16,i) = aZ;
    Eph(17,i) = E;
end
