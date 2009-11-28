function [Eph, iono] = RINEX_get_nav(file_nav)

% SYNTAX:
%   [Eph, iono] = RINEX_get_nav(file_nav);
%
% INPUT:
%   file_nav = RINEX navigation file
%
% OUTPUT:
%   Eph = matrix containing 21 ephemerides for each satellite
%   iono = matrix containing ionosphere parameters
%
% DESCRIPTION:
%   Parse a RINEX navigation file.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Media Center, Osaka City University, Japan
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

ioparam = 0;
Eph = [];
ion = zeros(8,1);

%open navigation file
fid = fopen(file_nav,'rt');

%find ionosphere parameters
while ((ioparam==0)&(~feof(fid)))

    %read the line and search the 'ION ALPHA' label
    lin = fgetl(fid);
    answer = findstr(lin,'ION ALPHA');

    %if the label was found
    if ~isempty(answer)
        %change flag
        ioparam = 1;
        %save the 8 ionosphere parameters
        [dato,lin] = strtok(lin);
        ion(1) = str2num(dato);
        [dato,lin] = strtok(lin);
        ion(2) = str2num(dato);
        [dato,lin] = strtok(lin);
        ion(3) = str2num(dato);
        [dato,lin] = strtok(lin);
        ion(4) = str2num(dato);
        lin = fgetl(fid);
        [dato,lin] = strtok(lin);
        ion(5) = str2num(dato);
        [dato,lin] = strtok(lin);
        ion(6) = str2num(dato);
        [dato,lin] = strtok(lin);
        ion(7) = str2num(dato);
        [dato,lin] = strtok(lin);
        ion(8) = str2num(dato);
        iono = ion;
    end

end

%if ionosphere parameters were not found
if (ioparam == 0)
    error('Warning: ionosphere parameters not found in navigation file %s',file_nav);
end

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
    lin3 = fgetl(fid);
    lin4 = fgetl(fid);
    lin5 = fgetl(fid);
    lin6 = fgetl(fid);
    lin7 = fgetl(fid);
    lin8 = fgetl(fid);

    svprn  = str2num(lin1(1:2));
    year   = str2num(lin1(3:6));
    month  = str2num(lin1(7:9));
    day    = str2num(lin1(10:12));
    hour   = str2num(lin1(13:15));
    minute = str2num(lin1(16:18));
    second = str2num(lin1(19:22));
    af0    = str2num(lin1(23:41));
    af1    = str2num(lin1(42:60));
    af2    = str2num(lin1(61:79));

    IODE   = str2num(lin2(4:22));
    crs    = str2num(lin2(23:41));
    deltan = str2num(lin2(42:60));
    M0     = str2num(lin2(61:79));

    cuc    = str2num(lin3(4:22));
    ecc    = str2num(lin3(23:41));
    cus    = str2num(lin3(42:60));
    roota  = str2num(lin3(61:79));

    toe    = str2num(lin4(4:22));
    cic    = str2num(lin4(23:41));
    Omega0 = str2num(lin4(42:60));
    cis    = str2num(lin4(61:79));

    i0     = str2num(lin5(4:22));
    crc    = str2num(lin5(23:41));
    omega  = str2num(lin5(42:60));
    Omegadot = str2num(lin5(61:79));

    idot   = str2num(lin6(4:22));
    codes  = str2num(lin6(23:41));
    weekno = str2num(lin6(42:60));
    L2flag = str2num(lin6(61:79));

    svaccur  = str2num(lin7(4:22));
    svhealth = str2num(lin7(23:41));
    tgd    = str2num(lin7(42:60));
    iodc   = lin7(61:79);

    tom = str2num(lin8(4:22));

    %debugging    
    %jd = julday(year+2000, month, day, 0);
    %[week, sec_of_week] = gps_time(jd);
    %time = sec_of_week + hour*3600+minute*60+second;
    %time-toe

    %save ephemerides
    Eph(1,i)  = svprn;
    Eph(2,i)  = af2;
    Eph(3,i)  = M0;
    Eph(4,i)  = roota;
    Eph(5,i)  = deltan;
    Eph(6,i)  = ecc;
    Eph(7,i)  = omega;
    Eph(8,i)  = cuc;
    Eph(9,i)  = cus;
    Eph(10,i) = crc;
    Eph(11,i) = crs;
    Eph(12,i) = i0;
    Eph(13,i) = idot;
    Eph(14,i) = cic;
    Eph(15,i) = cis;
    Eph(16,i) = Omega0;
    Eph(17,i) = Omegadot;
    Eph(18,i) = toe;
    Eph(19,i) = af0;
    Eph(20,i) = af1;
    Eph(21,i) = tom;
end
