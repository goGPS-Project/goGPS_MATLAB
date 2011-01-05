function [Eph, iono] = RINEX_get_nav(file_nav)

% SYNTAX:
%   [Eph, iono] = RINEX_get_nav(file_nav);
%
% INPUT:
%   file_nav = RINEX navigation file
%
% OUTPUT:
%   Eph = matrix containing 29 ephemerides for each satellite
%   iono = matrix containing ionosphere parameters
%
% DESCRIPTION:
%   Parse a RINEX navigation file.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.3 alpha
%
% Copyright (C) 2009-2010 Mirko Reguzzoni, Eugenio Realini
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
iono = zeros(8,1);

%open navigation file
fid = fopen(file_nav,'rt');

%read the header
header_end = [];
while (isempty(header_end))
    %read the line and search the 'ION ALPHA' label
    lin = fgetl(fid);
    iono_found = findstr(lin,'ION ALPHA');

    %if the label was found
    if ~isempty(iono_found)
        %change flag
        ioparam = 1;
        %save the 8 ionosphere parameters
        data = textscan(lin,'%f%f%f%f%*[^\n]');
        iono(1) = data{1};
        iono(2) = data{2};
        iono(3) = data{3};
        iono(4) = data{4};
        lin = [];
        while isempty(lin)
            lin = fgetl(fid);
        end
        data = textscan(lin,'%f%f%f%f%*[^\n]');
        iono(5) = data{1};
        iono(6) = data{2};
        iono(7) = data{3};
        iono(8) = data{4};
    end

    header_end = findstr(lin,'END OF HEADER');
end

%if ionosphere parameters were not found
if (ioparam == 0)
    fprintf('Warning: ionosphere parameters not found in navigation file\n');
end

i = 0;

%parse the rest of the file and store ephemerides
while (~feof(fid))
    
    lin1 = [];
    lin2 = [];
    lin3 = [];
    lin4 = [];
    lin5 = [];
    lin6 = [];
    lin7 = [];
    lin8 = [];

    i = i+1;
    while isempty(lin1)
        lin1 = fgetl(fid);
    end
    while isempty(lin2)
        lin2 = fgetl(fid);
    end
    while isempty(lin3)
        lin3 = fgetl(fid);
    end
    while isempty(lin4)
        lin4 = fgetl(fid);
    end
    while isempty(lin5)
        lin5 = fgetl(fid);
    end
    while isempty(lin6)
        lin6 = fgetl(fid);
    end
    while isempty(lin7)
        lin7 = fgetl(fid);
    end
    while isempty(lin8)
        lin8 = fgetl(fid);
    end
    
    if (lin1 == -1)
        break
    end

    svprn  = str2num(lin1(1:2)); %#ok<*ST2NM>
    year   = str2num(lin1(3:6)); %#ok<NASGU>
    month  = str2num(lin1(7:9)); %#ok<NASGU>
    day    = str2num(lin1(10:12)); %#ok<NASGU>
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

    i0       = str2num(lin5(4:22));
    crc      = str2num(lin5(23:41));
    omega    = str2num(lin5(42:60));
    Omegadot = str2num(lin5(61:79));

    idot       = str2num(lin6(4:22));
    code_on_L2 = str2num(lin6(23:41));
    weekno     = str2num(lin6(42:60));
    L2flag     = str2num(lin6(61:79));

    svaccur  = str2num(lin7(4:22));
    svhealth = str2num(lin7(23:41));
    tgd      = str2num(lin7(42:60));
    iodc     = str2num(lin7(61:79));

    tom     = str2num(lin8(4:22));
    if (length(lin8) > 22)
        fit_int = str2num(lin8(23:41));
    else
        fit_int = 0;
    end

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
    Eph(22,i) = IODE;
    Eph(23,i) = code_on_L2;
    Eph(24,i) = weekno;
    Eph(25,i) = L2flag;
    Eph(26,i) = svaccur;
    Eph(27,i) = svhealth;
    Eph(28,i) = tgd;
    Eph(29,i) = fit_int;
    
    %if IODC and IODE do not match, issue a warning
    if (iodc ~= IODE)
        fprintf('Warning: IODE and IODC values do not match (ephemerides for satellite %02d, time %dh %dm %.1fs)\n',svprn,hour,minute,second);
    end
end
