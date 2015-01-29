function [Eph, iono] = RINEX_get_nav(file_nav, constellations)

% SYNTAX:
%   [Eph, iono] = RINEX_get_nav(file_nav, constellations);
%
% INPUT:
%   file_nav = RINEX navigation file
%   constellations = struct with multi-constellation settings (see goGNSS.initConstellation)
%
% OUTPUT:
%   Eph = matrix containing 33 navigation parameters for each satellite
%   iono = matrix containing ionosphere parameters
%
% DESCRIPTION:
%   Parse a RINEX navigation file.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
% Portions of code contributed by Damiano Triglione (2012)
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

% ioparam = 0;
Eph = [];
iono = zeros(8,1);

if (nargin < 2 || isempty(constellations)) %then use only GPS as default
    constellations.GPS = struct('numSat', 32, 'enabled', 1, 'indexes', [1:32], 'PRN', [1:32]);
    constellations.nEnabledSat = 32;
    constellations.indexes = constellations.GPS.indexes;
    constellations.PRN     = constellations.GPS.PRN;
end

fprintf('%s',['Reading RINEX file ' file_nav ': ... ']);

%open navigation file
fid = fopen(file_nav,'rt');

%read the header
header_end = [];
while (isempty(header_end))
    %read the line and search the ionosphere labels
    lin = fgetl(fid);
    
    vers_found =  ~isempty(strfind(lin,'RINEX VERSION / TYPE'));
    iono_found = (~isempty(strfind(lin,'ION ALPHA')) || ~isempty(strfind(lin,'IONOSPHERIC CORR')));

    %if the ionosphere parameters label was found
    if (vers_found)
        version = str2num(lin(1:9));
    end
    
    %if the ionosphere parameters label was found
    if (iono_found)
        %change flag
        %         ioparam = 1;
        %save the 8 ionosphere parameters
        data = textscan(lin(5:end),'%f%f%f%f%*[^\n]');
        if ~isempty(data(4))
            iono(1) = data{1};
            iono(2) = data{2};
            iono(3) = data{3};
            iono(4) = data{4};
            lin = [];
            while isempty(lin)
                lin = fgetl(fid);
            end
            data = textscan(lin(5:end),'%f%f%f%f%*[^\n]');
            if ~isempty(data(4))
                iono(5) = data{1};
                iono(6) = data{2};
                iono(7) = data{3};
                iono(8) = data{4};
            else
                iono = zeros(8,1);
            end
        end
    end
    
    header_end = strfind(lin,'END OF HEADER');
end

% %if ionosphere parameters were not found
% if (ioparam == 0)
%     fprintf('... WARNING: ionosphere parameters not found in navigation file\n');
% end

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

    %read the first line (containing system and time information)
    while isempty(lin1)
        lin1 = fgetl(fid);
    end
    if (lin1 == -1)
        break
    end
    if (~isempty(strfind(lin1,'COMMENT')))
        continue
    end

    %character offset (to deal with various RINEX versions)
    sys_id   = lin1(1);
    sys_uint = uint8(sys_id);
    if (strcmp(sys_id,'G') || strcmp(sys_id,'R') || strcmp(sys_id,'E') || strcmp(sys_id,'C') || strcmp(sys_id,'J') || strcmp(sys_id,'S'))
        o = 1;                 %RINEX v2.12(not GPS) or v3.xx
        if (strcmp(sys_id,'G'))
            sys_index = constellations.GPS.indexes(1);
        elseif (strcmp(sys_id,'R'))
            sys_index = constellations.GLONASS.indexes(1);
        elseif (strcmp(sys_id,'E'))
            sys_index = constellations.Galileo.indexes(1);
        elseif (strcmp(sys_id,'C'))
            sys_index = constellations.BeiDou.indexes(1);
        elseif (strcmp(sys_id,'J'))
            sys_index = constellations.QZSS.indexes(1);
        elseif (strcmp(sys_id,'S'))
            %sys_index = constellations.SBAS.indexes(1);
        end
    elseif (sys_uint == 32 || (sys_uint >= 48 && sys_uint <= 57)) %if blank space or number
        if (strcmpi(file_nav(end),'g'))
            sys_id = 'R';
            o = 0;                 %RINEX v<=2.12, GLONASS
            sys_index = constellations.GLONASS.indexes(1);
        else
            sys_id = 'G';
            o = 0;                 %RINEX v<=2.12, GPS
            sys_index = constellations.GPS.indexes(1);
        end
    else
        return
    end
    
    %read the next 3 lines (common to all systems)
    while isempty(lin2)
        lin2 = fgetl(fid);
    end
    while isempty(lin3)
        lin3 = fgetl(fid);
    end
    while isempty(lin4)
        lin4 = fgetl(fid);
    end
    
    %if not a GLONASS or SBAS entry, read also the next 4 lines
    if (~strcmp(sys_id, 'R') && ~strcmp(sys_id, 'S'))
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
            %while lin8(end) == ' ', lin8 = lin8(1:end-1); end
            lin8 = RemoveUnwantedTrailingSpaces(lin8);
        end
    end
    
    switch sys_id
        case 'G'
            if (~constellations.GPS.enabled), continue, end
        case 'R'
            if (~constellations.GLONASS.enabled), continue, end
        case 'E'
            if (~constellations.Galileo.enabled), continue, end
        case 'C'
            if (~constellations.BeiDou.enabled), continue, end
        case 'J'
            if (~constellations.QZSS.enabled), continue, end
        case 'S'
            if (~constellations.SBAS.enabled), continue, end
    end

    svprn  = str2num(lin1(o+[1:2])); %When input is a scalar, str2double is better than str2num. But str2double does not support 'D'
    if (version < 3)
        year   = str2num(lin1(o+[3:6])); year = four_digit_year(year);
        month  = str2num(lin1(o+[7:9]));
        day    = str2num(lin1(o+[10:12]));
        hour   = str2num(lin1(o+[13:15]));
        minute = str2num(lin1(o+[16:18]));
        second = str2num(lin1(o+[19:22]));
    else
        year   = str2num(lin1(o+[3:6]+1));
        month  = str2num(lin1(o+[7:9]+1));
        day    = str2num(lin1(o+[10:12]+1));
        hour   = str2num(lin1(o+[13:15]+1));
        minute = str2num(lin1(o+[16:18]+1));
        second = str2num(lin1(o+[19:21]+1));
    end
    
    i = i+1;
    
    %if not GLONASS
    if (~strcmp(sys_id, 'R'))

        af0    = str2num(lin1(o+[23:41]));
        af1    = str2num(lin1(o+[42:60]));
        af2    = str2num(lin1(o+[61:79]));

        IODE   = str2num(lin2(o+[4:22]));
        crs    = str2num(lin2(o+[23:41]));
        deltan = str2num(lin2(o+[42:60]));
        M0     = str2num(lin2(o+[61:79]));
        
        cuc    = str2num(lin3(o+[4:22]));
        ecc    = str2num(lin3(o+[23:41]));
        cus    = str2num(lin3(o+[42:60]));
        roota  = str2num(lin3(o+[61:79]));
        
        toe    = str2num(lin4(o+[4:22]));
        cic    = str2num(lin4(o+[23:41]));
        Omega0 = str2num(lin4(o+[42:60]));
        cis    = str2num(lin4(o+[61:79]));
        
        i0       = str2num(lin5(o+[4:22]));
        crc      = str2num(lin5(o+[23:41]));
        omega    = str2num(lin5(o+[42:60]));
        Omegadot = str2num(lin5(o+[61:79]));
        
        idot       = str2num(lin6(o+[4:22]));
        code_on_L2 = str2num(lin6(o+[23:41]));
        weekno     = str2num(lin6(o+[42:60]));
        L2flag     = str2num(lin6(o+[61:79]));
        if (isempty(L2flag))
            L2flag = 0;
        end

        svaccur  = str2num(lin7(o+[4:22]));
        svhealth = str2num(lin7(o+[23:41]));
        tgd      = str2num(lin7(o+[42:60]));
        iodc     = str2num(lin7(o+[61:79]));

        tom      = str2num(lin8(o+[4:22])); %#ok<NASGU>
        if (length(lin8) > o+22)
            fit_int = str2num(lin8(o+[23:41]));
        else
            fit_int = 0;
        end
        
        [~, toc] = date2gps([year month day hour minute second]);
        
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
        Eph(21,i) = toc;
        Eph(22,i) = IODE;
        Eph(23,i) = code_on_L2;
        Eph(24,i) = weekno; if (strcmp(sys_id, 'E') && weekno > 2500); weekno = weekno - 1024; end;
        Eph(25,i) = L2flag;
        Eph(26,i) = svaccur;
        Eph(27,i) = svhealth;
        Eph(28,i) = tgd;
        Eph(29,i) = fit_int;
        Eph(30,i) = (sys_index-1) + svprn; %satellite index (consistent with other observation arrays)
        Eph(31,i) = int8(sys_id);
        Eph(32,i) = weektow2time(weekno, toe, sys_id);
        Eph(33,i) = weektow2time(weekno, toc, sys_id);
        
        %if IODC and IODE do not match, issue a warning
        if (iodc ~= IODE && ~strcmp(sys_id, 'C') && ~strcmp(sys_id, 'E'))
            fprintf('... WARNING: IODE and IODC values do not match (ephemerides for satellite %1s%02d, time %dh %dm %.1fs)\n',sys_id,svprn,hour,minute,second);
        end
        
    %if GLONASS
    else
        TauN   = -str2num(lin1(o+[23:41]));
        GammaN = str2num(lin1(o+[42:60]));
        tk     = str2num(lin1(o+[61:79]));
        
        X      = str2num(lin2(o+[4:22]));
        Xv     = str2num(lin2(o+[23:41]));
        Xa     = str2num(lin2(o+[42:60]));
        Bn     = str2num(lin2(o+[61:79])); %health flag
        
        Y      = str2num(lin3(o+[4:22]));
        Yv     = str2num(lin3(o+[23:41]));
        Ya     = str2num(lin3(o+[42:60]));
        freq_num = str2num(lin3(o+[61:79])); %frequency number
        
        Z      = str2num(lin4(o+[4:22]));
        Zv     = str2num(lin4(o+[23:41]));
        Za     = str2num(lin4(o+[42:60]));
        E      = str2num(lin4(o+[61:79])); %age of oper. information  (days)
        
        %frequencies on L1 and L2
        %freq_L1 = freq_num * 0.5625 + 1602.0;
        %freq_L2 = freq_num * 0.4375 + 1246.0;
        
        %convert GLONASS (UTC) date to GPS date
        date_GLO = datenum([year month day hour minute second]);
        date_GPS = utc2gps(date_GLO);
        
        %convert GPS date to seconds of week (used as GLONASS time-of-ephemeris)
        [week_toe, toe] = date2gps(datevec(date_GPS));
        
        %save ephemerides (position, velocity and acceleration vectors in ECEF system PZ-90.02)
        Eph(1,i)  = svprn;
        Eph(2,i)  = TauN;
        Eph(3,i)  = GammaN;
        Eph(4,i)  = tk;
        Eph(5,i)  = X*1e3;  %satellite X coordinate at ephemeris reference time [m]
        Eph(6,i)  = Y*1e3;  %satellite Y coordinate at ephemeris reference time [m]
        Eph(7,i)  = Z*1e3;  %satellite Z coordinate at ephemeris reference time [m]
        Eph(8,i)  = Xv*1e3; %satellite velocity along X at ephemeris reference time [m/s]
        Eph(9,i)  = Yv*1e3; %satellite velocity along Y at ephemeris reference time [m/s]
        Eph(10,i) = Zv*1e3; %satellite velocity along Z at ephemeris reference time [m/s]
        Eph(11,i) = Xa*1e3; %acceleration due to lunar-solar gravitational perturbation along X at ephemeris reference time [m/s^2]
        Eph(12,i) = Ya*1e3; %acceleration due to lunar-solar gravitational perturbation along Y at ephemeris reference time [m/s^2]
        Eph(13,i) = Za*1e3; %acceleration due to lunar-solar gravitational perturbation along Z at ephemeris reference time [m/s^2]
        Eph(14,i) = E;
        Eph(15,i) = freq_num;
        Eph(16,i) = 0;
        Eph(17,i) = 0;
        Eph(18,i) = toe;
        Eph(19,i) = 0;
        Eph(20,i) = 0;
        Eph(21,i) = 0;
        Eph(22,i) = 0;
        Eph(23,i) = 0;
        Eph(24,i) = week_toe;
        Eph(25,i) = 0;
        Eph(26,i) = 0;
        Eph(27,i) = Bn; %health flag
        Eph(28,i) = 0;
        Eph(29,i) = 0;
        Eph(30,i) = (sys_index-1) + svprn;
        Eph(31,i) = int8(sys_id);
        Eph(32,i) = weektow2time(week_toe, toe, sys_id);
        Eph(33,i) = 0;
    end
end

fclose(fid);

fprintf(['done\n']);
