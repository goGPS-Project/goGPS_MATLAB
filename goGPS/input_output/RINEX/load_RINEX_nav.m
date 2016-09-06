function [Eph, iono, flag_return] = load_RINEX_nav(filename, constellations, flag_SP3, iono_model, time, wait_dlg)

% SYNTAX:
%   [Eph, iono, flag_return] = load_RINEX_nav(filename, constellations, flag_SP3, iono_model, time, wait_dlg);
%
% INPUT:
%   filename = RINEX navigation file
%   constellations = struct with multi-constellation settings
%                   (see 'multi_constellation_settings.m' - empty if not available)
%   flag_SP3 = boolean flag to indicate SP3 availability
%   wait_dlg = optional handler to waitbar figure (optional)
%
% OUTPUT:
%   Eph = matrix containing 33 navigation parameters for each satellite
%   iono = vector containing ionosphere parameters
%   flag_return = notify the parent function that it should return 
%                 (downloaded navigation file still compressed).
%
% DESCRIPTION:
%   Parses RINEX navigation files.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2013 Mirko Reguzzoni,Eugenio Realini
% Portions of code contributed by Damiano Triglione (2012)
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

% Check the input arguments
if (nargin < 6)
    wait_dlg_PresenceFlag = false;
else
    wait_dlg_PresenceFlag = true;
end

flag_return = 0;

if (isempty(constellations)) %then use only GPS as default
    [constellations] = multi_constellation_settings(1, 0, 0, 0, 0, 0);
end
tic;

%number of satellite slots for enabled constellations
nSatTot = constellations.nEnabledSat;

%read navigation files
if (~flag_SP3)
    parse_file(0);
else
    Eph = zeros(33,nSatTot);
    iono = zeros(8,1);
end

%if Klobuchar ionospheric delay correction is requested but parameters are not available in the navigation file, try to download them
if (iono_model == 2 & ~any(iono))
    [week, sow] = time2weektow(time(1));
    [date, DOY] = gps2date(week, sow);
    
    filename_brdc = ['brdc' num2str(DOY,'%03d') '0.' num2str(two_digit_year(date(1,1)),'%02d') 'n'];
    filename_CGIM = ['CGIM' num2str(DOY,'%03d') '0.' num2str(two_digit_year(date(1,1)),'%02d') 'N'];
    
    pos = find(filename == '/'); if(isempty(pos)), pos = find(filename == '\'); end;
    nav_path = filename(1:pos(end));
    
    file_avail = 0;
    if (exist([nav_path filename_CGIM],'file'))
        filename = [nav_path filename_CGIM];
        file_avail = 1;
    elseif (exist([nav_path filename_brdc],'file'))
        filename = [nav_path filename_brdc];
        file_avail = 1;
    else
        filename = filename_brdc;
        [download_successful, compressed] = download_nav(filename);
        filename = [nav_path filename];
        if (download_successful)
            file_avail = 1;
        end
        if (compressed)
            flag_return = 1;
        end
    end
    
    if (file_avail)
        parse_file(1);
    end
end

    function parse_file(only_iono)
        
        if (wait_dlg_PresenceFlag)
            waitbar(0.5,wait_dlg,'Reading navigation files...')
        end
        
        Eph_G = []; iono_G = zeros(8,1);
        Eph_R = []; iono_R = zeros(8,1);
        Eph_E = []; iono_E = zeros(8,1);
        Eph_C = []; iono_C = zeros(8,1);
        Eph_J = []; iono_J = zeros(8,1);
        
        if (strcmpi(filename(end),'p'))
            flag_mixed = 1;
        else
            flag_mixed = 0;
        end
        
        if (constellations.GPS.enabled || flag_mixed || only_iono)
            if (exist(filename,'file'))
                %parse RINEX navigation file (GPS) NOTE: filename expected to
                %end with 'n' or 'N' (GPS) or with 'p' or 'P' (mixed GNSS)
                if(~only_iono), fprintf('%s',['Reading RINEX file ' filename ': ... ']); end
                [Eph_G, iono_G] = RINEX_get_nav(filename, constellations);
                if(~only_iono), fprintf('done\n'); end
            else
                fprintf('... WARNING: GPS navigation file not found. Disabling GPS positioning. \n');
                constellations.GPS.enabled = 0;
            end
        end
        
        if (constellations.GLONASS.enabled && ~only_iono)
            if (exist([filename(1:end-1) 'g'],'file'))
                %parse RINEX navigation file (GLONASS)
                if(~only_iono), fprintf('%s',['Reading RINEX file ' filename ': ... ']); end
                [Eph_R, iono_R] = RINEX_get_nav([filename(1:end-1) 'g'], constellations);
                if(~only_iono), fprintf('done\n'); end
            elseif (~flag_mixed)
                fprintf('... WARNING: GLONASS navigation file not found. Disabling GLONASS positioning. \n');
                constellations.GLONASS.enabled = 0;
            end
        end
        
        if (constellations.Galileo.enabled && ~only_iono)
            if (exist([filename(1:end-1) 'l'],'file'))
                %parse RINEX navigation file (Galileo)
                if(~only_iono), fprintf('%s',['Reading RINEX file ' filename ': ... ']); end
                [Eph_E, iono_E] = RINEX_get_nav([filename(1:end-1) 'l'], constellations);
                if(~only_iono), fprintf('done\n'); end
            elseif (~flag_mixed)
                fprintf('... WARNING: Galileo navigation file not found. Disabling Galileo positioning. \n');
                constellations.Galileo.enabled = 0;
            end
        end
        
        if (constellations.BeiDou.enabled && ~only_iono)
            if (exist([filename(1:end-1) 'c'],'file'))
                %parse RINEX navigation file (BeiDou)
                if(~only_iono), fprintf('%s',['Reading RINEX file ' filename ': ... ']); end
                [Eph_C, iono_C] = RINEX_get_nav([filename(1:end-1) 'c'], constellations);
                if(~only_iono), fprintf('done\n'); end
            elseif (~flag_mixed)
                fprintf('... WARNING: BeiDou navigation file not found. Disabling BeiDou positioning. \n');
                constellations.BeiDou.enabled = 0;
            end
        end
        
        if (constellations.QZSS.enabled && ~only_iono)
            if (exist([filename(1:end-1) 'q'],'file'))
                %parse RINEX navigation file (QZSS)
                if(~only_iono), fprintf('%s',['Reading RINEX file ' filename ': ... ']); end
                [Eph_J, iono_J] = RINEX_get_nav([filename(1:end-1) 'q'], constellations);
                if(~only_iono), fprintf('done\n'); end
            elseif (~flag_mixed)
                fprintf('... WARNING: QZSS navigation file not found. Disabling QZSS positioning. \n');
                constellations.QZSS.enabled = 0;
            end
        end
        
        if (~only_iono)
            Eph = [Eph_G Eph_R Eph_E Eph_C Eph_J];
        end
        
        if (any(iono_G))
            iono = iono_G;
        elseif (any(iono_R))
            iono = iono_R;
        elseif (any(iono_E))
            iono = iono_E;
        elseif (any(iono_C))
            iono = iono_C;
        elseif (any(iono_J))
            iono = iono_J;
        else
            iono = zeros(8,1);
            fprintf('... WARNING: Klobuchar ionosphere parameters not found in navigation file(s).\n');
        end
        
        if (wait_dlg_PresenceFlag)
            waitbar(1,wait_dlg)
        end
    end
end