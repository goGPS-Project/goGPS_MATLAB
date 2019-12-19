function [Eph, iono, flag_return] = load_RINEX_nav(filename, cc, flag_SP3, iono_model, time, wait_dlg)

% SYNTAX:
%   [Eph, iono, flag_return] = load_RINEX_nav(filename, constellations, flag_SP3, iono_model, time, wait_dlg);
%
% INPUT:
%   filename = RINEX navigation file
%   cc = Constellation_Collector object, contains the satus of the satellite systems in use
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
%  Contributors:     Damianop Triglione 2012, ...
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

% Check the input arguments
if (nargin < 6)
    wait_dlg_PresenceFlag = false;
else
    wait_dlg_PresenceFlag = true;
end

if (iscell(filename))
    filename = filename{1};
end

flag_return = 0;
log = Logger.getInstance();
state = Core.getCurrentSettings();

%number of satellite slots for enabled constellations
nSatTot = cc.getNumSat();

%read navigation files
if (~flag_SP3)
    parse_file(0);
else
    Eph = zeros(33,nSatTot);
    iono = zeros(8,1);
end

% Broadcast corrections in DD are currently causing problems (offset in UP) => not using them
if state.isModeSA()
    %if Klobuchar ionospheric delay correction is requested but parameters are not available in the navigation file, try to download them
    if ((iono_model == 2 && ~any(iono)) || (flag_SP3 && cc.getGLONASS().isActive()))
        [week, sow] = time2weektow(time(1));
        [date, DOY] = gps2date(week, sow);
        
        filename_brdm = ['brdm' num2str(DOY,'%03d') '0.' num2str(two_digit_year(date(1,1)),'%02d') 'p'];
        filename_brdc = ['brdc' num2str(DOY,'%03d') '0.' num2str(two_digit_year(date(1,1)),'%02d') 'n'];
        filename_CGIM = ['CGIM' num2str(DOY,'%03d') '0.' num2str(two_digit_year(date(1,1)),'%02d') 'N'];
        
        pos = find(filename == '/'); if(isempty(pos)), pos = find(filename == '\'); end
        nav_path = filename(1:pos(end));
        
        flag_GLO = flag_SP3 && cc.getGLONASS().isActive();
        
        file_avail = 0;
        if (exist([nav_path filename_brdm],'file') && flag_GLO)
            filename = [nav_path filename_brdm];
            file_avail = 1;
        elseif (exist([nav_path filename_CGIM],'file') && ~flag_GLO)
            filename = [nav_path filename_CGIM];
            file_avail = 1;
        elseif (exist([nav_path filename_brdc],'file') && ~flag_GLO)
            filename = [nav_path filename_brdc];
            file_avail = 1;
        else
            if (flag_GLO)
                filename = filename_brdm;
            else
                filename = filename_brdc;
            end
            [download_successful, compressed] = download_nav(filename, nav_path);
            filename = [nav_path filename];
            if (download_successful)
                file_avail = 1;
            end
            if (compressed)
                flag_return = 1;
            end
        end
        
        if (file_avail)
            if (flag_GLO)
                only_iono = 0;
            else
                only_iono = 1;
            end
            parse_file(only_iono);
        end
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
        Eph_I = []; iono_I = zeros(8,1);

        if (strcmpi(filename(end),'p'))
            flag_mixed = 1;
        else
            flag_mixed = 0;
        end

        if (cc.getGPS().isActive() || flag_mixed || only_iono)
            if (exist(filename,'file'))
                %parse RINEX navigation file (GPS) NOTE: filename expected to
                %end with 'n' or 'N' (GPS) or with 'p' or 'P' (mixed GNSS)
                if(~only_iono), log.addMessage(sprintf('%s',['Reading RINEX file ' filename ': ... '])); end
                [Eph_G, iono_G] = RINEX_get_nav(filename, cc);
                if(~only_iono), log.addStatusOk(); end
            else
                log.addWarning('GPS navigation file not found. Disabling GPS positioning. \n');
                cc.deactivateGPS();
            end
        end

        if (cc.getGLONASS().isActive() && ~only_iono)
            if (exist([filename(1:end-1) 'g'],'file'))
                %parse RINEX navigation file (GLONASS)
                if(~only_iono), log.addMessage(sprintf('%s',['Reading RINEX file ' filename(1:end-1) 'g: ... '])); end
                [Eph_R, iono_R] = RINEX_get_nav([filename(1:end-1) 'g'], cc);
                if(~only_iono), log.addStatusOk(); end
            elseif (~flag_mixed)
                log.addWarning('GLONASS navigation file not found. Disabling GLONASS positioning. \n');
                cc.deactivateGLONASS();
            end
        end

        if (cc.getGalileo().isActive() && ~only_iono)
            if (exist([filename(1:end-1) 'l'],'file'))
                %parse RINEX navigation file (Galileo)
                if(~only_iono), log.addMessage(sprintf('%s',['Reading RINEX file ' filename(1:end-1) 'l: ... '])); end
                [Eph_E, iono_E] = RINEX_get_nav([filename(1:end-1) 'l'], cc);
                if(~only_iono), log.addStatusOk(); end
            elseif (~flag_mixed)
                log.addWarning('Galileo navigation file not found. Disabling Galileo positioning. \n');
                cc.deactivateGalileo();
            end
        end

        if (cc.getBeiDou().isActive() && ~only_iono)
            if (exist([filename(1:end-1) 'c'],'file'))
                %parse RINEX navigation file (BeiDou)
                if(~only_iono), log.addMessage(sprintf('%s',['Reading RINEX file ' filename(1:end-1) 'c: ... '])); end
                [Eph_C, iono_C] = RINEX_get_nav([filename(1:end-1) 'c'], cc);
                if(~only_iono), log.addStatusOk(); end
            elseif (~flag_mixed)
                log.addWarning('BeiDou navigation file not found. Disabling BeiDou positioning. \n');
                cc.deactivateBeiDou();
            end
        end

        if (cc.getQZSS().isActive() && ~only_iono)
            if (exist([filename(1:end-1) 'q'],'file'))
                %parse RINEX navigation file (QZSS)
                if(~only_iono), log.addMessage(sprintf('%s',['Reading RINEX file ' filename(1:end-1) 'q: ... '])); end
                [Eph_J, iono_J] = RINEX_get_nav([filename(1:end-1) 'q'], cc);
                if(~only_iono), log.addStatusOk(); end
            elseif (~flag_mixed)
                log.addWarning('QZSS navigation file not found. Disabling QZSS positioning. \n');
                cc.deactivateQZSS();
            end
        end
        
        if (cc.getIRNSS().isActive() && ~only_iono)
            if (exist([filename(1:end-1) 'i'],'file'))
                %parse RINEX navigation file (IRNSS)
                if(~only_iono), log.addMessage(sprintf('%s',['Reading RINEX file ' filename(1:end-1) 'q: ... '])); end
                [Eph_I, iono_I] = RINEX_get_nav([filename(1:end-1) 'i'], cc);
                if(~only_iono), log.addStatusOk(); end
            elseif (~flag_mixed)
                log.addWarning('IRNSS navigation file not found. Disabling QZSS positioning. \n');
                cc.deactivateIRNSS();
            end
        end
        
        if (~only_iono)
            Eph = [Eph_G Eph_R Eph_E Eph_C Eph_J Eph_I];
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
        elseif (any(iono_I))
            iono = iono_I;
        else
            iono = zeros(8,1);
            log.addWarning('Klobuchar ionosphere parameters not found in navigation file(s).\n');
        end

        if (wait_dlg_PresenceFlag)
            waitbar(1,wait_dlg)
        end
    end
end
