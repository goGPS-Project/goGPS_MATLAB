%   OBJECT goO
% =========================================================================
%
% DESCRIPTION:
%   Object to read and collect observation of multiple receivers
%   (including master station)
%
% EXAMPLE
%   goObs = goO();
%
%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%----------------------------------------------------------------------------------------------
%
%    Code contributed by Andrea Gatti
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
%---------------------------------------------------------------------------------------------

classdef goO < handle
    
    % Constant values
    % => to discriminate them from function (in autocompletition) they are
    % written in capital letters
    properties (Constant)

    end

    properties (GetAccess = 'private', SetAccess = 'private')
        constellations;     % Structure to handle constellation parameters
    end

    % Creator (empty)
    methods
        % Object to manage a useful functions / standard parameters of goGPS
        function obj = goO()
        end
    end
    
    %   COMPATIBILITY FUNCTION (STATIC)
    % -------------------------------------------------------------------------
    % function to keep compatibility with the past
    methods (Access = 'public')
        function [pr1_R, pr1_M, ph1_R, ph1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
                dop1_R, dop1_M, dop2_R, dop2_M, snr1_R, snr1_M, ...
                snr2_R, snr2_M, time, time_R, time_M, week_R, week_M, ...
                date_R, date_M, pos_R, pos_M, Eph, iono, interval] = ...
                loadRINEX(obj, filename_nav, filename_R_obs, filename_M_obs, constellations, flag_SP3, wait_dlg)
            
            % Check the input arguments
            if (nargin < 6)
                wait_dlg_PresenceFlag = false;
            else
                wait_dlg_PresenceFlag = true;
            end
            if (isempty(filename_M_obs))
                filename_M_obs_PresenceFlag = false;
            else
                filename_M_obs_PresenceFlag = true;
            end
            [Eph, iono] = obj.readNavigational(filename_nav, constellations, flag_SP3); tic;
            
            %number of satellite slots for enabled constellations
            nSatTot = constellations.nEnabledSat;
            
            %fraction of INTERVAL (epoch-to-epoch timespan, as specified in the header)
            %that is allowed as maximum difference between rover and master timings
            %during synchronization
            max_desync_frac = 0.1;
            
            %-------------------------------------------------------------------------------
            
            % Read Master file if needed
            if (filename_M_obs_PresenceFlag)
                fprintf('Reading MASTER file:\n');
                [obsM intervalM knownPosM dateM] = obj.readRin(filename_M_obs, constellations);
            else
                intervalM = 1e10;
                obsM = 0;
            end
            [obsR intervalR knownPosR dateR] = obj.readRin(filename_R_obs, constellations);
            interval = min([intervalR, intervalM]);
           
            % Convert time in a format that goGPS is actually able to handle:
            timeGPS = datenum(dateR(1,:),dateR(2,:),dateR(3,:),dateR(4,:),dateR(5,:),dateR(6,:))' - datenum(1980,1,6,0,0,0);
            % In the future I just want to keep timeGPS
            [week, tow] = date2gps([dateR(1,:)',dateR(2,:)',dateR(3,:)',dateR(4,:)',dateR(5,:)',dateR(6,:)']);
            [time] = weektow2time(week, tow, 'G');
            
            % Compute empirical interval
            intervalEmp = datevec(median(diff(timeGPS))); intervalEmp = intervalEmp(6) + intervalEmp(5)*60 +  intervalEmp(4)*3600;
            if abs(intervalEmp-interval) > (max_desync_frac*interval)
                interval = intervalEmp;
            end
            intervalGPS = datenum(0,0,0,0,0,interval);
            
            % sSync data;
            startTime = timeGPS(1);
            stopTime = timeGPS(end);
            timeRef = (0:intervalGPS:stopTime-startTime)';                 % The starting epoch here is 0, this time is smoothed  
            nEpochs = length(timeRef);
            
            inRef = round((timeGPS-startTime)/intervalGPS)+1;
            time = zeros(length(timeRef),1); time(inRef) = timeGPS;        % Real time as seen by the receiver

            %-------------------------------------------------------------------------------
            
            % Allocate variables
            if (filename_M_obs_PresenceFlag)
                inRefM = round((timeGPS-startTime)/intervalGPS)+1;

                time_M = zeros(nEpochs,1);
                pr1_M = zeros(nSatTot,nEpochs);
                pr2_M = zeros(nSatTot,nEpochs);
                ph1_M = zeros(nSatTot,nEpochs);
                ph2_M = zeros(nSatTot,nEpochs);
                snr1_M = zeros(nSatTot,nEpochs);
                snr2_M = zeros(nSatTot,nEpochs);
                dop1_M = zeros(nSatTot,nEpochs);
                dop2_M = zeros(nSatTot,nEpochs);
                date_M = zeros(nEpochs,6);
            else
                time_M = []; %  zeros(nEpochs,1);
                pr1_M = [];  %  zeros(nSatTot,nEpochs);
                pr2_M = [];  %  zeros(nSatTot,nEpochs);
                ph1_M = [];  %  zeros(nSatTot,nEpochs);
                ph2_M = [];  %  zeros(nSatTot,nEpochs);
                snr1_M = []; %  zeros(nSatTot,nEpochs);
                snr2_M = []; %  zeros(nSatTot,nEpochs);
                dop1_M = []; %  zeros(nSatTot,nEpochs);
                dop2_M = []; %  zeros(nSatTot,nEpochs);
                date_M = []; %  zeros(nEpochs,6);                
            end
            
            time_R = zeros(nEpochs,1);
            pr1_R = zeros(nSatTot,nEpochs);
            pr2_R = zeros(nSatTot,nEpochs);
            ph1_R = zeros(nSatTot,nEpochs);
            ph2_R = zeros(nSatTot,nEpochs);
            dop1_R = zeros(nSatTot,nEpochs);
            dop2_R = zeros(nSatTot,nEpochs);
            snr1_R = zeros(nSatTot,nEpochs);
            snr2_R = zeros(nSatTot,nEpochs);
            date_R = zeros(nEpochs,6);
            
            %-------------------------------------------------------------------------------
            
            % Sync times
            
            %-------------------------------------------------------------------------------
                        
        end        
    end
    
    %   CONSTELLATION FUNCTION (STATIC)
    % -------------------------------------------------------------------------
    methods (Access = 'public')
        function [constellations] = initConstellation(obj, GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag)
            
            % SYNTAX:
            %   [constellations] = initConstellation(GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag);
            %
            % INPUT:
            %   GPS_flag = boolean flag for enabling/disabling GPS usage
            %   GLO_flag = boolean flag for enabling/disabling GLONASS usage
            %   GAL_flag = boolean flag for enabling/disabling Galileo usage
            %   BDS_flag = boolean flag for enabling/disabling BeiDou usage
            %   QZS_flag = boolean flag for enabling/disabling QZSS usage
            %   SBS_flag = boolean flag for enabling/disabling SBAS usage (for ranging)
            %
            % OUTPUT:
            %   constellations = struct with multi-constellation settings
            %
            % DESCRIPTION:
            %   Multi-constellation settings and initialization.
            
            GPS_PRN = 1:32;
            GLO_PRN = 1:24;
            GAL_PRN = 1:30;
            BDS_PRN = 1:37;
            QZS_PRN = 193:196;
            SBS_PRN = 0; %SBAS ranging not supported yet
            
            constellations.GPS     = struct('numSat', numel(GPS_PRN), 'enabled', GPS_flag, 'indexes', 0, 'PRN', GPS_PRN);
            constellations.GLONASS = struct('numSat', numel(GLO_PRN), 'enabled', GLO_flag, 'indexes', 0, 'PRN', GLO_PRN);
            constellations.Galileo = struct('numSat', numel(GAL_PRN), 'enabled', GAL_flag, 'indexes', 0, 'PRN', GAL_PRN);
            constellations.BeiDou  = struct('numSat', numel(BDS_PRN), 'enabled', BDS_flag, 'indexes', 0, 'PRN', BDS_PRN);
            constellations.QZSS    = struct('numSat', numel(QZS_PRN), 'enabled', QZS_flag, 'indexes', 0, 'PRN', QZS_PRN);
            constellations.SBAS    = struct('numSat', numel(SBS_PRN), 'enabled', 0,        'indexes', 0, 'PRN', SBS_PRN); %SBAS ranging not supported yet
            
            nSatTot = 0; %total number of satellites used given the enabled constellations
            q = 0;       %counter for enabled constellations
            
            systems = fieldnames(constellations);
            constellations.indexes = [];
            constellations.PRN = [];
            for i = 1 : numel(systems)
                if(constellations.(systems{i}).enabled)
                    nSatTot = nSatTot + constellations.(systems{i}).numSat;
                    q = q + 1;
                    if (q == 1)
                        indexes_tmp = [1 : constellations.(systems{i}).numSat];
                    else
                        indexes_tmp = [indexes_tmp(end) + 1 : indexes_tmp(end) + constellations.(systems{i}).numSat];
                    end
                    constellations.(systems{i}).indexes = indexes_tmp;
                    constellations.indexes = [constellations.indexes, indexes_tmp];
                    constellations.PRN = [constellations.PRN, constellations.(systems{i}).PRN];
                end
            end
            
            constellations.nEnabledSat = nSatTot;
            
            obj.constellations = constellations;
        end
    end
    
    %   RINEX LOADER
    % -------------------------------------------------------------------------
    methods (Access = 'public')
                
        % goGNSS.readRin(filename_R_obs, constellations);
        function [obs interval knownPos date] = readRin(obj, fileNameObs, constellations)
            t0=tic;            
            if (isempty(constellations)) %then use only GPS as default
                [constellations] = obj.initConstellation(1, 0, 0, 0, 0, 0);
            end
            % total number of satellites (according to enabled constellations)
            nSatTot = constellations.nEnabledSat;

            % starting index in the total array for the various constellations
            startSysId(goGNSS.ID_SBAS) = uint16(constellations.SBAS.indexes(1)*constellations.SBAS.enabled);
            startSysId(goGNSS.ID_GPS) = uint16(constellations.GPS.indexes(1)*constellations.GPS.enabled);
            startSysId(goGNSS.ID_GLONASS) = uint16(constellations.GLONASS.indexes(1)*constellations.GLONASS.enabled);
            startSysId(goGNSS.ID_GALILEO) = uint16(constellations.Galileo.indexes(1)*constellations.Galileo.enabled);
            startSysId(goGNSS.ID_BEIDOU) = uint16(constellations.BeiDou.indexes(1)*constellations.BeiDou.enabled);
            startSysId(goGNSS.ID_QZSS) = uint16(constellations.QZSS.indexes(1)*constellations.QZSS.enabled);

            %open RINEX observation file (ROVER)
            fprintf('Reading RINEX: %s\n', fileNameObs);
            fid = fopen(fileNameObs,'r');
            txtRin = textscan(fid,'%s','Delimiter','\n','whitespace','');
            fclose(fid);
            txtRin = txtRin{1};
            fprintf('The RINEX file has been read in %f seconds, start parsing...\n', toc(t0));
            
            %parse RINEX header
            [version obsTypes knownPos flagFoundTypes interval sysId line] = obj.RinParseHDR(txtRin);
            
            if ~flagFoundTypes
                error('Basic data is missing in the RINEX header');
            end
            
            [obsCol, nTypes] = obj.getObsTypes(obsTypes, sysId);
                        
            % find first and last epochs
            nObs = obj.getObsSize(txtRin, line, version, interval);
            fprintf('%d epochs should be present in the file\n', nObs);
            
            % init variables
            nSys = 1;   % Start with the hypotesis to have just one constellation available             
            date = zeros(6,nObs);
            obs(nObs) = struct('L1', zeros(nSatTot,1), 'L2', zeros(nSatTot,1), ...
                'C1', zeros(nSatTot,1), ...
                'P1', zeros(nSatTot,1), 'P2', zeros(nSatTot,1), ...
                'S1', zeros(nSatTot,1), 'S2', zeros(nSatTot,1), ...
                'D1', zeros(nSatTot,1), 'D2', zeros(nSatTot,1));
            
            % obs = zeros(nTypes, nObs, nSatTot);
            snr = zeros(nTypes, nObs, nSatTot);
            
            % Loops for reading epochs
            if (version < 3) % RINEX 2
                i=1;
                fprintf('\n           ');
                
                while (line <= length(txtRin)) && (i <= nObs)
                    [date(:,i), nSat, satList, satTypes, line] = obj.RinGetEpoch2(txtRin, line);
                    
                    %[obs(i) line1] = obj.RINEX_get_obs(txtRin, line, nSat, satList, satTypes, obsCol, nTypes, constellations);
                    % c=struct2cell(obs(i)); ccc = [c{1} c{2} c{3} c{4} c{5} c{6} c{7} c{8}];
                    [obs(i) line] = obj.RinGetObs2(txtRin, line, nSat, satList, satTypes, obsCol, nTypes, startSysId, nSatTot, obs(i));
                    % c2=struct2cell(obs(i)); ccc = [c{1} c{2} c{3} c{4} c{5} c{6} c{7} c{8}]-[c2{1} c2{2} c2{3} c2{4} c2{5} c2{6} c2{7} c2{8}];
                    fprintf('\b\b\b\b\b\b\b\b\b\b\b%05d/%05d',i,nObs);
                    i=i+1;
                end
            else
                i=1;
                fprintf('\n           ');
                while (line <= length(txtRin)) && (i <= nObs)
                    [date(:,i), nSat, line] = RinGetEpoch3(txtRin, line);
                    
                    line = line + nSat;
                    fprintf('\b\b\b\b\b\b\b\b\b\b\b%05d/%05d',i,nObs);
                    i=i+1;
                end
            end
            fprintf('\n');           
            toc(t0)            
        end
        
        function [Eph, iono] = readNavigational(obj, filename_nav, constellations, flag_SP3)
            
            % Check the input arguments
            if (isempty(constellations)) % then use only GPS as default
                constellations = obj.initConstellation(1, 0, 0, 0, 0, 0);
            end
            tic;
                        
            %read navigation files
            if (~flag_SP3)
                fprintf('Load navigation file...\n');
                
                Eph_G = []; iono_G = zeros(8,1);
                Eph_R = []; iono_R = zeros(8,1);
                Eph_E = []; iono_E = zeros(8,1);
                Eph_C = []; iono_C = zeros(8,1);
                Eph_J = []; iono_J = zeros(8,1);
                
                if (strcmpi(filename_nav(end),'p'))
                    flag_mixed = 1;
                else
                    flag_mixed = 0;
                end
                
                if (constellations.GPS.enabled || flag_mixed)
                    if (exist(filename_nav,'file'))
                        %parse RINEX navigation file (GPS) NOTE: filename expected to
                        %end with 'n' or 'N' (GPS) or with 'p' or 'P' (mixed GNSS)
                        [Eph_G, iono_G] = RINEX_get_nav(filename_nav, constellations);
                    else
                        fprintf('... WARNING: GPS navigation file not found. Disabling GPS positioning. \n');
                        constellations.GPS.enabled = 0;
                    end
                end
                
                if (constellations.GLONASS.enabled)
                    if (exist([filename_nav(1:end-1) 'g'],'file'))
                        %parse RINEX navigation file (GLONASS)
                        [Eph_R, iono_R] = RINEX_get_nav([filename_nav(1:end-1) 'g'], constellations);
                    elseif (~flag_mixed)
                        fprintf('... WARNING: GLONASS navigation file not found. Disabling GLONASS positioning. \n');
                        constellations.GLONASS.enabled = 0;
                    end
                end
                
                if (constellations.Galileo.enabled)
                    if (exist([filename_nav(1:end-1) 'l'],'file'))
                        %parse RINEX navigation file (Galileo)
                        [Eph_E, iono_E] = RINEX_get_nav([filename_nav(1:end-1) 'l'], constellations);
                    elseif (~flag_mixed)
                        fprintf('... WARNING: Galileo navigation file not found. Disabling Galileo positioning. \n');
                        constellations.Galileo.enabled = 0;
                    end
                end
                
                if (constellations.BeiDou.enabled)
                    if (exist([filename_nav(1:end-1) 'b'],'file'))
                        parse RINEX navigation file (BeiDou)
                        [Eph_C, iono_C] = RINEX_get_nav([filename_nav(1:end-1) 'b'], constellations);
                    elseif (~flag_mixed)
                        fprintf('... WARNING: BeiDou navigation file not found. Disabling BeiDou positioning. \n');
                        constellations.BeiDou.enabled = 0;
                    end
                end
                
                if (constellations.QZSS.enabled)
                    if (exist([filename_nav(1:end-1) 'q'],'file'))
                        %parse RINEX navigation file (QZSS)
                        [Eph_J, iono_J] = RINEX_get_nav([filename_nav(1:end-1) 'q'], constellations);
                    elseif (~flag_mixed)
                        fprintf('... WARNING: QZSS navigation file not found. Disabling QZSS positioning. \n');
                        constellations.QZSS.enabled = 0;
                    end
                end
                
                Eph = [Eph_G Eph_R Eph_E Eph_C Eph_J];
                
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
                    fprintf('... WARNING: ionosphere parameters not found in navigation file(s).\n');
                end
                
                if (wait_dlg_PresenceFlag)
                    waitbar(1,wait_dlg)
                end
            else
                fprintf('Loading SP3 files')
                Eph = zeros(33,nSatTot);
                iono = zeros(8,1);
            end
            obs.eph = Eph;
            obs.iono = iono;
        end
    end
    
    methods (Static, Access = 'public')
        function [version, obsTypes, knownPos, flagFoundTypes, interval, sysId, line] = RinParseHDR(txtRin)
            flagFoundTypes = 0;
            obsTypes = cell(0,0);
            sysId = cell(0,0);
            knownPos = [];
            interval = 1; %default to 1 second (1 Hz observations)
            version = 2;
            
            %parse first line
            line = 1; lin = txtRin{line};
            
            %constellation counter for RINEX v3.xx
            c = 1;
            
            %check if the end of the header or the end of the file has been reached
            while isempty(strfind(lin,'END OF HEADER'))

                answer = strfind(lin,'RINEX VERSION'); %RINEX v2.xx
                if ~isempty(answer)
                    version = floor(sscanf(lin(1:15),'%f'));
                end    
                
                answer = strfind(lin,'# / TYPES OF OBSERV'); %RINEX v2.xx
                if ~isempty(answer)
                    obsTypes{1} = [];
                    nObs = sscanf(lin(1:6),'%d');
                    nLinObs = ceil(nObs/9);
                    for i = 1 : nLinObs
                        if (i > 1)
                            line = 1; lin = txtRin{line};
                        end
                        n = min(nObs,9);
                        for k = 1 : n
                            ot = sscanf(lin(k*6+1:k*6+6),'%s');
                            obsTypes{1} = [obsTypes{1} ot];
                        end
                        nObs = nObs - 9;
                    end
                    
                    flagFoundTypes = 1;
                end
                
                answer = strfind(lin,'SYS / # / OBS TYPES'); %RINEX v3.xx
                if ~isempty(answer)
                    sysId{c} = sscanf(lin(1),'%s');
                    nObs = sscanf(lin(2:6),'%d');
                    obsTypes.(sysId{c}) = [];
                    nLinObs = ceil(nObs/13);
                    for i = 1 : nLinObs
                        if (i > 1)
                            line = 1; lin = txtRin{line};
                        end
                        n = min(nObs,13);
                        for k = 0 : n-1
                            ot = sscanf(lin(6+k*4+1:6+k*4+4),'%s');
                            obsTypes.(sysId{c}) = [obsTypes.(sysId{c}) ot];
                        end
                        nObs = nObs - 13;
                    end
                    
                    c = c + 1;
                    flagFoundTypes = 1;
                end
                
                answer = strfind(lin,'APPROX POSITION XYZ');
                if ~isempty(answer)
                    X = sscanf(lin(1:14),'%f');
                    Y = sscanf(lin(15:28),'%f');
                    Z = sscanf(lin(29:42),'%f');
                    knownPos = [X; Y; Z];
                end
                answer = strfind(lin,'INTERVAL');
                if ~isempty(answer)
                    interval = sscanf(lin(1:10),'%f');
                end
                
                %parse next line
                line = line +1; lin = txtRin{line};
            end
            
            %check RINEX version
            if (~isempty(sysId) && (version == 2))
                version = 3;
            end
        end
        
        function [obsCol, nTypes] = getObsTypes(obsTypes, sysId)
            
            if (isempty(sysId)) %RINEX v2.xx does not have sysId
                
                nTypes = size(obsTypes{1},2)/2;
                
                %search L1 column
                s1 = strfind(obsTypes{1}, 'L1');
                s2 = strfind(obsTypes{1}, 'LA');
                col_L1 = [s1 s2];
                
                %search L2 column
                s1 = strfind(obsTypes{1}, 'L2');
                s2 = strfind(obsTypes{1}, 'LC');
                col_L2 = [s1 s2];
                
                %search C1 column
                s1 = strfind(obsTypes{1}, 'C1');
                s2 = strfind(obsTypes{1}, 'CA');
                col_C1 = [s1 s2];
                
                %search P1 column
                s1 = strfind(obsTypes{1}, 'P1');
                s2 = strfind(obsTypes{1}, 'CA'); %QZSS does not use P1
                col_P1 = [s1 s2];
                
                %if RINEX v2.12 and GPS/GLONASS P1 observations are not available
                if (length(col_P1) ~= 2 && ~isempty(s2))
                    %keep QZSS CA observations as C1
                    col_P1 = [];
                end
                
                %search P2 column
                s1 = strfind(obsTypes{1}, 'P2');
                s2 = strfind(obsTypes{1}, 'CC');
                col_P2 = [s1 s2];
                
                %search S1 column
                s1 = strfind(obsTypes{1}, 'S1');
                s2 = strfind(obsTypes{1}, 'SA');
                col_S1 = [s1 s2];
                
                %search S2 column
                s1 = strfind(obsTypes{1}, 'S2');
                s2 = strfind(obsTypes{1}, 'SC');
                col_S2 = [s1 s2];
                
                %search D1 column
                s1 = strfind(obsTypes{1}, 'D1');
                s2 = strfind(obsTypes{1}, 'DA');
                col_D1 = [s1 s2];
                
                %search D2 column
                s1 = strfind(obsTypes{1}, 'D2');
                s2 = strfind(obsTypes{1}, 'DC');
                col_D2 = [s1 s2];
                
                obsCol.L1 = uint8((col_L1+1)/2);
                obsCol.L2 = uint8((col_L2+1)/2);
                obsCol.C1 = uint8((col_C1+1)/2);
                obsCol.P1 = uint8((col_P1+1)/2);
                obsCol.P2 = uint8((col_P2+1)/2);
                obsCol.S1 = uint8((col_S1+1)/2);
                obsCol.S2 = uint8((col_S2+1)/2);
                obsCol.D1 = uint8((col_D1+1)/2);
                obsCol.D2 = uint8((col_D2+1)/2);
                
            else %RINEX v3.xx
                for c = 1 : length(sysId)
                    
                    nTypes.(sysId{c}) = size(obsTypes.(sysId{c}),2)/3;
                    
                    switch sysId{c}
                        case 'G' %GPS
                            idL1 = 'L1C';
                            idL2 = 'L2W';
                            idC1 = 'C1C';
                            idP1 = 'C1P';
                            idP2 = 'C2W';
                            idS1 = 'S1C';
                            idS2 = 'S2W';
                            idD1 = 'D1C';
                            idD2 = 'D2W';
                        case 'R' %GLONASS
                            idL1 = 'L1C';
                            idL2 = 'L2P';
                            idC1 = 'C1C';
                            idP1 = 'C1P';
                            idP2 = 'C2P';
                            idS1 = 'S1C';
                            idS2 = 'S2P';
                            idD1 = 'D1C';
                            idD2 = 'D2P';
                        case 'E' %Galileo
                            idL1 = 'L1X';
                            idL2 = 'L5X';
                            idC1 = 'C1X';
                            idP1 = '...'; % <-- ?
                            idP2 = 'C5X';
                            idS1 = 'S1X';
                            idS2 = 'S5X';
                            idD1 = 'D1X';
                            idD2 = 'D5X';
                        case 'C' %Compass/Beidou
                            idL1 = 'L2I';
                            idL2 = 'L7I';
                            idC1 = 'C2I';
                            idP1 = '...'; % <-- ?
                            idP2 = 'C7I';
                            idS1 = 'S2I';
                            idS2 = 'S7I';
                            idD1 = 'D2I';
                            idD2 = 'D7I';
                        case 'Q' %QZSS
                            idL1 = 'L1C';
                            idL2 = 'L2C';
                            idC1 = 'C1C';
                            idP1 = '...'; %QZSS does not use P1
                            idP2 = 'C2C';
                            idS1 = 'S1C';
                            idS2 = 'S2C';
                            idD1 = 'D1C';
                            idD2 = 'D2C';
                    end
                    
                    %search L1, L2, C1, P1, P2, S1, S2, D1, D2 columns
                    col_L1 = strfind(obsTypes.(sysId{c}), idL1);
                    col_L2 = strfind(obsTypes.(sysId{c}), idL2);
                    col_C1 = strfind(obsTypes.(sysId{c}), idC1);
                    col_P1 = strfind(obsTypes.(sysId{c}), idP1);
                    col_P2 = strfind(obsTypes.(sysId{c}), idP2);
                    col_S1 = strfind(obsTypes.(sysId{c}), idS1);
                    col_S2 = strfind(obsTypes.(sysId{c}), idS2);
                    col_D1 = strfind(obsTypes.(sysId{c}), idD1);
                    col_D2 = strfind(obsTypes.(sysId{c}), idD2);

                    obsCol.(sysId{c}).L1 = (col_L1+2)/3;
                    obsCol.(sysId{c}).L2 = (col_L2+2)/3;
                    obsCol.(sysId{c}).C1 = (col_C1+2)/3;
                    obsCol.(sysId{c}).P1 = (col_P1+2)/3;
                    obsCol.(sysId{c}).P2 = (col_P2+2)/3;
                    obsCol.(sysId{c}).S1 = (col_S1+2)/3;
                    obsCol.(sysId{c}).S2 = (col_S2+2)/3;
                    obsCol.(sysId{c}).D1 = (col_D1+2)/3;
                    obsCol.(sysId{c}).D2 = (col_D2+2)/3;
                end
            end
        end
        
        function nObs = getObsSize(txtRin, line, version, interval)
            % get RINEX first epoch
            timeGPS1 = [];
            while isempty(timeGPS1) && (line < length(txtRin))
                line = line+1; lin = txtRin{line};
                timeGPS1 = obj.getRinGPSTime(lin, version);
            end
            line = length(txtRin);
            % get RINEX last epoch
            timeGPS2 = [];
            while isempty(timeGPS2) && (line > 1)
                line = line-1; lin = txtRin{line};
                timeGPS2 = obj.getRinGPSTime(lin, version);
            end
            
            % Estimate the number of epoch to consider
            nObs = round(timeGPS2-timeGPS1)/interval;
        end

        function  timeGPS = getRinGPSTime(lin, version)
            % get time epoch for RINEX (format: seconds since 1980/01/06)
            if version >= 3
                timeGPS = obj.getRinGPSTime3(lin);
            else
                timeGPS = obj.getRinGPSTime2(lin);
            end
        end
        
        function  timeGPS = getRinGPSTime2(lin)
            % get time epoch for RINEX v2 (format: seconds since 1980/01/06)       
            data = sscanf(lin(1:min(26,end)),'%f%f%f%f%f%f');
            if length(data) ~= 6
                timeGPS = [];
                return; % This is not a line containing a day
            end            
            %computation of the GPS time in weeks and seconds of week
            year = data(1);
            year = 2000 + year - ((year > 79)*100);
            GPS000 = 723186; % datenum(1980,1,6,0,0,0);
            timeGPS =  (datenum(year,data(2),data(3),data(4),data(5),data(6))-GPS000)*(3600*24);
        end
        
        function  timeGPS = getRinGPSTime3(lin)
            % get time epoch for RINEX v3 (format: seconds since 1980/01/06)
            data = sscanf(lin(2:min(29,end)),'%f%f%f%f%f%f');
            if length(data) ~= 6
                timeGPS = [];
                return; % This is not a line containing a day
            end            
            %computation of the GPS time in weeks and seconds of week
            year = data{1};
            year = 2000 + year - ((year > 79)*100);
            GPS000 = 723186; % datenum(1980,1,6,0,0,0);
            timeGPS =  (datenum(year,data(2),data(3),data(4),data(5),data(6))-GPS000)*(3600*24);
        end
        
        function  date = getRinFullTime(lin, version)
            % get time epoch for RINEX (format (struct): year, month, dat, hour, minute, seconds)
            if version >= 3
                date = obj.getRinFullTime3(lin);
            else
                date = obj.getRinFullTime2(lin);
            end
        end
        
        function  date = getRinFullTime2(lin)
            % get time epoch for RINEX v2 (format: year, month, dat, hour, minute, seconds)
            date = sscanf(lin(1:min(26,end)),'%f%f%f%f%f%f');
            if length(date) ~= 6
                date = [];
                return; % This is not a line containing a day
            end    
            date(1)   = 2000 + date(1) - ((date(1) > 79)*100);
        end
        
        function  date = getRinFullTime3(lin)
            % get time epoch for RINEX v3 (format: year, month, dat, hour, minute, seconds)
            date = sscanf(lin(2:min(29,end)),'%f%f%f%f%f%f');
            if length(date) ~= 6
                date = [];
                return; % This is not a line containing a day
            end    
            date(1)   = 2000 + date(1) - ((date(1) > 79)*100);
        end
        
        function [date, nSat, satList, satTypes, line] = RinGetEpoch2(txtRin, line)
            date = [];
            while isempty(date) && (line < length(txtRin))
                line = line+1; lin = txtRin{line};
                date = obj.getRinFullTime2(lin);
            end
            if (~isempty(date) && (line < length(txtRin)))
                %number of visible satellites
                nSat = sscanf(lin(30:32),'%d');
                
                % The last character to represent satellites is 68 (12 satellites in view)
                % whenever I see more satellitea they are written in the following lines
                nlines = ceil(nSat/12);
                tmp = char(zeros(1,nSat*3));
                n = 1;
                while ( n < (nlines - 2))
                    start = (n-1)*36;
                    tmp(start+1, start+36) = lin(33:68);
                    line = line+1; lin = txtRin{line};
                    n = n+1;
                end
                line = line - (nlines > 1);
                start = (n-1)*36;
                stop = mod(nSat,12)*3;
                tmp(start+1: start + stop) = lin(33:32+stop);
                % tmp now contains the list of satellite seen
                
                pos = 1;
                satList = uint16(zeros(nSat,1));
                satTypes = char(zeros(nSat,1));
                for i = 1 : nSat
                    %check if GPS satellites are labeled 'G' or not labeled
                    if (strcmp(tmp(pos),' '))
                        type = 'G';
                    else
                        type = tmp(pos);
                    end
                    % sat_types = [sat_types; type];
                    satTypes(i) = type;
                    % sat(i) = sscanf(lin(pos+1:pos+2),'%d');
                    satList(i) = mod((tmp(pos+1)-48)*10+(tmp(pos+2)-48),160);
                    pos = pos + 3;
                end
            else
                nSat = [];
                satList = [];
                satTypes = [];
                line = [];
            end
        end
        
        function [date, nSat, line] = RinGetEpoch3(txtRin, line)
            date = [];
            while isempty(date) && (line < length(txtRin))
                line = line+1; lin = txtRin{line};
                date = obj.getRinFullTime3(lin);
            end
            if (~isempty(date) && (line < length(txtRin)))
                %number of visible satellites
                [nSat] = sscanf(lin(33:35),'%d');
            else
                nSat = [];
            end
        end
                        
        function [obs_struct, line] = RinGetObs2(txtRin, line, nSat, sat, sat_types, obs_col, nObsTypes, startSysId, nSatTot, obs_struct)
                                    
            obs_tmp = struct('TMP1', zeros(nSatTot,1), 'TMP2', zeros(nSatTot,1));
            
            obs_struct = struct('L1', zeros(nSatTot,1), 'L2', zeros(nSatTot,1), ...
                'C1', zeros(nSatTot,1), ...
                'P1', zeros(nSatTot,1), 'P2', zeros(nSatTot,1), ...
                'S1', zeros(nSatTot,1), 'S2', zeros(nSatTot,1), ...
                'D1', zeros(nSatTot,1), 'D2', zeros(nSatTot,1));
            
            % array to contain the starting index for each constellation in the total array (with length = nSatTot)
            sat_types_id = zeros(size(sat_types));
            
            % convert constellations letter to starting index in the total array
            satNum = sat_types == 'G';
            sat_types_id(satNum) = startSysId(goGNSS.ID_GPS);
            satNum = sat_types == 'R';
            if(sum(satNum)>0), sat_types_id(satNum) = startSysId(goGNSS.ID_GLONASS); end;
            satNum = sat_types == 'E';
            if(sum(satNum)>0), sat_types_id(satNum) = startSysId(goGNSS.ID_GALILEO); end;
            satNum = sat_types == 'C';
            if(sum(satNum)>0), sat_types_id(satNum) = startSysId(goGNSS.ID_BEIDOU); end;
            satNum = sat_types == 'J';
            if(sum(satNum)>0), sat_types_id(satNum) = startSysId(goGNSS.ID_QZSS); end;
            satNum = sat_types == 'S';
            if(sum(satNum)>0), sat_types_id(satNum) = startSysId(goGNSS.ID_SBAS); end;
            
            % observation types

            nLinesToRead = ceil(nObsTypes/5);  % I read a maximum of 5 obs per line => this is the number of lines to read
            nObsToRead = nLinesToRead * 5;     % Each line contains 5 observations
            
            %data read and assignment
            lin = char(32*uint8(ones(16*nObsToRead,1))'); % preallocate the line to read
            
            % Mask to filter all the possible observations (max 15)
            mask = false(16,nObsToRead);
            mask(2:14,:) = true;
            % preallocate a matrix of 15 strings (of length 14 characters)
            % notice that each observation element has a max length of 13 char,
            % the first character is added as a padding to separate the strings for
            % the command sscanf that now can be launched once for each satellite
            strObs = char(ones(14,nObsToRead)*32);
            
            for s = 1 : nSat
                
                % DEBUG
                if (sat(s) > 32)
                    sat(s) = 32; %this happens only with SBAS; it's already fixed in the multi-constellation version
                end
                
                %lin = char(lin*0+32); % clear line -> fill with spaces
                lin = char(32*uint8(ones(size(lin)))); % clear line -> fill with spaces
                
                if (sat_types_id(s) ~= 0)
                    % read all the lines containing the observations needed
                    for l = 1 : (nLinesToRead)
                        line = line+1; linTmp = txtRin{line};
                        linLengthTmp = length(linTmp);
                        lin((80*(l-1))+(1:linLengthTmp)) = linTmp;  %each line has a maximum lenght of 80 characters
                    end
                    linLength = 80*(nLinesToRead-1)+linLengthTmp;
                    
                    % convert the lines read from the RINEX file to a single matrix
                    % containing all the observations
                    strObs(1:13,:) = (reshape(lin(mask(:)),13,nObsToRead));
                    strSum = sum(strObs);
                    fltObs = sscanf(strObs, '%f'); % read all the observations in the string
                    obsId = 0; % index of the current observation
                    % start parsing the observation string
                    satId = sat_types_id(s)+sat(s)-1;
                    for k = 1 : min(nObsTypes, floor(linLength/16))
                        % check if the element is empty
                        if (strSum(k) ~= 448) % if the current val is not missing (full of spaces)
                            obsId = obsId+1;
                            %obs = sscanf(lin(mask(:,k)), '%f');
                            
                            %check and assign the observation type
                            if k == obs_col.L1
                                obs_struct.L1(satId) = fltObs(obsId);
                                if (linLength>=16*k)
                                    %convert signal-to-noise ratio
                                    % faster conversion of a single ASCII character into an int
                                    snr = mod((lin(16*k)-48),16);
                                    obs_tmp.TMP1(satId) = 6 * snr;
                                end
                            elseif k == obs_col.L2
                                obs_struct.L2(satId) = fltObs(obsId);
                                if (linLength>=16*k)
                                    %convert signal-to-noise ratio
                                    % faster conversion of a single ASCII character into an int
                                    snr = mod((lin(16*k)-48),16);
                                    obs_tmp.TMP2(satId) = 6 * snr;
                                end
                            elseif k == obs_col.C1
                                obs_struct.C1(satId) = fltObs(obsId);
                            elseif k == obs_col.P1
                                obs_struct.P1(satId) = fltObs(obsId);
                            elseif k == obs_col.P2
                                obs_struct.P2(satId) = fltObs(obsId);
                            elseif k == obs_col.D1
                                obs_struct.D1(satId) = fltObs(obsId);
                            elseif k == obs_col.D2
                                obs_struct.D2(satId) = fltObs(obsId);
                            elseif k == obs_col.S1
                                obs_struct.S1(satId) = fltObs(obsId);
                            elseif k == obs_col.S2
                                obs_struct.S2(satId) = fltObs(obsId);
                            end
                        end
                    end
                    if (isempty(obs_struct.S1(satId)))
                        obs_struct.S1(satId) = obs_tmp.TMP1(satId);
                    end
                    if (isempty(obs_struct.S2(satId)))
                        obs_struct.S2(satId) = obs_tmp.TMP2(satId);
                    end
                else
                    % skip all the observation lines for the unused satellite
                    for l = 1 : (nLinesToRead)
                        line = line+1;
                    end
                end
            end
        end
    end
    
end

% load loadRINEXtmp
% tic; [pr1_R, pr1_M, ph1_R, ph1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
% dop1_R, dop1_M, dop2_R, dop2_M, snr1_R, snr1_M, ...
% snr2_R, snr2_M, time, time_R, time_M, week_R, week_M, ...
% date_R, date_M, pos_R, pos_M, Eph, iono, interval] = ...
% load_RINEX(filename_nav, filename_R_obs, '', constellations, flag_SP3); toc;
% %%
% tic; [pr1_R, pr1_M, ph1_R, ph1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
% dop1_R, dop1_M, dop2_R, dop2_M, snr1_R, snr1_M, ...
% snr2_R, snr2_M, time, time_R, time_M, week_R, week_M, ...
% date_R, date_M, pos_R, pos_M, Eph, iono, interval] = ...
% goGNSS.loadRINEX(filename_nav, filename_R_obs, '', constellations, flag_SP3); toc;
