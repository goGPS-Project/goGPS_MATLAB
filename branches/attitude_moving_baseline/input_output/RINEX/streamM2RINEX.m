function streamM2RINEX(fileroot, path, week, rinex_metadata, constellations, wait_dlg)

% SYNTAX:
%   streamM2RINEX(fileroot, path, week, rinex_metadata, constellations, wait_dlg);
%
% INPUT:
%   fileroot = input file root (master data, binary stream)
%   path = output path
%   week = GPS week
%   rinex_metadata = struct with RINEX metadata
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%   wait_dlg = optional handler to waitbar figure
%
% OUTPUT:
%
% DESCRIPTION:
%   File conversion from master stream (RTCM 3.x) to RINEX format.

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

nSatTot = constellations.nEnabledSat;
if (nSatTot == 0)
    fprintf('No constellations selected, setting default: GPS-only\n');
    [constellations] = goGNSS.initConstellation(1, 0, 0, 0, 0, 0);
    nSatTot = constellations.nEnabledSat;
end

%check requested RINEX version
if (strcmp(rinex_metadata.version(1),'3'))
    rin_ver_id = 3;
else
    rin_ver_id = 2;
end

if (nargin == 6)
    waitbar(0.5,wait_dlg,'Reading master stream files...')
end

%check if the system uses 3-digit exponential notation
three_digit_exp = (length(sprintf('%1.1E',1)) == 8);

%MASTER stream reading
data_master_all = [];                                                %overall stream
hour = 0;                                                            %hour index (integer)
hour_str = num2str(hour,'%03d');                                     %hour index (string)
d = dir([fileroot '_master_' hour_str '.bin']);                      %file to be read
while ~isempty(d)
    if (nargin == 5)
        fprintf(['Reading: ' fileroot '_master_' hour_str '.bin\n']);
    end
    num_bytes = d.bytes;                                             %file size (number of bytes)
    fid_master = fopen([fileroot '_master_' hour_str '.bin']);       %file opening
    data_master = fread(fid_master,num_bytes,'uint8');               %file reading
    data_master = dec2bin(data_master,8);                            %conversion into binary number (N x 8bits matrix)
    data_master = data_master';                                      %transposed (8bits x N matrix)
    data_master = data_master(:)';                                   %conversion into a string (8N bits vector)
    fclose(fid_master);                                              %file closing
    data_master_all = [data_master_all data_master];                 %stream concatenation
    hour = hour+1;                                                   %hour increase
    hour_str = num2str(hour,'%03d');
    d = dir([fileroot '_master_' hour_str '.bin']);                  %file to be read
end

%MASTER stream reading (non-goGPS binary format)
d = dir(fileroot);                                                   %file to be read
if ~isempty(d)
    if (nargin == 1)
        fprintf(['Reading: ' fileroot '\n']);
    end
    num_bytes = d.bytes;                                             %file size (number of bytes)
    fid_master = fopen(fileroot);                                    %file opening
    data_master_all = fread(fid_master,num_bytes,'uint8');           %file reading
    data_master_all = dec2bin(data_master_all,8);                    %conversion in binary number (N x 8bits matrix)
    data_master_all = data_master_all';                              %transposed (8bits x N matrix)
    data_master_all = data_master_all(:)';                           %conversion into a string (8N bits vector)
    fclose(fid_master);                                              %file closing
end

clear hour hour_str d
clear data_master fid_master

if (nargin == 6)
    waitbar(1,wait_dlg)
end

%----------------------------------------------------------------------------------------------

if (~isempty(data_master_all))
    
    %displaying
    if (nargin == 5)
        fprintf('Decoding master data \n');
    end
    
    pos = 1;
    sixofeight = [];
    is_rtcm2 = 1;
    
    while (pos + 7 <= length(data_master_all))
        if (~strcmp(data_master_all(pos:pos+1),'01'))
            is_rtcm2 = 0;
            break
        end
        sixofeight = [sixofeight fliplr(data_master_all(pos+2:pos+7))];
        pos = pos + 8;
    end
    
    %stream decodification
    if(is_rtcm2)
        error('RTCM2.x conversion not supported yet!');
        %     [cell_master] = decode_rtcm2(sixofeight);
    else
        if (nargin == 6)
            [cell_master] = decode_rtcm3(data_master_all, constellations, wait_dlg);
        else
            [cell_master] = decode_rtcm3(data_master_all, constellations);
        end
    end
    clear data_master_all
    
    if (~isempty(cell_master))
        
        %initialization (to make the writing faster)
        Ncell   = size(cell_master,2);                        %number of read RTCM packets
        time_M  = zeros(Ncell,1);                             %GPS time
        pr1_M   = zeros(nSatTot,Ncell);                       %code observations
        ph1_M   = zeros(nSatTot,Ncell);                       %phase observations
        snr1_M  = zeros(nSatTot,Ncell);                       %signal-to-noise ratio
        pr2_M   = zeros(nSatTot,Ncell);                       %code observations
        ph2_M   = zeros(nSatTot,Ncell);                       %phase observations
        snr2_M  = zeros(nSatTot,Ncell);                       %signal-to-noise ratio
        Eph_M   = zeros(33,nSatTot,Ncell);                    %ephemerides
        pos_M   = zeros(3,1);                                 %master station position
        
        if (nargin == 6)
            waitbar(0,wait_dlg,'Reading master data...')
        end
        
        flag_L2 = 0;
        
        i = 1;
        toe_old = -1;
        time_M_old = -1;
        flag_no_obs = 1;
        week_M = ones(Ncell,1)*week;
        week_cycle = 0;
        for j = 1 : Ncell
            if (nargin == 6)
                waitbar(j/Ncell,wait_dlg)
            end
            
            if (cell_master{1,j} == 1002)                     %RTCM 1002 message
                
                idx = constellations.GPS.indexes;
                flag_no_obs = 1;
                
                time_M(i)     = cell_master{2,j}(2);          %GPS time logging
                pr1_M(idx,i)  = cell_master{3,j}(idx,2);      %code observations logging
                ph1_M(idx,i)  = cell_master{3,j}(idx,3);      %phase observations logging
                snr1_M(idx,i) = cell_master{3,j}(idx,5);      %signal-to-noise ratio logging
                
                if (time_M(i) < time_M_old)
                    week_cycle = week_cycle + 1;
                end
                
                week_M(i) = week_M(i) + week_cycle;
                time_M_old = time_M(i);
                
                i = i+1;                                      %epoch counter increase
                
                if (cell_master{3,j}(idx,1) == 0)
                    code_type = 'C1';
                else
                    code_type = 'P1';
                end
                
            elseif (cell_master{1,j} == 1004)                 %RTCM 1004 message
                
                idx = constellations.GPS.indexes;
                flag_no_obs = 1;
                
                time_M(i)     = cell_master{2,j}(2);          %GPS time logging
                pr1_M(idx,i)  = cell_master{3,j}(idx,2);      %code observations logging (L1)
                ph1_M(idx,i)  = cell_master{3,j}(idx,3);      %phase observations logging (L1)
                snr1_M(idx,i) = cell_master{3,j}(idx,5);      %signal-to-noise ratio logging (L1)
                pr2_M(idx,i)  = cell_master{3,j}(idx,7);      %code observations logging (L2)
                ph2_M(idx,i)  = cell_master{3,j}(idx,8);      %phase observations logging (L2)
                snr2_M(idx,i) = cell_master{3,j}(idx,10);     %signal-to-noise ratio logging (L2)
                
                flag_L2 = 1;
                
                if (time_M(i) < time_M_old)
                    week_cycle = week_cycle + 1;
                end
                
                week_M(i) = week_M(i) + week_cycle;
                time_M_old = time_M(i);
                
                i = i+1;
                
                if (cell_master{3,j}(idx,1) == 0)
                    code_type = 'C1';
                else
                    code_type = 'P1';
                end
                
            elseif ((cell_master{1,j} == 1005) | (cell_master{1,j} == 1006)) & (pos_M == 0)
                
                coordX_M = cell_master{2,j}(8);
                coordY_M = cell_master{2,j}(9);
                coordZ_M = cell_master{2,j}(10);
                
                pos_M(:,1) = [coordX_M; coordY_M; coordZ_M];
                
            elseif (cell_master{1,j} == 1010)                 %RTCM 1010 message
                
                idx = constellations.GLONASS.indexes;
                
                %time_M(i)     = cell_master{2,j}(2);          %GLONASS time logging
                pr1_M(idx,i)  = cell_master{3,j}(idx,2);      %code observations logging (L1)
                ph1_M(idx,i)  = cell_master{3,j}(idx,3);      %phase observations logging (L1)
                snr1_M(idx,i) = cell_master{3,j}(idx,5);      %signal-to-noise ratio logging (L1)
                
                %flag_L2 = 1;
                
                %i = i+1;
                
                %if (cell_master{3,j}(:,1) == 0)
                %    code_type = 'C1';
                %else
                %    code_type = 'P1';
                %end
                
            elseif (cell_master{1,j} == 1012)                 %RTCM 1012 message
                
                idx = constellations.GLONASS.indexes;
                
                %time_M(i)     = cell_master{2,j}(2);          %GLONASS time logging
                pr1_M(idx,i)  = cell_master{3,j}(idx,2);      %code observations logging (L1)
                ph1_M(idx,i)  = cell_master{3,j}(idx,3);      %phase observations logging (L1)
                snr1_M(idx,i) = cell_master{3,j}(idx,5);      %signal-to-noise ratio logging (L1)
                pr2_M(idx,i)  = cell_master{3,j}(idx,7);      %code observations logging (L2)
                ph2_M(idx,i)  = cell_master{3,j}(idx,8);      %phase observations logging (L2)
                snr2_M(idx,i) = cell_master{3,j}(idx,10);     %signal-to-noise ratio logging (L2)
                
                %flag_L2 = 1;
                
                %i = i+1;
                
                %if (cell_master{3,j}(:,1) == 0)
                %    code_type = 'C1';
                %else
                %    code_type = 'P1';
                %end
                
            elseif (cell_master{1,j} == 1019)                 %RTCM 1019 message
                
                %satellite number
                sat = cell_master{2,j}(30);                   %satellite number
                toe = cell_master{2,j}(18);                   %time of ephemeris
                
                %update the epoch counter (in case no observation messages
                %are available)
                if (flag_no_obs && (toe ~= toe_old))
                    toe_old = toe;
                    i = i + 1;
                end
                
                %if the ephemerides are not already available
                if (isempty(find(Eph_M(18,sat,:) ==  toe, 1)))
                    Eph_M(:,sat,i) = cell_master{2,j}(:);     %single satellite ephemerides logging
                end
                
            elseif (cell_master{1,j} == 1020)                 %RTCM 1020 message
                
                %satellite number
                sat = cell_master{2,j}(30);                   %satellite number
                toe = cell_master{2,j}(18);                   %time of ephemeris
                
                %if the ephemerides are not already available
                if (isempty(find(Eph_M(18,sat,:) ==  toe, 1)))
                    Eph_M(:,sat,i) = cell_master{2,j}(:);     %single satellite ephemerides logging
                end
            end
        end
        clear Ncell pos sat toe
        
        %residual data erase (after initialization)
        time_M(i:end)   = [];
        week_M(i:end)   = [];
        pr1_M(:,i:end)  = [];
        ph1_M(:,i:end)  = [];
        snr1_M(:,i:end) = [];
        pr2_M(:,i:end)  = [];
        ph2_M(:,i:end)  = [];
        snr2_M(:,i:end) = [];
        Eph_M(:,:,i:end) = [];
        
        %manage "nearly null" data
        ph1_M(ph1_M < 1e-100) = 0;
        ph2_M(ph2_M < 1e-100) = 0;
        
        if (any(time_M))
            %date decoding
            [date, DOY] = gps2date(week_M, time_M);
            
            [DOYs, DOY_pos] = unique(DOY, 'first');
            n_doy = length(DOYs);
        else
            %displaying
            if (nargin == 6)
                msgbox('No raw data acquired.');
            else
                fprintf('No raw data acquired.\n');
            end
            
            return
        end
        
        for d = 1 : n_doy
            
            %epoch indexes for each DOY
            first_epoch = DOY_pos(d);
            if (d < n_doy)
                last_epoch  = DOY_pos(d+1)-1;
            else
                last_epoch = length(time_M);
            end
            
            %----------------------------------------------------------------------------------------------
            % OBSERVATION RATE (INTERVAL)
            %----------------------------------------------------------------------------------------------

            interval = median(time_M(first_epoch+1:last_epoch) - time_M(first_epoch:last_epoch-1));
            
            %----------------------------------------------------------------------------------------------
            % MULTI-CONSTELLATION SETTINGS
            %----------------------------------------------------------------------------------------------
            
            GPSenabled = constellations.GPS.enabled;
            GLOenabled = constellations.GLONASS.enabled;
            GALenabled = constellations.Galileo.enabled;
            BDSenabled = constellations.BeiDou.enabled;
            QZSenabled = constellations.QZSS.enabled;
            SBSenabled = constellations.SBAS.enabled;
            
            if (GPSenabled), GPSavailable = any(any(pr1_M(constellations.GPS.indexes,:))); else GPSavailable = 0; end
            if (GLOenabled), GLOavailable = any(any(pr1_M(constellations.GLONASS.indexes,:))); else GLOavailable = 0; end
            if (GALenabled), GALavailable = any(any(pr1_M(constellations.Galileo.indexes,:))); else GALavailable = 0; end
            if (BDSenabled), BDSavailable = any(any(pr1_M(constellations.BeiDou.indexes,:))); else BDSavailable = 0; end
            if (QZSenabled), QZSavailable = any(any(pr1_M(constellations.QZSS.indexes,:))); else QZSavailable = 0; end
            if (SBSenabled), SBSavailable = any(any(pr1_M(constellations.SBAS.indexes,:))); else SBSavailable = 0; end
            
            Eph_ind_array = Eph_M(30,:,:);
            Eph_ind_array = Eph_ind_array(:);
            if (GPSenabled), GPSavailableEPH = ~isempty(intersect(Eph_ind_array,constellations.GPS.indexes)); else GPSavailableEPH = 0; end
            if (GLOenabled), GLOavailableEPH = ~isempty(intersect(Eph_ind_array,constellations.GLONASS.indexes)); else GLOavailableEPH = 0; end
            if (GALenabled), GALavailableEPH = ~isempty(intersect(Eph_ind_array,constellations.Galileo.indexes)); else GALavailableEPH = 0; end
            if (BDSenabled), BDSavailableEPH = ~isempty(intersect(Eph_ind_array,constellations.BeiDou.indexes)); else BDSavailableEPH = 0; end
            if (QZSenabled), QZSavailableEPH = ~isempty(intersect(Eph_ind_array,constellations.QZSS.indexes)); else QZSavailableEPH = 0; end
            if (SBSenabled), SBSavailableEPH = ~isempty(intersect(Eph_ind_array,constellations.SBAS.indexes)); else SBSavailableEPH = 0; end
            
            i = 1;
            if (GPSenabled && ~GPSavailable), warning_msg{i} = 'GPS observations not available'; i = i + 1; end
            if (GLOenabled && ~GLOavailable), warning_msg{i} = 'GLONASS observations not available'; i = i + 1; end
            if (GALenabled && ~GALavailable), warning_msg{i} = 'Galileo observations not available'; i = i + 1; end
            if (BDSenabled && ~BDSavailable), warning_msg{i} = 'BeiDou observations not available'; i = i + 1; end
            if (QZSenabled && ~QZSavailable), warning_msg{i} = 'QZSS observations not available'; i = i + 1; end
            if (SBSenabled && ~SBSavailable), warning_msg{i} = 'SBAS observations not available'; i = i + 1; end
            
            if (GPSenabled && GPSavailable && ~GPSavailableEPH), warning_msg{i} = 'GPS ephemerides not available'; i = i + 1; end
            if (GLOenabled && GLOavailable && ~GLOavailableEPH), warning_msg{i} = 'GLONASS ephemerides not available'; i = i + 1; end
            if (GALenabled && GALavailable && ~GALavailableEPH), warning_msg{i} = 'Galileo ephemerides not available'; i = i + 1; end
            if (BDSenabled && BDSavailable && ~BDSavailableEPH), warning_msg{i} = 'BeiDou ephemerides not available'; i = i + 1; end
            if (QZSenabled && QZSavailable && ~QZSavailableEPH), warning_msg{i} = 'QZSS ephemerides not available'; i = i + 1; end
            if (SBSenabled && SBSavailable && ~SBSavailableEPH), warning_msg{i} = 'SBAS ephemerides not available'; i = i + 1; end
            
            if (i > 1 && n_doy == 1)
                if (nargin == 5)
                    msgbox(warning_msg);
                else
                    for j = 1 : (i-1); fprintf(['WARNING: ' warning_msg{j} '.\n']); end
                end
            end
            
            GPSactive = (GPSenabled && GPSavailable);
            GLOactive = (GLOenabled && GLOavailable);
            GALactive = (GALenabled && GALavailable);
            BDSactive = (BDSenabled && BDSavailable);
            QZSactive = (QZSenabled && QZSavailable);
            SBSactive = (SBSenabled && SBSavailable);
            
            GPSactiveEPH = (GPSenabled && GPSavailable && GPSavailableEPH);
            GLOactiveEPH = (GLOenabled && GLOavailable && GLOavailableEPH);
            GALactiveEPH = (GALenabled && GALavailable && GALavailableEPH);
            BDSactiveEPH = (BDSenabled && BDSavailable && BDSavailableEPH);
            QZSactiveEPH = (QZSenabled && QZSavailable && QZSavailableEPH);
            SBSactiveEPH = (SBSenabled && SBSavailable && SBSavailableEPH);
            
            SYSactive = [GPSactive, GLOactive, GALactive, ...
                BDSactive, QZSactive, SBSactive];
            
            SYSactiveEPH = [GPSactiveEPH, GLOactiveEPH, GALactiveEPH, ...
                BDSactiveEPH, QZSactiveEPH, SBSactiveEPH];
            
            mixed_sys  = (sum(SYSactive) >= 2);
            single_sys = (sum(SYSactive) == 1);
            
            mixed_sys_eph  = (sum(SYSactiveEPH) >= 2);
            single_sys_eph = (sum(SYSactiveEPH) == 1);
            
            %----------------------------------------------------------------------------------------------
            % RINEX OBSERVATION FILE
            %----------------------------------------------------------------------------------------------
            
            %displaying
            if (nargin == 5)
                fprintf('Writing master observation file...\n');
            end
            
            if (nargin == 6)
                waitbar(0,wait_dlg,'Writing master observation file...')
            end
            
            %maximum length for MARKER NAME in the RINEX filename
            mn_len = min(4,length(rinex_metadata.marker_name));
            marker = rinex_metadata.marker_name(1:mn_len);
            for r = 1 : (4-mn_len)
                %fill in remaining letters with 'X'
                marker = [marker 'X']; %#ok<AGROW>
            end
            
            %create RINEX observation file
            fid_obs = fopen([path marker sprintf('%03d', DOYs(d)) '0' sprintf('.%2do', two_digit_year(date(first_epoch,1)))],'wt');
            
            %file type
            if (mixed_sys)
                file_type = 'M: Mixed';
            elseif(single_sys && GPSactive)
                file_type = 'G: GPS';
            elseif(single_sys && GLOactive)
                file_type = 'R: GLONASS';
            elseif(single_sys && GALactive)
                file_type = 'E: Galileo';
            elseif(single_sys && SBSactive)
                file_type = 'S: SBAS Payload';
            end
            
            %current date and time in UTC
            date_str_UTC = local_time_to_utc(now, 30);
            [date_str_UTC, time_str_UTC] = strtok(date_str_UTC,'T');
            time_str_UTC = time_str_UTC(2:end-1);
            
            %maximum filename length for RUN BY field
            rb_len = min(20,length(rinex_metadata.agency));
            
            %maximum filename length for MARKER NAME field
            mn_len = min(60,length(rinex_metadata.marker_name));
            
            %maximum filename length for OBSERVER field
            ob_len = min(20,length(rinex_metadata.observer));
            
            %maximum filename length for AGENCY field
            oa_len = min(40,length(rinex_metadata.observer_agency));
            
            %system id
            sys = [];
            band1 = [];
            band2 = [];
            attrib = [];
            if (GPSactive), sys = [sys 'G']; band1 = [band1 '1']; attrib = [attrib 'C']; end
            if (GLOactive), sys = [sys 'R']; band1 = [band1 '1']; attrib = [attrib 'C']; end
            if (GALactive), sys = [sys 'E']; band1 = [band1 '1']; attrib = [attrib 'X']; end
            if (BDSactive), sys = [sys 'C']; band1 = [band1 '2']; attrib = [attrib 'I']; end
            if (QZSactive), sys = [sys 'J']; band1 = [band1 '1']; attrib = [attrib 'C']; end
            if (SBSactive), sys = [sys 'S']; band1 = [band1 '1']; attrib = [attrib 'C']; end
            
            if (flag_L2)
                if (GPSactive), band2 = [band2 '2']; end
                if (GLOactive), band2 = [band2 '2']; end
                if (GALactive), band2 = [band2 '2']; end
                if (BDSactive), band2 = [band2 '7']; end
                if (QZSactive), band2 = [band2 '2']; end
                if (SBSactive), band2 = [band2 '2']; end
            end
            
            %find first available epoch for the active constellations
            [~, t] = find(pr1_M(:,first_epoch:last_epoch) ~= 0, 1);
            
            %write header
            fprintf(fid_obs,'%9.2f           OBSERVATION DATA    %-19s RINEX VERSION / TYPE\n', str2num(rinex_metadata.version), file_type);
            fprintf(fid_obs,'%-20s%-20s%-20sPGM / RUN BY / DATE \n', 'goGPS', rinex_metadata.agency(1:rb_len), [date_str_UTC ' ' time_str_UTC ' UTC']);
            fprintf(fid_obs,'%-60sMARKER NAME         \n',rinex_metadata.marker_name(1:mn_len));
            if (rin_ver_id == 3)
                fprintf(fid_obs,'%-20s                                        MARKER TYPE         \n', rinex_metadata.marker_type);
            end
            fprintf(fid_obs,'%-20s%-40sOBSERVER / AGENCY   \n', rinex_metadata.observer(1:ob_len), rinex_metadata.observer_agency(1:oa_len));
            fprintf(fid_obs,'                                                            REC # / TYPE / VERS \n');
            fprintf(fid_obs,'                                                            ANT # / TYPE        \n');
            fprintf(fid_obs,'%14.4f%14.4f%14.4f                  APPROX POSITION XYZ \n', pos_M(1,1), pos_M(2,1), pos_M(3,1));
            fprintf(fid_obs,'        0.0000        0.0000        0.0000                  ANTENNA: DELTA H/E/N\n');
            fprintf(fid_obs,'     1     1                                                WAVELENGTH FACT L1/2\n');
            if (rin_ver_id == 2)
                if (flag_L2)
                    fprintf(fid_obs,['     6    ' code_type '    P2    L1    L2    S1    S2                  # / TYPES OF OBSERV \n']);
                else
                    fprintf(fid_obs,['     3    ' code_type '    L1    S1                                    # / TYPES OF OBSERV \n']);
                end
            else
                if (flag_L2)
                    for s = 1 : length(sys)
                        fprintf(fid_obs,['%-1s    6 ' code_type(1) '%c%c ' code_type(1) '%c%c L%c%c L%c%c S%c%c S%c%c                              SYS / # / OBS TYPES \n'],sys(s),band1(s),attrib(s),band1(s),attrib(s),band1(s),attrib(s),band2(s),attrib(s),band2(s),attrib(s),band2(s),attrib(s));
                    end
                else
                    for s = 1 : length(sys)
                        fprintf(fid_obs,['%-1s    3 ' code_type(1) '%c%c L%c%c S%c%c                                          SYS / # / OBS TYPES \n'],sys(s),band1(s),attrib(s),band1(s),attrib(s),band1(s),attrib(s));
                    end
                end
                fprintf(fid_obs,'DBHZ                                                        SIGNAL STRENGTH UNIT\n');
            end
            fprintf(fid_obs,'%10.3f                                                  INTERVAL            \n', interval);
            fprintf(fid_obs,'%6d%6d%6d%6d%6d%13.7f     GPS         TIME OF FIRST OBS   \n', ...
                date(first_epoch+t-1,1), date(first_epoch+t-1,2), date(first_epoch+t-1,3), date(first_epoch+t-1,4), date(first_epoch+t-1,5), date(first_epoch+t-1,6));
            if (rin_ver_id == 3)
                for s = 1 : length(sys)
                    fprintf(fid_obs,'%c                                                           SYS / PHASE SHIFTS  \n',sys(s));
                end
            end
            fprintf(fid_obs,'                                                            END OF HEADER       \n');
            
            clear cell_master
            
            %-------------------------------------------------------------------------------
            
            %number of records
            N = length(time_M(first_epoch:last_epoch));
            
            %epoch counter for the current DOY
            e = 1;
            
            %cycle through epochs
            for i = first_epoch : last_epoch
                if (nargin == 6)
                    waitbar(e/N,wait_dlg)
                end
                
                sat = find(pr1_M(:,i) ~= 0);
                n = length(sat);
                
                %assign system and PRN code to each satellite
                [sys, prn] = find_sat_system(sat, constellations);
                
                if (e > 1)
                    current_epoch = datenum([date(i,1), date(i,2), date(i,3), date(i,4), date(i,5), date(i,6)]);
                    if (current_epoch == previous_epoch)
                        continue %some binary streams contained the same epoch twice (some kind of bug...)
                    end
                end
                
                previous_epoch = datenum([date(i,1), date(i,2), date(i,3), date(i,4), date(i,5), date(i,6)]);
                
                %if no observations are available, do not write anything
                if (n > 0)
                    
                    %RINEX v2.xx
                    if (rin_ver_id == 2)
                        %epoch record
                        fprintf(fid_obs,' %02d %2d %2d %2d %2d %10.7f  0 %2d', ...
                            two_digit_year(date(i,1)), date(i,2), date(i,3), date(i,4), date(i,5), date(i,6), n);
                        if (n>12)
                            lines = floor(n/12);
                            for l = 1 : lines
                                for j = 1+(l-1)*12 : 12+(l-1)*12
                                    fprintf(fid_obs,'%c%02d',sys(j),prn(j));
                                end
                                if (l ~= lines || rem(n,12) ~= 0)
                                    fprintf(fid_obs,'\n');
                                    fprintf(fid_obs,'%32s','');
                                end
                            end
                            for j = 13+(l-1)*12 : n
                                fprintf(fid_obs,'%c%02d',sys(j),prn(j));
                            end
                        else
                            for j = 1 : n
                                fprintf(fid_obs,'%c%02d',sys(j),prn(j));
                            end
                        end
                        fprintf(fid_obs,'\n');
                        %observation record(s)
                        for j = 1 : n
                            fprintf(fid_obs,'%14.3f %1d',pr1_M(sat(j),i),floor(snr1_M(sat(j),i)/6));
                            if (flag_L2)
                                fprintf(fid_obs,'%14.3f %1d',pr2_M(sat(j),i),floor(snr2_M(sat(j),i)/6));
                            end
                            if (abs(ph1_M(sat(j),i)) > 1e-100)
                                fprintf(fid_obs,'%14.3f %1d',ph1_M(sat(j),i),floor(snr1_M(sat(j),i)/6));
                            else
                                fprintf(fid_obs,'                ');
                            end
                            if (flag_L2)
                                if (abs(ph2_M(sat(j),i)) > 1e-100)
                                    fprintf(fid_obs,'%14.3f %1d',ph2_M(sat(j),i),floor(snr2_M(sat(j),i)/6));
                                else
                                    fprintf(fid_obs,'                ');
                                end
                            end
                            fprintf(fid_obs,'%14.3f %1d',snr1_M(sat(j),i),floor(snr1_M(sat(j),i)/6));
                            if (flag_L2)
                                fprintf(fid_obs,'\n');
                                fprintf(fid_obs,'%14.3f %1d',snr2_M(sat(j),i),floor(snr2_M(sat(j),i)/6));
                            end
                            fprintf(fid_obs,'\n');
                        end
                        %RINEX v3.xx
                    else
                        %epoch record
                        fprintf(fid_obs,'> %4d %2d %2d %2d %2d %10.7f  0 %2d\n', ...
                            date(i,1), date(i,2), date(i,3), date(i,4), date(i,5), date(i,6), n);
                        %observation record(s)
                        for j = 1 : n
                            fprintf(fid_obs,'%c%02d',sys(j),prn(j));
                            fprintf(fid_obs,'%14.3f %1d',pr1_M(sat(j),i),floor(snr1_M(sat(j),i)/6));
                            if (flag_L2)
                                fprintf(fid_obs,'%14.3f %1d',pr2_M(sat(j),i),floor(snr2_M(sat(j),i)/6));
                            end
                            if (abs(ph1_M(sat(j),i)) > 1e-100)
                                fprintf(fid_obs,'%14.3f %1d',ph1_M(sat(j),i),floor(snr1_M(sat(j),i)/6));
                            else
                                fprintf(fid_obs,'                ');
                            end
                            if (flag_L2)
                                if (abs(ph2_M(sat(j),i)) > 1e-100)
                                    fprintf(fid_obs,'%14.3f %1d',ph2_M(sat(j),i),floor(snr2_M(sat(j),i)/6));
                                else
                                    fprintf(fid_obs,'                ');
                                end
                            end
                            fprintf(fid_obs,'%14.3f %1d',snr1_M(sat(j),i),floor(snr1_M(sat(j),i)/6));
                            if (flag_L2)
                                fprintf(fid_obs,'%14.3f %1d',snr2_M(sat(j),i),floor(snr1_M(sat(j),i)/6));
                            end
                            fprintf(fid_obs,'\n');
                        end
                    end
                end
                
                e = e + 1;
            end
            
            %close RINEX observation file
            fclose(fid_obs);
            
            %----------------------------------------------------------------------------------------------
            % RINEX NAVIGATION FILE
            %----------------------------------------------------------------------------------------------
            
            %if ephemerides are available
            if (~isempty(find(Eph_M(1,:,first_epoch:last_epoch) ~= 0, 1)))
                
                N = size(Eph_M, 3);
                
                %displaying
                if (nargin == 5)
                    fprintf(['Writing rover navigation file...\n']);
                end
                
                if (nargin == 6)
                    waitbar(0,wait_dlg,'Writing rover navigation file...')
                end
                
                ext = [];
                if (rin_ver_id == 2)
                    sat_sys = cell(0);
                    m = 1;
                    if (GPSactiveEPH), sat_sys{m} = 'N: GPS NAV DATA';  m = m + 1; ext = [ext 'n']; end
                    if (GLOactiveEPH), sat_sys{m} = 'GLONASS NAV DATA';            ext = [ext 'g']; end
                else
                    if (mixed_sys_eph)
                        ext = 'p';
                    elseif(single_sys_eph && GPSactiveEPH)
                        ext = 'n';
                    elseif(single_sys_eph && GLOactiveEPH)
                        ext = 'g';
                    elseif(single_sys_eph && GALactiveEPH)
                        ext = 'l';
                    end
                    sat_sys{1} = file_type;
                end
                
                for f = 1 : length(ext)
                    
                    %create RINEX navigation file
                    fid_nav = fopen([path marker sprintf('%03d', DOYs(d)) '0' sprintf('.%2d%c', two_digit_year(date(1,1)), ext(f))], 'wt');
                    
                    %write header
                    if (rin_ver_id == 2)
                        fprintf(fid_nav,'%9.2f           %-39s RINEX VERSION / TYPE\n', str2num(rinex_metadata.version), sat_sys{f});
                    else
                        fprintf(fid_nav,'%9.2f           N: GNSS NAV DATA    %-19s RINEX VERSION / TYPE\n', str2num(rinex_metadata.version), sat_sys{f});
                    end
                    fprintf(fid_nav,'%-20s%-20s%-20sPGM / RUN BY / DATE \n', 'goGPS', rinex_metadata.agency(1:rb_len), [date_str_UTC ' ' time_str_UTC ' UTC']);
                    fprintf(fid_nav,'                                                            END OF HEADER       \n');
                    
                    %epoch counter for the current DOY
                    e = 1;
                    
                    for i = first_epoch : last_epoch
                        if (nargin == 6)
                            waitbar(e/N,wait_dlg)
                        end
                        
                        idxEph = find(Eph_M(1,:,i) ~= 0);
                        
                        toe_old = -1;
                        week_E = week_M(first_epoch);
                        
                        for j = 1 : length(idxEph)
                            
                            lineE = [];
                            linesE = [];
                            
                            %if GPS/Galileo/BeiDou/QZSS
                            if ((strcmp(ext(f),'n') &&  strcmp(char(Eph_M(31,idxEph(j),i)),'G')) || ...
                                    (strcmp(ext(f),'p') && (strcmp(char(Eph_M(31,idxEph(j),i)),'G')  || ...
                                    strcmp(char(Eph_M(31,idxEph(j),i)),'E')  || ...
                                    strcmp(char(Eph_M(31,idxEph(j),i)),'C')  || ...
                                    strcmp(char(Eph_M(31,idxEph(j),i)),'J'))))
                                
                                prn      = Eph_M(1,idxEph(j),i);
                                af2      = Eph_M(2,idxEph(j),i);
                                M0       = Eph_M(3,idxEph(j),i);
                                roota    = Eph_M(4,idxEph(j),i);
                                deltan   = Eph_M(5,idxEph(j),i);
                                ecc      = Eph_M(6,idxEph(j),i);
                                omega    = Eph_M(7,idxEph(j),i);
                                cuc      = Eph_M(8,idxEph(j),i);
                                cus      = Eph_M(9,idxEph(j),i);
                                crc      = Eph_M(10,idxEph(j),i);
                                crs      = Eph_M(11,idxEph(j),i);
                                i0       = Eph_M(12,idxEph(j),i);
                                idot     = Eph_M(13,idxEph(j),i);
                                cic      = Eph_M(14,idxEph(j),i);
                                cis      = Eph_M(15,idxEph(j),i);
                                Omega0   = Eph_M(16,idxEph(j),i);
                                Omegadot = Eph_M(17,idxEph(j),i);
                                toe      = Eph_M(18,idxEph(j),i);
                                af0      = Eph_M(19,idxEph(j),i);
                                af1      = Eph_M(20,idxEph(j),i);
                                toc      = Eph_M(21,idxEph(j),i);
                                IODE     = Eph_M(22,idxEph(j),i);
                                codes    = Eph_M(23,idxEph(j),i);
                                weekno   = Eph_M(24,idxEph(j),i); %#ok<NASGU>
                                L2flag   = Eph_M(25,idxEph(j),i);
                                svaccur  = Eph_M(26,idxEph(j),i);
                                svhealth = Eph_M(27,idxEph(j),i);
                                tgd      = Eph_M(28,idxEph(j),i);
                                fit_int  = Eph_M(29,idxEph(j),i);
                                sys      = Eph_M(31,idxEph(j),i);
                                
                                if (toe < toe_old)
                                    week_E = week_E + 1;
                                end
                                
                                toe_old = toe;
                                
                                %time of measurement decoding
                                date_toc = gps2date(week_E, toc);
                                
                                if (rin_ver_id == 2)
                                    
                                    date_toc(1) = two_digit_year(date_toc(1));
                                    
                                    lineE(1,:) = sprintf('%2d %02d %2d %2d %2d %2d%5.1f% 18.12E% 18.12E% 18.12E\n', ...
                                        prn ,date_toc(1), date_toc(2), date_toc(3), date_toc(4), date_toc(5), date_toc(6), af0, af1, af2);
                                    linesE(1,:) = sprintf('   % 18.12E% 18.12E% 18.12E% 18.12E\n', IODE , crs, deltan, M0);
                                    linesE(2,:) = sprintf('   % 18.12E% 18.12E% 18.12E% 18.12E\n', cuc, ecc, cus, roota);
                                    linesE(3,:) = sprintf('   % 18.12E% 18.12E% 18.12E% 18.12E\n', toe, cic, Omega0, cis);
                                    linesE(4,:) = sprintf('   % 18.12E% 18.12E% 18.12E% 18.12E\n', i0, crc, omega, Omegadot);
                                    linesE(5,:) = sprintf('   % 18.12E% 18.12E% 18.12E% 18.12E\n', idot, codes, week, L2flag);
                                    linesE(6,:) = sprintf('   % 18.12E% 18.12E% 18.12E% 18.12E\n', svaccur, svhealth, tgd, IODE);
                                    linesE(7,:) = sprintf('   % 18.12E% 18.12E% 18.12E% 18.12E\n', toc, fit_int, 0, 0);  %here "toc" should be "tom" (e.g. derived from Z-count in Hand Over Word)
                                    %if RINEX v3.xx
                                else
                                    lineE(1,:) = sprintf('%c%02d %04d %2d %2d %2d %2d %2d% 18.12E% 18.12E% 18.12E\n', ...
                                        char(sys), prn, date_toc(1), date_toc(2), date_toc(3), date_toc(4), date_toc(5), date_toc(6), af0, af1, af2);
                                    linesE(1,:) = sprintf('    % 18.12E% 18.12E% 18.12E% 18.12E\n', IODE , crs, deltan, M0);
                                    linesE(2,:) = sprintf('    % 18.12E% 18.12E% 18.12E% 18.12E\n', cuc, ecc, cus, roota);
                                    linesE(3,:) = sprintf('    % 18.12E% 18.12E% 18.12E% 18.12E\n', toe, cic, Omega0, cis);
                                    linesE(4,:) = sprintf('    % 18.12E% 18.12E% 18.12E% 18.12E\n', i0, crc, omega, Omegadot);
                                    linesE(5,:) = sprintf('    % 18.12E% 18.12E% 18.12E% 18.12E\n', idot, codes, week, L2flag);
                                    linesE(6,:) = sprintf('    % 18.12E% 18.12E% 18.12E% 18.12E\n', svaccur, svhealth, tgd, IODE);
                                    linesE(7,:) = sprintf('    % 18.12E% 18.12E% 18.12E% 18.12E\n', toc, fit_int, 0, 0);  %here "toc" should be "tom" (e.g. derived from Z-count in Hand Over Word)
                                end
                            %if GLONASS
                            elseif ((strcmp(ext(f),'g') || strcmp(ext(f),'p')) && strcmp(char(Eph_M(31,idxEph(j),i)),'R'))
                                prn      = Eph_M(1,idxEph(j),i);
                                TauN     = Eph_M(2,idxEph(j),i);
                                GammaN   = Eph_M(3,idxEph(j),i);
                                tk       = Eph_M(4,idxEph(j),i);
                                X        = Eph_M(5,idxEph(j),i)/1e3;  %satellite X coordinate at ephemeris reference time [km]
                                Y        = Eph_M(6,idxEph(j),i)/1e3;  %satellite Y coordinate at ephemeris reference time [km]
                                Z        = Eph_M(7,idxEph(j),i)/1e3;  %satellite Z coordinate at ephemeris reference time [km]
                                Xv       = Eph_M(8,idxEph(j),i)/1e3;  %satellite velocity along X at ephemeris reference time [m/s]
                                Yv       = Eph_M(9,idxEph(j),i)/1e3;  %satellite velocity along Y at ephemeris reference time [m/s]
                                Zv       = Eph_M(10,idxEph(j),i)/1e3; %satellite velocity along Z at ephemeris reference time [m/s]
                                Xa       = Eph_M(11,idxEph(j),i)/1e3; %acceleration due to lunar-solar gravitational perturbation along X at ephemeris reference time [m/s^2]
                                Ya       = Eph_M(12,idxEph(j),i)/1e3; %acceleration due to lunar-solar gravitational perturbation along Y at ephemeris reference time [m/s^2]
                                Za       = Eph_M(13,idxEph(j),i)/1e3; %acceleration due to lunar-solar gravitational perturbation along Z at ephemeris reference time [m/s^2]
                                E        = Eph_M(14,idxEph(j),i);
                                freq_num = Eph_M(15,idxEph(j),i);
                                tb       = Eph_M(16,idxEph(j),i); %#ok<NASGU>
                                Bn       = Eph_M(27,idxEph(j),i);
                                sys      = Eph_M(31,idxEph(j),i);
                                
                                %time of measurement decoding and
                                %conversion from GPS date to GLONASS (UTC) date
                                %date_GPS = gps2date(week, time_M(i));
                                date_GPS = date(i,:);
                                date_GLO = datevec(gps2utc(datenum(date_GPS)));
                                %sec_of_day = date_GLO(3)*3600 + date_GLO(2)*60 + date_GLO(1);
                                datenum_GLO = datenum([date_GLO(1) date_GLO(2) date_GLO(3) 0 0 0]) + datenum(tk/86400);
                                date_GLO = datevec(datenum_GLO);
                                
                                if (rin_ver_id == 2)
                                    
                                    date_GLO(1) = two_digit_year(date_GLO(1));
                                    
                                    lineE(1,:) = sprintf('%2d %02d %2d %2d %2d %2d%5.1f% 18.12E% 18.12E% 18.12E\n', ...
                                        prn, date_GLO(1), date_GLO(2), date_GLO(3), date_GLO(4), date_GLO(5), date_GLO(6), -TauN, GammaN, tk);
                                    linesE(1,:) = sprintf('   % 18.12E% 18.12E% 18.12E% 18.12E\n', X, Xv, single(Xa), Bn);
                                    linesE(2,:) = sprintf('   % 18.12E% 18.12E% 18.12E% 18.12E\n', Y, Yv, single(Ya), freq_num);
                                    linesE(3,:) = sprintf('   % 18.12E% 18.12E% 18.12E% 18.12E\n', Z, Zv, single(Za), E);
                                    %if RINEX v3.xx
                                else
                                    lineE(1,:) = sprintf('%c%02d %04d %2d %2d %2d %2d %2d% 18.12E% 18.12E% 18.12E\n', ...
                                        char(sys), prn, date_GLO(1), date_GLO(2), date_GLO(3), date_GLO(4), date_GLO(5), date_GLO(6), -TauN, GammaN, tk);
                                    linesE(1,:) = sprintf('    % 18.12E% 18.12E% 18.12E% 18.12E\n', X, Xv, single(Xa), Bn);
                                    linesE(2,:) = sprintf('    % 18.12E% 18.12E% 18.12E% 18.12E\n', Y, Yv, single(Ya), freq_num);
                                    linesE(3,:) = sprintf('    % 18.12E% 18.12E% 18.12E% 18.12E\n', Z, Zv, single(Za), E);
                                end
                            end
                            
                            if (~isempty(lineE))
                                %if running on Windows, convert three-digits exponential notation
                                %to two-digits; in any case, replace 'E' with 'D' and print the string
                                n = size(linesE,1);
                                if (~isunix)
                                    lineD = strrep(lineE(1,:),'E+0','D+');
                                    lineD = strrep(lineD,'E-0','D-');
                                    fprintf(fid_nav,'%s',lineD);
                                    for k = 1 : n
                                        lineD = strrep(linesE(k,:),'E+0','D+');
                                        lineD = strrep(lineD,'E-0','D-');
                                        fprintf(fid_nav,'%s',lineD);
                                    end
                                else
                                    lineD = strrep(lineE(1,:),'E+','D+');
                                    lineD = strrep(lineD,'E-','D-');
                                    fprintf(fid_nav,'%s',lineD);
                                    for k = 1 : n
                                        lineD = strrep(linesE(k,:),'E+','D+');
                                        lineD = strrep(lineD,'E-','D-');
                                        fprintf(fid_nav,'%s',lineD);
                                    end
                                end
                            end
                        end
                        
                        e = e + 1;
                    end
                    
                    %close RINEX navigation file
                    fclose(fid_nav);
                end
            end
        end
    else
        %displaying
        if (nargin == 6)
            msgbox('Master: no navigation data acquired.');
        else
            fprintf('Master: no navigation data acquired! \n');
        end
    end
else
    %displaying
    if (nargin == 6)
        msgbox('No master data acquired.');
    else
        fprintf('No master data acquired! \n');
    end
end
