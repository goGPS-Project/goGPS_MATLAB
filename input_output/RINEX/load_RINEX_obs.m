function [pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2, ...
          time_ref, time, week, date, pos, interval, antoff, antmod] = ...
          load_RINEX_obs(filename, constellations, wait_dlg)

% SYNTAX:
%   [pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2, ...
%    time_ref, time, week, date, pos, interval, antoff, antmod] = ...
%    load_RINEX_obs(filename, constellations, wait_dlg);
%
% INPUT:
%   filename = RINEX observation file(s)
%   constellations = struct with multi-constellation settings
%                   (see 'multi_constellation_settings.m' - empty if not available)
%   wait_dlg = optional handler to waitbar figure (optional)
%
% OUTPUT:
%   pr1 = code observation (L1 carrier)
%   ph1 = phase observation (L1 carrier)
%   pr2 = code observation (L2 carrier)
%   ph2 = phase observation (L2 carrier)
%   dop1 = Doppler observation (L1 carrier)
%   dop2 = Doppler observation (L2 carrier)
%   snr1 = signal-to-noise ratio (L1 carrier)
%   snr2 = signal-to-noise ratio (L2 carrier)
%   time = receiver seconds-of-week
%   week = GPS week
%   date = date (year,month,day,hour,minute,second)
%   pos = rover approximate position
%   interval = observation time interval [s]
%   antoff = antenna offset [m]
%   antmod = antenna model [string]
%
% DESCRIPTION:
%   Parses RINEX observation files.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.2 beta
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

global fout_report

% Check the input arguments
if (nargin < 3)
    wait_dlg_PresenceFlag = false;
else
    wait_dlg_PresenceFlag = true;
end

if (nargin < 2 || isempty(constellations)) %then use only GPS as default
    [constellations] = goGNSS.initConstellation(1, 0, 0, 0, 0, 0);
end

%number of satellite slots for enabled constellations
nSatTot = constellations.nEnabledSat;

%number of RINEX files to be read
if (iscell(filename))
    nFiles = size(filename,1);
else
    nFiles = 1;
end

%variable initialization
nEpochs = 86400;
time = NaN(nEpochs,1,nFiles);
tow = NaN(nEpochs,1,nFiles);
week = NaN(nEpochs,1,nFiles);
pr1 = NaN(nSatTot,nEpochs,nFiles);
pr2 = NaN(nSatTot,nEpochs,nFiles);
ph1 = NaN(nSatTot,nEpochs,nFiles);
ph2 = NaN(nSatTot,nEpochs,nFiles);
dop1 = NaN(nSatTot,nEpochs,nFiles);
dop2 = NaN(nSatTot,nEpochs,nFiles);
snr1 = NaN(nSatTot,nEpochs,nFiles);
snr2 = NaN(nSatTot,nEpochs,nFiles);
date = NaN(nEpochs,6,nFiles);
pos = zeros(3,1,nFiles);
interval = zeros(1,1,nFiles);
antoff = zeros(3,1,nFiles);
antmod = cell(1,1,nFiles);

for f = 1 : nFiles

    if (iscell(filename))
        current_file = filename{f,1};
    else
        current_file = filename;
    end
    
    fprintf('%s',['Reading RINEX file ' current_file ': ... ']);
    
    %open RINEX observation file
    fid = fopen(current_file,'r');
    
    if (wait_dlg_PresenceFlag)
        waitbar(0.5,wait_dlg,['RINEX file ' current_file ': parsing header...'])
    end
    
    %parse RINEX header
    [obs_type, pos(:,1,f), basic_info, interval(1,1,f), sysId, antoff(:,1,f), antmod{1,1,f}] = RINEX_parse_hdr(fid);
    
    %check the availability of basic data to parse the RINEX file
    if (basic_info == 0)
        error(['RINEX file ' current_file ': basic data is missing in the file header'])
    end
    
    %find observation type columns
    [obsColumns, nObsTypes] = obs_type_find(obs_type, sysId);
    
    if (wait_dlg_PresenceFlag)
        waitbar(1,wait_dlg)
    end
    
    %-------------------------------------------------------------------------------
    
    if (wait_dlg_PresenceFlag)
        waitbar(0.5,wait_dlg,['RINEX file ' current_file ': reading observations...'])
    end
    
    k = 1;
    while (~feof(fid))
        
        %read data for the current epoch (ROVER)
        [time(k,1,f), date(k,:,f), num_sat, sat, sat_types, tow(k,1,f)] = RINEX_get_epoch(fid);
        
        if (k > nEpochs)
            %variable initialization (GPS)
            pr1(:,k,f) = zeros(nSatTot,1);
            pr2(:,k,f) = zeros(nSatTot,1);
            ph1(:,k,f) = zeros(nSatTot,1);
            ph2(:,k,f) = zeros(nSatTot,1);
            dop1(:,k,f) = zeros(nSatTot,1);
            dop2(:,k,f) = zeros(nSatTot,1);
            snr1(:,k,f) = zeros(nSatTot,1);
            snr2(:,k,f) = zeros(nSatTot,1);
            
            nEpochs = nEpochs  + 1;
        end
        
        %read ROVER observations
        obs = RINEX_get_obs(fid, num_sat, sat, sat_types, obsColumns, nObsTypes, constellations);
        
        %read ROVER observations
        if (~any(obs.C1) || sum(obs.P1 ~= 0) == sum(obs.C1 ~= 0))
            pr1(:,k,f) = obs.P1;
        else
            pr1(:,k,f) = obs.C1;
        end
        pr2(:,k,f) = obs.P2;
        ph1(:,k,f) = obs.L1;
        ph2(:,k,f) = obs.L2;
        dop1(:,k,f) = obs.D1;
        dop2(:,k,f) = obs.D2;
        snr1(:,k,f) = obs.S1;
        snr2(:,k,f) = obs.S2;
        
        k = k + 1;
    end
    
    if (wait_dlg_PresenceFlag)
        waitbar(1,wait_dlg)
    end
    
    %GPS week number
    week(:,1,f) = date2gps(date(:,:,f));
    
    %observation rate
    if (interval(:,1,f) == 0)
        interval(:,1,f) = median(time(2:k-1,1,f) - time(1:k-2,1,f));
    end
    
    %-------------------------------------------------------------------------------
    
    %close RINEX files
    fclose(fid);

    if exist('fout_report','var') &&  fout_report~=-1
        % write report file
        if f==1
            % table header
            fprintf(fout_report,'\n\nSTATION INFORMATION\n');
            fprintf(fout_report,'-------------------\n\n');
            fprintf(fout_report,'Observations (RAW)              Start time               End time                 Rate  #Sat   #Epoch    #Frq   #C1/P1  #C2/P2     #L1     #L2   #DOP1   #DOP2  %%Epoch %%L2/L1\n');
            fprintf(fout_report,'-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
        end
        j=strfind(current_file,'\');
        if isempty(j)
            j=strfind(current_file,'/');
        end
        if isempty(j)
            j=0;
        end
        % create statistics on observations
        stat_sat = ((ph1(:,:,f)~=0 & isfinite(ph1(:,:,f))) + (ph2(:,:,f)~=0 & isfinite(ph2(:,:,f))) + ...
            (pr1(:,:,f)~=0 & isfinite(pr1(:,:,f))) + (pr2(:,:,f)~=0 & isfinite(pr2(:,:,f))) + ...
            (dop1(:,:,f)~=0 & isfinite(dop1(:,:,f))) + (dop2(:,:,f)~=0 & isfinite(dop2(:,:,f))))~=0;
        n_sat=sum(sum(stat_sat,2)~=0);
        n_epoch=sum(sum(stat_sat,1)~=0);
        n_ph1 = sum(sum((ph1(:,:,f)~=0 & isfinite(ph1(:,:,f)))));
        n_ph2 = sum(sum((ph2(:,:,f)~=0 & isfinite(ph2(:,:,f)))));
        n_pr1 = sum(sum((pr1(:,:,f)~=0 & isfinite(pr1(:,:,f)))));
        n_pr2 = sum(sum((pr2(:,:,f)~=0 & isfinite(pr2(:,:,f)))));
        n_dop1 = sum(sum((dop1(:,:,f)~=0 & isfinite(dop1(:,:,f)))));
        n_dop2 = sum(sum((dop2(:,:,f)~=0 & isfinite(dop2(:,:,f)))));
        time_start=sprintf('%04d-%02d-%02d %02d:%02d:%06.3f',date(1,1,f),date(1,2,f),date(1,3,f),date(1,4,f),date(1,5,f),date(1,6,f));
        time_end=sprintf('%04d-%02d-%02d %02d:%02d:%06.3f',date(k-1,1,f),date(k-1,2,f),date(k-1,3,f),date(k-1,4,f),date(k-1,5,f),date(k-1,6,f));
        
        stat_sat = ((ph2(:,:,f)~=0 & isfinite(ph2(:,:,f))) + ((pr2(:,:,f)~=0 & isfinite(pr2(:,:,f)))) + (dop2(:,:,f)~=0 & isfinite(dop2(:,:,f))))~=0;
        if any(stat_sat(:))
            nfreq=2;
        else
            nfreq=1;
        end
		
        n_epoch_expected = length((roundmod(time(1,1,f),interval(1,1,f)) : interval(1,1,f) : roundmod(time(k-1,1,f),interval(1,1,f))));

        epoch_completeness=n_epoch/n_epoch_expected*100;
        if nfreq == 2
            L1L2_completeness = sprintf('%6.1f', n_ph2/n_ph1*100);
        else
            L1L2_completeness = sprintf('%6s', '0');
        end
        
        fprintf(fout_report,'%-30s  %23s  %23s  %4.1f    %2d   %6d  %6d  %7d %7d %7d %7d %7d %7d  %6.1f %-6s\n', current_file(j(end)+1:end), time_start, time_end, interval(1,1,f), ...
            n_sat, n_epoch, nfreq, n_pr1, n_pr2, n_ph1, n_ph2, n_dop1, n_dop2, epoch_completeness, L1L2_completeness);
    end
    
    
    fprintf('done\n');
end

%sync observations
[time_ref, time, week, date, pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2, interval] = ...
sync_obs(time, week, date, pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2, interval);

for f = 1 : nFiles
    holes = find(week(:,1,f) == 0);
    for h = holes'
        if (h > 1)
            week(h,1,f) = week(h-1,1,f);
            date(h,:,f) = datevec(datenum(date(h-1,:,f)) + datenum([0 0 0 0 0 interval]));
        end
    end
end


if exist('fout_report','var') &&  fout_report~=-1
    % write report file
    fprintf(fout_report,'\nObservations (after sync)       Start time               End time                 Rate  #Sat   #Epoch    #Frq   #C1/P1  #C2/P2     #L1     #L2   #DOP1   #DOP2  %%Epoch %%L2/L1\n');
    fprintf(fout_report,'-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
    
    for f = 1 : nFiles        
        if (iscell(filename))
            current_file = filename{f,1};
        else
            current_file = filename;
        end
        
                j=strfind(current_file,'\');
        if isempty(j)
            j=strfind(current_file,'/');
        end
        if isempty(j)
            j=0;
        end
        % create statistics on observations
        stat_sat = ((ph1(:,:,f)~=0 & isfinite(ph1(:,:,f))) + (ph2(:,:,f)~=0 & isfinite(ph2(:,:,f))) + ...
            (pr1(:,:,f)~=0 & isfinite(pr1(:,:,f))) + (pr2(:,:,f)~=0 & isfinite(pr2(:,:,f))) + ...
            (dop1(:,:,f)~=0 & isfinite(dop1(:,:,f))) + (dop2(:,:,f)~=0 & isfinite(dop2(:,:,f))))~=0;
        n_sat=sum(sum(stat_sat,2)~=0);
        n_epoch=sum(sum(stat_sat,1)~=0);
        n_ph1 = sum(sum((ph1(:,:,f)~=0 & isfinite(ph1(:,:,f)))));
        n_ph2 = sum(sum((ph2(:,:,f)~=0 & isfinite(ph2(:,:,f)))));
        n_pr1 = sum(sum((pr1(:,:,f)~=0 & isfinite(pr1(:,:,f)))));
        n_pr2 = sum(sum((pr2(:,:,f)~=0 & isfinite(pr2(:,:,f)))));
        n_dop1 = sum(sum((dop1(:,:,f)~=0 & isfinite(dop1(:,:,f))))); 
        n_dop2 = sum(sum((dop2(:,:,f)~=0 & isfinite(dop2(:,:,f)))));
        time_start=sprintf('%04d-%02d-%02d %02d:%02d:%06.3f',date(1,1,f),date(1,2,f),date(1,3,f),date(1,4,f),date(1,5,f),date(1,6,f));
        time_end=sprintf('%04d-%02d-%02d %02d:%02d:%06.3f',date(end,1,f),date(end,2,f),date(end,3,f),date(end,4,f),date(end,5,f),date(end,6,f));
        
        stat_sat = ((ph2(:,:,f)~=0 & isfinite(ph2(:,:,f))) + ((pr2(:,:,f)~=0 & isfinite(pr2(:,:,f)))) + (dop2(:,:,f)~=0 & isfinite(dop2(:,:,f))))~=0;
        if any(stat_sat(:))
            nfreq=2;
        else
            nfreq=1;
        end
        
        n_epoch_expected = size(time,1);

        epoch_completeness=n_epoch/n_epoch_expected*100;
        if nfreq == 2
            L1L2_completeness = sprintf('%6.1f', n_ph2/n_ph1*100);
        else
            L1L2_completeness = sprintf('%6s', '0');
        end
        
        fprintf(fout_report,'%-30s  %23s  %23s  %4.1f    %2d   %6d  %6d  %7d %7d %7d %7d %7d %7d  %6.1f %-6s\n', current_file(j(end)+1:end), time_start, time_end, interval, ...
            n_sat, n_epoch, nfreq, n_pr1, n_pr2, n_ph1, n_ph2, n_dop1, n_dop2, epoch_completeness, L1L2_completeness);
    end
end
        
