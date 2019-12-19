function [pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2, ...
          time_zero, time_GPS, time, week, date, pos, interval, antoff, antmod, codeC1, marker] = ...
          load_RINEX_obs(filename, cc, processing_interval, wait_dlg)

% SYNTAX:
%   [pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2, ...
%    time_GPS, time, week, date, pos, interval, antoff, antmod, codeC1] = ...
%    load_RINEX_obs(filename, cc, processing_interval, wait_dlg);
%
% INPUT:
%   filename = RINEX observation file(s)
%   cc = Constellation_Collector object, contains the satus of the satellite systems in use
%   processing_interval = user-requested processing interval
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
%   codeC1 = boolean variable to notify if the C1 code is used instead of P1
%   marker = marker name [string]
%
% DESCRIPTION:
%   Parses RINEX observation files.

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
%  Contributors:     Damianop Triglione,
%                    Stefano Caldera,
%                    Andrea Gatti, ...
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

global report

log = Logger.getInstance();
state = Core.getCurrentSettings();

% Check the input arguments
if (nargin < 4)
    wait_dlg_PresenceFlag = false;
else
    wait_dlg_PresenceFlag = true;
end

if (nargin < 3 || processing_interval < 0)
    processing_interval = 0;
end

%number of satellite slots for enabled constellations
nSatTot = cc.getNumSat();

%number of RINEX files to be read
if (iscell(filename))
    nFiles = size(filename,1);
else
    nFiles = 1;
end

%variable initialization
nEpochs = 87000;
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
codeC1 = zeros(nSatTot,nEpochs,nFiles);
marker = cell(1,1,nFiles);

max_k = 0;
for f = 1 : nFiles
    if (iscell(filename))
        current_file = filename{f,1};
    else
        current_file = filename;
    end

    log.addMessage(sprintf('Reading file %s', current_file));
    File_Rinex(current_file,9);

    %open RINEX observation file
    fid = fopen(current_file,'r');

    if (wait_dlg_PresenceFlag)
        waitbar(0.5,wait_dlg, ['RINEX file ' current_file ': parsing header...'])
    end

    %parse RINEX header
    [obs_type, pos(:,1,f), basic_info, interval(1,1,f), sysId, antoff(:,1,f), antmod{1,1,f}, marker{1,1,f}] = RINEX_parse_hdr_old(fid);

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
        [date(k,:,f), num_sat, sat, sat_types] = RINEX_get_epoch_old(fid);
        if ~isempty(date(k,1,f))
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
            obs = RINEX_get_obs_old(fid, num_sat, sat, sat_types, obsColumns, nObsTypes, cc);

            idx_P1 = obs.P1 ~= 0;
            idx_C1 = obs.C1 ~= 0;
            idx_codeC1 = idx_P1 - idx_C1;
            codeC1(idx_codeC1 < 0,k,f) = 1;
            pr1(:,k,f) = zeros(size(pr1(:,k,f)));
            pr1(idx_P1,k,f) = obs.P1(idx_P1);
            pr1(find(codeC1(:,k,f)),k,f) = obs.C1(find(codeC1(:,k,f))); %#ok<FNDSB>

            %         %read ROVER observations
            %         if (~any(obs.C1) || sum(obs.P1 ~= 0) == sum(obs.C1 ~= 0))
            %             pr1(:,k,f) = obs.P1;
            %         else
            %             pr1(:,k,f) = obs.C1;
            %             codeC1(:,:,f) = 1;
            %         end
            pr2(:,k,f) = obs.P2;
            ph1(:,k,f) = obs.L1;
            ph2(:,k,f) = obs.L2;
            dop1(:,k,f) = obs.D1;
            dop2(:,k,f) = obs.D2;
            snr1(:,k,f) = obs.S1;
            snr2(:,k,f) = obs.S2;
            k = k + 1;
        end
    end

    max_k = max(max_k, k-1);

    if (wait_dlg_PresenceFlag)
        waitbar(1,wait_dlg)
    end

    time(f) = GPS_Time(date(1:k-1,:,f));
    % try to guess observation rate when not read from header
    if (interval(1,1,f) == 0)
        interval(1,1,f) = time(f).getRate();
    end

    %-------------------------------------------------------------------------------

    fclose(fid);

%     if (processing_interval > interval(:,1,f))
%         interval(:,1,f) = processing_interval;
%     end

    if (~isempty(report) && report.opt.write == 1)
        % extract quality parameters for report
        j=strfind(current_file,'\');
        if isempty(j)
            j=strfind(current_file,'/');
        end
        if isempty(j)
            j=0;
        end
        report.obs.filename(f)=cellstr(current_file(j(end)+1:end));
        % create statistics on observations
        stat_sat = ((ph1(:,:,f)~=0 & isfinite(ph1(:,:,f))) + (ph2(:,:,f)~=0 & isfinite(ph2(:,:,f))) + ...
            (pr1(:,:,f)~=0 & isfinite(pr1(:,:,f))) + (pr2(:,:,f)~=0 & isfinite(pr2(:,:,f))) + ...
            (dop1(:,:,f)~=0 & isfinite(dop1(:,:,f))) + (dop2(:,:,f)~=0 & isfinite(dop2(:,:,f))))~=0;
        report.obs_raw.n_sat(f)=sum(sum(stat_sat,2)~=0);
        report.obs_raw.n_epoch(f)=sum(sum(stat_sat,1)~=0);
        report.obs_raw.n_ph1(f) = sum(sum((ph1(:,:,f)~=0 & isfinite(ph1(:,:,f)))));
        report.obs_raw.n_ph2(f) = sum(sum((ph2(:,:,f)~=0 & isfinite(ph2(:,:,f)))));
        report.obs_raw.n_pr1(f) = sum(sum((pr1(:,:,f)~=0 & isfinite(pr1(:,:,f)))));
        report.obs_raw.n_pr2(f) = sum(sum((pr2(:,:,f)~=0 & isfinite(pr2(:,:,f)))));
        report.obs_raw.n_dop1(f) = sum(sum((dop1(:,:,f)~=0 & isfinite(dop1(:,:,f)))));
        report.obs_raw.n_dop2(f) = sum(sum((dop2(:,:,f)~=0 & isfinite(dop2(:,:,f)))));
        report.obs_raw.interval(f) = interval(1,1,f);
        report.obs_raw.time_start(f)=cellstr(sprintf('%04d-%02d-%02d %02d:%02d:%06.3f',date(1,1,f),date(1,2,f),date(1,3,f),date(1,4,f),date(1,5,f),date(1,6,f)));
        report.obs_raw.time_end(f)=cellstr(sprintf('%04d-%02d-%02d %02d:%02d:%06.3f',date(k-1,1,f),date(k-1,2,f),date(k-1,3,f),date(k-1,4,f),date(k-1,5,f),date(k-1,6,f)));

        stat_sat = ((ph2(:,:,f)~=0 & isfinite(ph2(:,:,f))) + ((pr2(:,:,f)~=0 & isfinite(pr2(:,:,f)))) + (dop2(:,:,f)~=0 & isfinite(dop2(:,:,f))))~=0;
        if any(stat_sat(:))
            report.obs_raw.nfreq(f)=2;
        else
            report.obs_raw.nfreq(f)=1;
        end
        report.obs_raw.n_epoch_expected(f) = time(f).getExpectedLen();

        report.obs_raw.epoch_completeness(f)=report.obs_raw.n_epoch(f)/report.obs_raw.n_epoch_expected(f)*100;
        if report.obs_raw.nfreq(f) == 2
            report.obs_raw.L1L2_completeness(f) = cellstr(sprintf('%6.1f', report.obs_raw.n_ph2(f)/report.obs_raw.n_ph1(f)*100));
        else
            report.obs_raw.L1L2_completeness(f) = cellstr(sprintf('%6s', '0'));
        end
    end

    log.addMessage('Done reading current RINEX');
end

% trim output (it have been pre-allocated bigger)
pr1 = pr1(:,(1 : max_k),:);
pr2 = pr2(:,(1 : max_k),:);
ph1 = ph1(:,(1 : max_k),:);
ph2 = ph2(:,(1 : max_k),:);
dop1 = dop1(:,(1 : max_k),:);
dop2 = dop2(:,(1 : max_k),:);
snr1 = snr1(:,(1 : max_k),:);
snr2 = snr2(:,(1 : max_k),:);
date = date((1 : max_k),:,:);
codeC1 = codeC1(:,(1 : max_k),:);

log.newLine();
log.addMessage('Syncing observations if needed');

%sync observations
[time_zero, time_GPS, time, week, date, pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2, codeC1, interval] = ...
sync_obs(time, date, pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2, codeC1, interval, processing_interval);

log.addMessage('Trimming short arcs if needed');
% remove short arcs
min_arc = state.getMinArc();
for f = 1 : nFiles
    pr1(:,:,f) = remove_short_arcs(pr1(:,:,f), min_arc);
    pr2(:,:,f) = remove_short_arcs(pr2(:,:,f), min_arc);
    ph1(:,:,f) = remove_short_arcs(ph1(:,:,f), min_arc);
    ph2(:,:,f) = remove_short_arcs(ph2(:,:,f), min_arc);
end

% Find when all the dataset have at least one good observations in common
resync_flag_ok = any(pr1,1) & any(pr1,1) & any(pr1,1) & any(ph1,1);
max_sync = find(sum(resync_flag_ok,3) == size(resync_flag_ok,3), 1, 'last');

pr1 = pr1(:,(1 : max_sync),:);
pr2 = pr2(:,(1 : max_sync),:);
ph1 = ph1(:,(1 : max_sync),:);
ph2 = ph2(:,(1 : max_sync),:);
dop1 = dop1(:,(1 : max_sync),:);
dop2 = dop2(:,(1 : max_sync),:);
snr1 = snr1(:,(1 : max_sync),:);
snr2 = snr2(:,(1 : max_sync),:);
date = date((1 : max_sync),:,:);
codeC1 = codeC1(:,(1 : max_sync),:);
time_GPS = time_GPS(1 : max_sync);
time = time((1 : max_sync),:,:);
week = week((1 : max_sync),:,:);
date = date((1 : max_sync),:,:);

for f = 1 : nFiles
    holes = find(week(:,1,f) == 0);
    for h = holes'
        if (h > 1)
            time(h,:,f) = time(h-1,1,f) + interval;
            week(h,1,f) = week(h-1,1,f);
            date(h,:,f) = datevec(datenum(date(h-1,:,f)) + datenum([0 0 0 0 0 interval]));
        elseif (holes(end)+1 <= length(week(:,1,f)))
            time(h,1,f) = time(holes(end)+1,1,f) - interval*holes(end);
            week(h,1,f) = week(holes(end)+1,1,f);
            date(h,:,f) = datevec(datenum(date(holes(end)+1,:,f)) - datenum([0 0 0 0 0 interval*holes(end)]));
        end
    end
end
log.newLine();

if (~isempty(report) && report.opt.write == 1)
    % extract quality parameters for report
    for f = 1 : nFiles
         % create statistics on observations
        stat_sat = ((ph1(:,:,f)~=0 & isfinite(ph1(:,:,f))) + (ph2(:,:,f)~=0 & isfinite(ph2(:,:,f))) + ...
            (pr1(:,:,f)~=0 & isfinite(pr1(:,:,f))) + (pr2(:,:,f)~=0 & isfinite(pr2(:,:,f))) + ...
            (dop1(:,:,f)~=0 & isfinite(dop1(:,:,f))) + (dop2(:,:,f)~=0 & isfinite(dop2(:,:,f))))~=0;
        report.obs_sync.n_sat(f)=sum(sum(stat_sat,2)~=0);
        report.obs_sync.n_epoch(f)=sum(sum(stat_sat,1)~=0);
        report.obs_sync.n_ph1(f) = sum(sum((ph1(:,:,f)~=0 & isfinite(ph1(:,:,f)))));
        report.obs_sync.n_ph2(f) = sum(sum((ph2(:,:,f)~=0 & isfinite(ph2(:,:,f)))));
        report.obs_sync.n_pr1(f) = sum(sum((pr1(:,:,f)~=0 & isfinite(pr1(:,:,f)))));
        report.obs_sync.n_pr2(f) = sum(sum((pr2(:,:,f)~=0 & isfinite(pr2(:,:,f)))));
        report.obs_sync.n_dop1(f) = sum(sum((dop1(:,:,f)~=0 & isfinite(dop1(:,:,f)))));
        report.obs_sync.n_dop2(f) = sum(sum((dop2(:,:,f)~=0 & isfinite(dop2(:,:,f)))));
        report.obs_sync.interval(f) = interval;
        report.obs_sync.time_start(f)=cellstr(sprintf('%04d-%02d-%02d %02d:%02d:%06.3f',date(1,1,f),date(1,2,f),date(1,3,f),date(1,4,f),date(1,5,f),date(1,6,f)));
        report.obs_sync.time_end(f)=cellstr(sprintf('%04d-%02d-%02d %02d:%02d:%06.3f',date(size(date,1),1,f),date(size(date,1),2,f),date(size(date,1),3,f),date(size(date,1),4,f),date(size(date,1),5,f),date(size(date,1),6,f)));

        stat_sat = ((ph2(:,:,f)~=0 & isfinite(ph2(:,:,f))) + ((pr2(:,:,f)~=0 & isfinite(pr2(:,:,f)))) + (dop2(:,:,f)~=0 & isfinite(dop2(:,:,f))))~=0;
        if any(stat_sat(:))
            report.obs_sync.nfreq(f)=2;
        else
            report.obs_sync.nfreq(f)=1;
        end
        report.obs_sync.n_epoch_expected(f) = length((roundmod(time(1,1,f),interval) : interval : roundmod(time(size(date,1),1,f),interval)));

        report.obs_sync.epoch_completeness(f)=report.obs_sync.n_epoch(f)/report.obs_sync.n_epoch_expected(f)*100;
        if report.obs_sync.nfreq(f) == 2
            report.obs_sync.L1L2_completeness(f) = cellstr(sprintf('%6.1f', report.obs_sync.n_ph2(f)/report.obs_sync.n_ph1(f)*100));
        else
            report.obs_sync.L1L2_completeness(f) = cellstr(sprintf('%6s', '0'));
        end
    end
end
