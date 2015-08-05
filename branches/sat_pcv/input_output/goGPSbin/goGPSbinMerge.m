function goGPSbinMerge(filerootR, filerootM, filerootOUT, wait_dlg)

% SYNTAX:
%   goGPSbinMerge(filerootR, filerootM, filerootOUT, wait_dlg);
%
% INPUT:
%   filerootR  = input rover file root
%   filerootM  = input master file root
%   filerootOUT = output file root
%   wait_dlg = optional handler to waitbar figure
%
% OUTPUT:
%
% DESCRIPTION:
%   Merge two distinct observations datasets saved by goGPS (*_obs_* and
%   *_eph_* files).

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

if ~isempty(dir([filerootR '_obs_*'])) & ~isempty(dir([filerootM '_obs_*'])) ...
    & (~isempty(dir([filerootR '_eph_*'])) | ~isempty(dir([filerootM '_eph_*'])))
    
    %ROVER and MASTER datasets reading
    if (nargin == 4)
        [time_GPS, week_R, time_R, time_M, pr1_R, pr1_M, ph1_R, ph1_M, dop1_R, snr_R, snr_M, pos_M, Eph, ...
            iono, loss_R, loss_M] = load_observ(filerootR, filerootM, wait_dlg);
    else
        [time_GPS, week_R, time_R, time_M, pr1_R, pr1_M, ph1_R, ph1_M, dop1_R, snr_R, snr_M, pos_M, Eph, ...
            iono, loss_R, loss_M] = load_observ(filerootR, filerootM);
    end
    
    num_sat = size(pr1_R,1);

    EphAvailable = find(Eph(30,:,:)~=0, 1);
    %if the dataset has ephemerides available at least for one epoch
    if (~isempty(EphAvailable))
        satEph = find(sum(abs(Eph(:,:,1)))~=0);
        satObs = find( (pr1_R(:,1) ~= 0) & (pr1_M(:,1) ~= 0));
        while (length(satEph) < length(satObs)) | (length(satObs) < 4)
            
            time_GPS(1) = [];
            week_R(1)   = [];
            time_R(1)   = [];
            time_M(1)   = [];
            pr1_R(:,1)  = [];
            pr1_M(:,1)  = [];
            ph1_R(:,1)  = [];
            ph1_M(:,1)  = [];
            dop1_R(:,1) = [];
            snr_R(:,1)  = [];
            snr_M(:,1)  = [];
            pos_M(:,1)  = [];
            Eph(:,:,1)  = [];
            iono(:,1)   = [];
            loss_R(1)   = [];
            loss_M(1)   = [];
            
            %satObs_R = find( (pr1_R(:,1) ~= 0) & (ph1_R(:,1) ~= 0) );
            %satObs_M = find( (pr1_M(:,1) ~= 0) & (ph1_M(:,1) ~= 0) );
            satObs = find( (pr1_R(:,1) ~= 0) & (pr1_M(:,1) ~= 0));
            satEph = find(sum(abs(Eph(:,:,1)))~=0);
        end
        
        %remove observations without ephemerides
        for i = 1 : length(time_GPS)
            satEph = find(sum(abs(Eph(:,:,i)))~=0);
            delsat = setdiff(1:num_sat,satEph);
            pr1_R(delsat,i)  = 0;
            pr1_M(delsat,i)  = 0;
            ph1_R(delsat,i)  = 0;
            ph1_M(delsat,i)  = 0;
            dop1_R(delsat,i) = 0;
            snr_R(delsat,i)  = 0;
            snr_M(delsat,i)  = 0;
        end
    else
        if (nargin == 4)
            msgbox('Ephemerides are required to create goGPS binary data.');
        else
            fprintf('Ephemerides are required to create goGPS binary data.\n');
        end
        return
    end

    %complete/partial path
    tMin = 1;
    tMax = 1e30;
    tMin = max(tMin,1);
    tMax = min(tMax,length(time_GPS));
    time_GPS = time_GPS(tMin:tMax);
    week_R = week_R(tMin:tMax);
    time_R = time_R(tMin:tMax);
    time_M = time_M(tMin:tMax);
    pr1_R = pr1_R(:,tMin:tMax);
    pr1_M = pr1_M(:,tMin:tMax);
    ph1_R = ph1_R(:,tMin:tMax);
    ph1_M = ph1_M(:,tMin:tMax);
    dop1_R = dop1_R(:,tMin:tMax);
    snr_R = snr_R(:,tMin:tMax);
    snr_M = snr_M(:,tMin:tMax);
    pos_M = pos_M(:,tMin:tMax);
    Eph = Eph(:,:,tMin:tMax);
    iono = iono(:,tMin:tMax);
    
    %do not overwrite existing files
    i = 1;
    j = length(filerootOUT);
    while (~isempty(dir([filerootOUT '_obs*.bin'])) | ...
            ~isempty(dir([filerootOUT '_eph*.bin'])) )
        
        filerootOUT(j+1:j+4) = ['_' num2str(i,'%03d')];
        i = i + 1;
    end
    
    %open output files
    fid_obs = fopen([filerootOUT '_obs_000.bin'],'w+');
    fid_eph = fopen([filerootOUT '_eph_000.bin'],'w+');
    
    %write number of satellites
    fwrite(fid_obs, num_sat, 'int8');
    fwrite(fid_eph, num_sat, 'int8');
    
    %"file hour" variable
    hour = 0;
    
    if (nargin == 4)
        waitbar(0,wait_dlg,'Writing goGPS binary data...')
    end
    
    %write output files
    for t = 1 : length(time_GPS)
        
        if (nargin == 4)
            waitbar(t/length(time_GPS),wait_dlg)
        end
        
        %-------------------------------------
        % file management
        %-------------------------------------
        
        if (floor(t/3600) > hour)
            
            hour = floor(t/3600);
            hour_str = num2str(hour,'%03d');
            
            fclose(fid_obs);
            fclose(fid_eph);
            
            fid_obs    = fopen([filerootOUT '_obs_'    hour_str '.bin'],'w+');
            fid_eph    = fopen([filerootOUT '_eph_'    hour_str '.bin'],'w+');
            
        end
        
        Eph_t = Eph(:,:,t);
        fwrite(fid_obs, [time_GPS(t); time_M(t); time_R(t); week_R(t); pr1_M(:,t); pr1_R(:,t); ph1_M(:,t); ph1_R(:,t); dop1_R(:,t); snr_M(:,t); snr_R(:,t); pos_M(:,t); iono(:,t)], 'double');
        fwrite(fid_eph, [time_GPS(t); Eph_t(:)], 'double');
    end
    
    %close files
    fclose(fid_obs);
    fclose(fid_eph);
else
    if (nargin == 4)
        msgbox('Both rover and master datasets are required to merge goGPS binary data.');
    else
        fprintf('Both rover and master datasets are required to merge goGPS binary data.\n');
    end
end
