function RINEX2goGPSbin(filename_R_obs, filename_R_nav, filename_M_obs, filename_M_nav, filerootOUT, wait_dlg)

% SYNTAX:
%   RINEX2goGPSbin(filename_R_obs, filename_R_nav, filename_M_obs, filename_M_nav, filerootOUT, wait_dlg);
%
% INPUT:
%   filename_R_obs = rover observation RINEX file
%   filename_R_nav = rover navigation RINEX file
%   filename_M_obs = master observation RINEX file
%   filename_M_nav = master navigation RINEX file
%   filerootOUT = output file root
%   wait_dlg = optional handler to waitbar figure
%
% OUTPUT:
%
% DESCRIPTION:
%   File conversion from rover and master RINEX files to goGPS binary
%   format (*_obs_* and *_eph_* files).

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 beta
%
% Copyright (C) 2009-2010 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Graduate School for Creative Cities, Osaka City University, Japan
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

if (~isempty(dir(filename_R_obs)) & ~isempty(dir(filename_M_obs)) & (~isempty(dir(filename_R_nav)) | ~isempty(dir(filename_M_nav))))
    
    if (isempty(dir(filename_R_nav)))
        filename_R_nav = filename_M_nav;
    elseif (isempty(dir(filename_M_nav)))
        filename_M_nav = filename_R_nav;
    end
    
    %ROVER and MASTER RINEX files reading
    if (nargin == 6)
        [pr1_R, pr1_M, ph1_R, ph1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
            Eph_R, Eph_M, iono_R, iono_M, snr_R, snr_M, ...
            pr1_RR, pr1_MR, ph1_RR, ph1_MR, pr2_RR, pr2_MR, ph2_RR, ph2_MR, ...
            Eph_RR, Eph_MR, snr_RR, snr_MR, ...
            time_GPS, date, pos_M] = ...
            load_RINEX(filename_R_obs, filename_R_nav, filename_M_obs, filename_M_nav, wait_dlg); %#ok<ASGLU>
    else
        [pr1_R, pr1_M, ph1_R, ph1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
            Eph_R, Eph_M, iono_R, iono_M, snr_R, snr_M, ...
            pr1_RR, pr1_MR, ph1_RR, ph1_MR, pr2_RR, pr2_MR, ph2_RR, ph2_MR, ...
            Eph_RR, Eph_MR, snr_RR, snr_MR, ...
            time_GPS, date, pos_M] = ...
            load_RINEX(filename_R_obs, filename_R_nav, filename_M_obs, filename_M_nav); %#ok<ASGLU>
    end
    
    %GPS time
    time_R = time_GPS;
    time_M = time_GPS;
    
    %GPS week number
    date(:,1) = date(:,1) + 2000;
    week_R = floor((datenum(date) - datenum([1980,1,6,0,0,0]))/7);
    
    %select ephemerides source
    if (~Eph_M)
        Eph_tmp = Eph_R;
    else
        Eph_tmp = Eph_M;
    end

    %reconstruct ephemerides complete set
    Eph = zeros(29,32,length(time_GPS));
    for i = 1 : length(time_GPS)
        Eph(:,:,i) = rt_find_eph(Eph_tmp, time_GPS(i));
    end
    
    %select ionosphere parameters source
    if (~iono_M)
        iono(1,1:length(time_GPS)) = iono_R(1);
        iono(2,1:length(time_GPS)) = iono_R(1);
        iono(3,1:length(time_GPS)) = iono_R(1);
        iono(4,1:length(time_GPS)) = iono_R(1);
        iono(5,1:length(time_GPS)) = iono_R(1);
        iono(6,1:length(time_GPS)) = iono_R(1);
        iono(7,1:length(time_GPS)) = iono_R(1);
        iono(8,1:length(time_GPS)) = iono_R(1);
    else
        iono(1,1:length(time_GPS)) = iono_M(1);
        iono(2,1:length(time_GPS)) = iono_M(1);
        iono(3,1:length(time_GPS)) = iono_M(1);
        iono(4,1:length(time_GPS)) = iono_M(1);
        iono(5,1:length(time_GPS)) = iono_M(1);
        iono(6,1:length(time_GPS)) = iono_M(1);
        iono(7,1:length(time_GPS)) = iono_M(1);
        iono(8,1:length(time_GPS)) = iono_M(1);
    end
    
    %master position
    pos_M(1,1:length(time_GPS)) = pos_M(1);
    pos_M(2,1:length(time_GPS)) = pos_M(2);
    pos_M(3,1:length(time_GPS)) = pos_M(3);
    
    %remove epochs in which not all ephemerides are available or available
    %satellites are < 4
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
        snr_R(:,1)  = [];
        snr_M(:,1)  = [];
        pos_M(:,1)  = [];
        Eph(:,:,1)  = [];
        iono(:,1)   = [];
        
        %satObs_R = find( (pr1_R(:,1) ~= 0) & (ph1_R(:,1) ~= 0) );
        %satObs_M = find( (pr1_M(:,1) ~= 0) & (ph1_M(:,1) ~= 0) );
        satObs = find( (pr1_R(:,1) ~= 0) & (pr1_M(:,1) ~= 0));
        satEph = find(sum(abs(Eph(:,:,1)))~=0);
    end
    
    %remove observations without ephemerides
    for i = 1 : length(time_GPS)
        satEph = find(sum(abs(Eph(:,:,i)))~=0);
        delsat = setdiff(1:32,satEph);
        pr1_R(delsat,i) = 0;
        pr1_M(delsat,i) = 0;
        ph1_R(delsat,i) = 0;
        ph1_M(delsat,i) = 0;
        snr_R(delsat,i) = 0;
        snr_M(delsat,i) = 0;
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
        
        filerootOUT(j+1:j+3) = ['_' num2str(i,'%02d')];
        i = i + 1;
    end
    
    %open output files
    fid_obs = fopen([filerootOUT '_obs_00.bin'],'w+');
    fid_eph = fopen([filerootOUT '_eph_00.bin'],'w+');
    
    %"file hour" variable
    hour = 0;
    
    if (nargin == 6)
        waitbar(0,wait_dlg,'Writing goGPS binary data...')
    end
    
    %write output files
    for t = 1 : length(time_GPS)
        
        if (nargin == 6)
            waitbar(t/length(time_GPS),wait_dlg)
        end
        
        %-------------------------------------
        % file management
        %-------------------------------------
        
        if (floor(t/3600) > hour)
            
            hour = floor(t/3600);
            hour_str = num2str(hour,'%02d');
            
            fclose(fid_obs);
            fclose(fid_eph);
            
            fid_obs    = fopen([filerootOUT '_obs_'    hour_str '.bin'],'w+');
            fid_eph    = fopen([filerootOUT '_eph_'    hour_str '.bin'],'w+');
            
        end
        
        Eph_t = Eph(:,:,t);
        fwrite(fid_obs, [time_GPS(t); time_M(t); time_R(t); week_R(t); pr1_M(:,t); pr1_R(:,t); ph1_M(:,t); ph1_R(:,t); snr_M(:,t); snr_R(:,t); pos_M(:,t); iono(:,t)], 'double');
        fwrite(fid_eph, [time_GPS(t); Eph_t(:)], 'double');
    end
    
    %close files
    fclose(fid_obs);
    fclose(fid_eph);
else
    if (nargin == 6)
        msgbox('Both rover and master RINEX files are required to create goGPS binary data.');
    else
        fprintf('Both rover and master RINEX files are required to create goGPS binary data.\n');
    end
end
