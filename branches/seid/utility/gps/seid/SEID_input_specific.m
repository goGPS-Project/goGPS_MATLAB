function [ L1_series,L2_series,P2_series,elev_series,azim_series,time_series,name_series,L1_sta] = SEID_input_specific(staname,navfile,constellations,L1_flag)
%UNTITLED この関数の概要をここに記述
%   詳細説明をここに記述

%check the number of files to be read
nmax = size(staname,1);

flag_SP3_found = 0;

[Eph, iono] = load_RINEX_nav(navfile, constellations, flag_SP3_found);

%retrieve multi-constellation wavelengths
lambda = goGNSS.getGNSSWavelengths(Eph, constellations.nEnabledSat);

%load the RINEX observation of target stations
[pr1_R, ph1_R, pr2_R, ph2_R, dop1_R, dop2_R, snr1_R, snr2_R, ...
    time_GPS, time_R, week_R, date_R, pos_R, interval, flag_P1] = ...
    load_RINEX_obs(staname, constellations);

%   [pr1_R, ~ ,ph1_R, ~, pr2_R, ~, ph2_R, ~, ...
%           dop1_R, ~, dop2_R, ~, snr1_R, ~, ...
%           snr2_R, ~, ~, time_R, ~, week_R, ~, ...
%           date_R, ~, pos_R, ~, Eph, iono, interval, flag_P1] = ...
%           load_RINEX2(filename_nav, file_target, [], constellations, flag_SP3_found);

%       [pr1_R, pr1_M, ph1_R, ph1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
%           dop1_R, dop1_M, dop2_R, dop2_M, snr1_R, snr1_M, ...
%           snr2_R, snr2_M, time, time_R, time_M, week_R, week_M, ...
%           date_R, date_M, pos_R, pos_M, Eph, iono, interval, flag_P1] = ...
%           load_RINEX2(filename_nav, file_target, [], constellations, flag_SP3_found);

[pr1_R, ph1_R, pr2_R, ph2_R,~,~,~]=fix_clock_resets_time_code_phase0930(pr1_R, ph1_R, pr2_R, ph2_R, Eph, iono, snr1_R, time_R, lambda);
%save parameters for target station as stracture

if L1_flag==1
    
    phase = 1;
    
    %code and phase
    L1_sta.pr1_R=pr1_R;
    L1_sta.pr2_R=pr2_R;
    L1_sta.ph1_R=ph1_R;
    L1_sta.ph2_R=ph2_R;
    
    %dop and snr
    L1_sta.dop1_R=dop1_R;
    L1_sta.dop2_R=dop2_R;
    L1_sta.snr1_R=snr1_R;
    L1_sta.snr2_R=snr2_R;
    
    %time, intervel and flags
    L1_sta.time_R=time_R;
    L1_sta.week_R=week_R;
    L1_sta.date_R=date_R;
    L1_sta.pos_R=pos_R;
    L1_sta.interval=interval;
    L1_sta.flag_P1=flag_P1;
    %        tar_pr1_R = pr1_R; tar_pr2_R = pr2_R; tar_ph1_R = ph1_R; tar_ph2_R = ph2_R;
    %         tar_dop1_R = dop1_R; tar_dop2_R = dop2_R; tar_snr_1 = snr1_R; tar_snr2_R = snr2_R;
    %         tar_time_R = time_R; tar_week_R=week_R; tar_date_R = date_R; tar_pos_R = pos_R; tar_interval = interval;
    %         tar_flag_P1=flag_P1;
else
    phase = 2;
end

%  load_RINEX(filename_nav, filename_R_obs, filename_M_obs, constellations, flag_SP3, wait_dlg)

%     SP3_time = [];
%     SP3_coor = [];
%     SP3_clck = [];
%
%     if (flag_SP3_found)
%         %display message
%         fprintf('Reading SP3 file...\n');
%
%         [SP3_time, SP3_coor, SP3_clck] = load_SP3(filename_nav, time_R, week_R);
%     end

[date, doy] = gps2date(week_R, time_R);

n_epochs = length(time_R);

name_series     = staname; %#ok<*SAGROW>
time_series     = datenum(date);
%         stations(n_sta).interval = interval;
%         stations(n_sta).pos      = pos_R;

L1_series   = NaN(32,n_epochs);
L2_series   = NaN(32,n_epochs);
P2_series = NaN(32,n_epochs);
azim_series = NaN(32,n_epochs);
elev_series = NaN(32,n_epochs);
%         stations(n_sta).PDOP = NaN(1,n_epochs);

for t = 1 : length(time_R)
    
    Eph_t = rt_find_eph(Eph, time_R(t),constellations.nEnabledSat);
    
    %available satellites
    sat0 = find(pr1_R(:,t) ~= 0);
    
    if (numel(sat0) >= 4)
        
        if (any(pos_R))
            flag_XR = 1;
        else
            flag_XR = 0;
        end
        
        %compute satellite azimuth and elevation
        [~, ~, XS, ~, ~, ~, ~, ~, ~, sat, el, az, ~, ~, ~, PDOP] = init_positioning(time_R(t), pr1_R(sat0,t), snr1_R(sat0,t), Eph_t, [], iono, [], pos_R, [], [], sat0, [], lambda(sat0,:), 0, 0, phase, flag_XR, 0);
        
        L1_series(sat,t)   = ph1_R(sat,t);
        L2_series(sat,t)   = ph2_R(sat,t);
        P2_series(sat,t) = pr2_R(sat,t);
        azim_series(sat,t) = az;
        elev_series(sat,t) = el;
        %                 stations(n_sta).PDOP(t)     = PDOP;
    end
end
