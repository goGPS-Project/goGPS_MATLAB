function [ L1_series,L2_series,P2_series,elev_series,azim_series,time_series,name_series,L1_sta] = SEID_input_specific(staname,navfile,constellations,L1_flag)

%check the number of files to be read
nmax = size(staname,1);

flag_SP3_found = 0;
SP3 = [];

[Eph, iono] = load_RINEX_nav(navfile, constellations, flag_SP3_found);

%retrieve multi-constellation wavelengths
lambda = goGNSS.getGNSSWavelengths(Eph, SP3, constellations.nEnabledSat);

%load the RINEX observation of target stations
[pr1_R, ph1_R, pr2_R, ph2_R, dop1_R, dop2_R, snr1_R, snr2_R, ...
    time_GPS, time_R, week_R, date_R, pos_R, interval, antoff, antmod, flag_C1] = ...
    load_RINEX_obs(staname, constellations);

date_GPS = gps2date(week_R, weektime2tow(week_R, time_GPS));
date_R   = gps2date(week_R, weektime2tow(week_R, time_R));
date = date_GPS;

%[pr1_R, ph1_R, pr2_R, ph2_R,~,~,~]=fix_clock_resets_time_code_phase0930(pr1_R, ph1_R, pr2_R, ph2_R, Eph, iono, snr1_R, time_R, lambda);
fprintf('Observation pre-processing and cycle-slip detection...\n');
[pr1_R, ph1_R, pr2_R, ph2_R] = pre_processing(time_GPS, time_R, pos_R, pr1_R, ph1_R, pr2_R, ph2_R, dop1_R, dop2_R, snr1_R, Eph, SP3, iono, lambda, 1, 'NONE', size(pr1_R,1), [], 0, []);
fprintf('... done.\n');

%save parameters for target station as structure
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
    L1_sta.time_R=time_GPS;
    L1_sta.week_R=week_R;
    L1_sta.date_R=date_R;
    L1_sta.pos_R=pos_R;
    L1_sta.interval=interval;
    L1_sta.flag_C1=flag_C1;
else
    phase = 2;
end

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

n_epochs = length(time_R);

name_series = staname; %#ok<*SAGROW>
time_series = datenum(date);
L1_series   = NaN(32,n_epochs);
L2_series   = NaN(32,n_epochs);
P2_series   = NaN(32,n_epochs);
azim_series = NaN(32,n_epochs);
elev_series = NaN(32,n_epochs);

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
        [~, ~, XS, ~, ~, ~, ~, ~, ~, sat, el, az, ~, ~, ~, PDOP] = init_positioning(time_R(t), pr1_R(sat0,t), snr1_R(sat0,t), Eph_t, [], iono, [], pos_R, [], [], sat0, [], lambda(sat0,:), 0, 0, phase, flag_XR, 0, 0);
        
        L1_series(sat,t)   = ph1_R(sat,t);
        L2_series(sat,t)   = ph2_R(sat,t);
        P2_series(sat,t)   = pr2_R(sat,t);
        azim_series(sat,t) = az;
        elev_series(sat,t) = el;
    end
end
