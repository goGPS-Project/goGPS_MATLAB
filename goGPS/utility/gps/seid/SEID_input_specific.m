function [ L1_series,L2_series,P2_series,elev_series,azim_series,time_series,name_series,antenna_series,L1_sta] = SEID_input_specific(staname,navfile,constellations,L1_flag)

%check the number of files to be read
nmax = size(staname,1);

flag_SP3 = 1;
SP3 = [];

[Eph, iono] = load_RINEX_nav(navfile, constellations, flag_SP3);

%retrieve multi-constellation wavelengths
lambda = goGNSS.getGNSSWavelengths(Eph, SP3, constellations.nEnabledSat);

%load the RINEX observation of target stations
[pr1_R, ph1_R, pr2_R, ph2_R, dop1_R, dop2_R, snr1_R, snr2_R, ...
    time_GPS, time_R, week_R, date_R, pos_R, interval, antoff_R, antmod_R, codeC1_R, marker_R] = ...
    load_RINEX_obs(staname, constellations);

date_GPS = gps2date(week_R, weektime2tow(week_R, time_GPS));
date_R   = gps2date(week_R, weektime2tow(week_R, time_R));
date = date_GPS;

%read satellite antenna phase center offset (NOTE: reading only L1 offset for now)
antmod_S = sat_antenna_ID(constellations);
antenna_PCV_S = read_antenna_PCV('../data/stations/I08.ATX', antmod_S, date_R);

if (flag_SP3)
    %display message
    fprintf('Reading SP3 file...\n');
    
    %----------------------------------------------------------------------------------------------
    % LOAD SP3 DATA
    %----------------------------------------------------------------------------------------------
    SP3 = load_SP3(navfile, time_GPS, week_R, constellations);
    
    %store satellite antenna PCO/PCV and satellite type
    SP3.antPCO = zeros(1,3,size(antenna_PCV_S,2));
    SP3.satType = cell(1,size(antenna_PCV_S,2));
    for sat = 1 : size(antenna_PCV_S,2)
        if (antenna_PCV_S(sat).n_frequency ~= 0)
            SP3.antPCO(:,:,sat) = antenna_PCV_S(sat).offset(:,:,1);
            SP3.satType{1,sat} = antenna_PCV_S(sat).type;
        else
            SP3.avail(sat) = 0;
        end
    end
    
    %compute sun and moon position
    fprintf('Computing Sun and Moon position...');
    [X_sun, X_moon] = sun_moon_pos(datevec(gps2utc(datenum(date_R))));
    fprintf(' done\n');
    
    %store the position of Sun and Moon
    SP3.t_sun  = time_GPS;
    SP3.X_sun  = X_sun;
    SP3.X_moon = X_moon;
    
    %----------------------------------------------------------------------------------------------
    % LOAD DCB DATA (DIFFERENTIAL CODE BIASES)
    %----------------------------------------------------------------------------------------------
    
    %NOTE: if not using SP3 ephemeris or if DCB files are not available, the
    %      'SP3.DCB' structure will be initialized to zero/empty arrays and it will not
    %      have any effect on the positioning
    
    %if (~strcmp(obs_comb, 'IONO_FREE'))
    %try first to read already available DCB files
    DCB = load_dcb('../data/DCB', week_R, time_R, codeC1_R, constellations);
    
    %if DCB files are not available or not sufficient, try to download them
    if (isempty(DCB))
        
        %download
        [file_dcb, compressed] = download_dcb([week_R(1) week_R(end)], [time_R(1) time_R(end)]);
        
        if (compressed)
            return
        end
        
        %try again to read DCB files
        DCB = load_dcb('../data/DCB', week_R, time_R, codeC1_R, constellations);
    end
    
    SP3.DCB = DCB;
    %else
    %SP3.DCB = [];
    %end
end

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
    L1_sta.flag_C1=codeC1_R;
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

name_series = cell2mat(marker_R); %#ok<*SAGROW>
antenna_series = cell2mat(antmod_R);
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
        [~, ~, XS, ~, ~, ~, ~, ~, ~, sat, el, az, ~, ~, ~, PDOP] = init_positioning(time_R(t), pr1_R(sat0,t), snr1_R(sat0,t), Eph_t, SP3, iono, [], pos_R, [], [], sat0, [], lambda(sat0,:), 0, 0, phase, flag_XR, 0, 0);
        
        L1_series(sat,t)   = ph1_R(sat,t);
        L2_series(sat,t)   = ph2_R(sat,t);
        P2_series(sat,t)   = pr2_R(sat,t);
        azim_series(sat,t) = az;
        elev_series(sat,t) = el;
    end
end
