%plot satellite orbits

%RINEX files
filename_obs = '../data/data_RINEX/test_ephemeris/RISH1930.12o';
filename_nav = '../data/data_RINEX/test_ephemeris/brdc1930.12n';
filename_sp3 = '../data/data_RINEX/test_ephemeris/igs16963.sp3';
filename_pco = '../data/antenna/ATX/I08.ATX';

obs_comb = 'NONE';
frequencies = [1];

%enabled constellations
GPS_flag = 1;  GLO_flag = 0;  GAL_flag = 0;  BDS_flag = 0;  QZS_flag = 0; SBS_flag = 0;
GPS_col = 'b'; GLO_col = 'r'; GAL_col = 'g'; BDS_col = 'c'; QZS_col = 'm';

color = GPS_col;

sec_day = 86400;

[constellations] = goGNSS.initConstellation(GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag);
nSatTot = constellations.nEnabledSat;

[pr1_R, ph1_R, pr2_R, ph2_R, dop1_R, dop2_R, snr1_R, snr2_R, time_GPS, time_R, week_R, date_R, pos_R, interval, antoff_R, antmod_R, codeC1_R, marker_R] = load_RINEX_obs(filename_obs, constellations); %#ok<*ASGLU>
Eph = load_RINEX_nav(filename_nav, constellations, 0, 0, 0);
% SP3 = load_SP3(filename_sp3, time_GPS, week_R, constellations);
SP3 = [];

%read receiver antenna phase center offset (PCO) and variation (PCV)
antenna_PCV = read_antenna_PCV(filename_pco, antmod_R);

%read satellite antenna phase center offset (NOTE: reading only L1 offset for now)
antmod_S = sat_antenna_ID(constellations);
antenna_PCV_S = read_antenna_PCV(filename_pco, antmod_S, date_R);

lambda = goGNSS.getGNSSWavelengths(Eph, SP3, nSatTot);

% %store satellite antenna PCO/PCV and satellite type
% SP3.antPCO = zeros(1,3,size(antenna_PCV_S,2));
% SP3.satType = cell(1,size(antenna_PCV_S,2));
% for sat = 1 : size(antenna_PCV_S,2)
%     if (antenna_PCV_S(sat).n_frequency ~= 0)
%         SP3.antPCO(:,:,sat) = antenna_PCV_S(sat).offset(:,:,1);
%         SP3.satType{1,sat} = antenna_PCV_S(sat).type;
%     else
%         SP3.avail(sat) = 0;
%     end
% end
%
% %compute sun and moon position
% fprintf('Computing Sun and Moon position...');
% [X_sun, X_moon] = sun_moon_pos(datevec(gps2utc(datenum(date_R))));
% fprintf(' done\n');
%
% %store the position of Sun and Moon
% SP3.t_sun  = time_GPS;
% SP3.X_sun  = X_sun;
% SP3.X_moon = X_moon;
%
% DCB = load_dcb(state.dcb_dir, week_R, time_R, codeC1_R, constellations);
%
% %if DCB files are not available or not sufficient, try to download them
% if ((~any(DCB.P1C1.value(:)) | ~any(DCB.P1P2.value(:))) && constellations.GPS.enabled)
%
%     %download
%     [file_dcb, compressed] = download_dcb([week_R(1) week_R(end)], [time_R(1) time_R(end)]);
%
%     if (compressed)
%         return
%     end
%
%     %try again to read DCB files
%     DCB = load_dcb(state.dcb_dir, week_R, time_R, codeC1_R, constellations);
% end
%
% SP3.DCB = DCB;

%computation interval [s]
% interval = 30;

%approximate satellite-receiver ranges
% pr1_R = ones(nSatTot,1)*2.5e7;

week = median(Eph(24,:));
days = zeros(size(Eph,2),1);
for d = 1 : size(Eph,2);
    dd = gps2date(week,Eph(21,d));
    days(d) = dd(3);
end
day = median(days);
idx = find(abs(days-day) < 2);
% week = median(Eph(24,:));
% idx = find(Eph(24,:) == week);
fit_interval = median(Eph(29,idx));
if (fit_interval == 0)
    fit_interval = 4;
end
time_start = min(Eph(32,idx)) - fit_interval*3600/2;
time_end = max(Eph(32,idx)) + fit_interval*3600/2;

% figure('Color','white')
% hold on
% grid on
% xlim([-180 180])
% ylim([-90 90])

time = time_start : interval : time_end;
nEpochs = length(time);

XS = zeros(nSatTot,3,nEpochs);
phi  = zeros(nSatTot,nEpochs);
lam  = zeros(nSatTot,nEpochs);
h    = zeros(nSatTot,nEpochs);

pr1_R(pr1_R==0) = NaN;
dist = mean(pr1_R,2,'omitnan');

for e = 1 : nEpochs

    XS(:,:,e) = satellite_positions(time(e), dist, 1:nSatTot, Eph, SP3, [], zeros(nSatTot,1), zeros(nSatTot,1), 0, frequencies, obs_comb, lambda);

    for s = 1 : nSatTot
        [phi(s,e), lam(s,e), h(s,e)] = cart2geod(XS(s,1,e), XS(s,2,e), XS(s,3,e));
    end
end

% lam = lam*180/pi;
% phi = phi*180/pi;
% lam(lam == 0) = nan;
% phi(phi == 0) = nan;
%
% % coltab = jet(2*nSatTot);
% m_proj('miller','lat',82);
% m_coast('color',[0 .6 0]);
% hold on
% for s = 1 : nSatTot
%     if (any(~isnan(lam(s,:))))
% %         if (max(lam(s,:)) - min(lam(s,:)) > 180)
% %             [lam(s,:), idx] = sort(lam(s,:));
% %             phi(s,:) = phi(s,idx);
% %         end
% if (constellations.BeiDou.enabled && nSatTot == 37 && s < 6)
%     marker_size = 10;
% else
%     marker_size = 4;
% end
%
%         m_plot(lam(s,:),phi(s,:),'Marker','.','MarkerSize',marker_size,'Color',color,'LineStyle','none')
% %     h = plot(lam(s,:),phi(s,:),'.','MarkerSize',4);
% %     set(h,'Color',coltab(2*s-1,:));
%     end
% end
% m_grid('box','fancy','tickdir','in');
