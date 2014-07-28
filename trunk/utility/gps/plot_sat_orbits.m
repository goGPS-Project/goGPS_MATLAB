%plot satellite orbits

%RINEX files
% filename_obs = '../data/data_RINEX/multiGNSS_test/gmsd1800.14o';
filename_nav = '../data/data_RINEX/multiGNSS_test/brdm1800.14p';

%enabled constellations
GPS_flag = 1;  GLO_flag = 0;  GAL_flag = 0;  BDS_flag = 0;  QZS_flag = 0; SBS_flag = 0;
GPS_col = 'b'; GLO_col = 'r'; GAL_col = 'g'; BDS_col = 'c'; QZS_col = 'm';

color = GPS_col;

sec_day = 86400;

[constellations] = goGNSS.initConstellation(GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag);
nSatTot = constellations.nEnabledSat;

% [pr1_R, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, interval, ~, ~] = load_RINEX_obs(filename_obs, constellations);
Eph = load_RINEX_nav(filename_nav, constellations, 0);

%computation interval [s]
interval = 30;

%approximate satellite-receiver ranges
pr1_R = ones(nSatTot,1)*2.5e7;

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

for e = 1 : nEpochs
    
    dist = mean(pr1_R,2);

    XS(:,:,e) = satellite_positions(time(e), dist, 1:nSatTot, Eph, [], [], zeros(nSatTot,1), zeros(nSatTot,1), 0);
    
    for s = 1 : nSatTot
        [phi(s,e), lam(s,e), h(s,e)] = cart2geod(XS(s,1,e), XS(s,2,e), XS(s,3,e));
    end
end

lam = lam*180/pi;
phi = phi*180/pi;
lam(lam == 0) = nan;
phi(phi == 0) = nan;

% coltab = jet(2*nSatTot);
m_proj('miller','lat',82);
m_coast('color',[0 .6 0]);
hold on
for s = 1 : nSatTot
    if (any(~isnan(lam(s,:))))
%         if (max(lam(s,:)) - min(lam(s,:)) > 180)
%             [lam(s,:), idx] = sort(lam(s,:));
%             phi(s,:) = phi(s,idx);
%         end
if (constellations.BeiDou.enabled && nSatTot == 37 && s < 6)
    marker_size = 10;
else
    marker_size = 4;
end

        m_plot(lam(s,:),phi(s,:),'Marker','.','MarkerSize',marker_size,'Color',color,'LineStyle','none')
%     h = plot(lam(s,:),phi(s,:),'.','MarkerSize',4);
%     set(h,'Color',coltab(2*s-1,:));
    end
end
m_grid('box','fancy','tickdir','in');
