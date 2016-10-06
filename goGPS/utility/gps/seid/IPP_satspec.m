function [ipp_lat, ipp_lon, elS] = IPP_satspec(elev_series, azim_series, commontime, stations_idx, PRN, pos_RM)

n_sta=length(elev_series);

ipp_lat = NaN(n_sta,length(commontime));
ipp_lon = NaN(n_sta,length(commontime));

for s = 1 : n_sta
    %extract azimuth
    azS = azim_series{s}(PRN,stations_idx(s,:));
    pos = find(azS > 180);
    azS(pos) = azS(pos) - 360;
    
    %extract elevation
    elS = elev_series{s}(PRN,stations_idx(s,:));
    
    azS = azS*pi/180;
    elS = elS*pi/180;
    
    [latR, lonR] = cart2geod(pos_RM(1,1,s), pos_RM(2,1,s), pos_RM(3,1,s));
    
    for e = 1 : length(commontime)
        [latpp, lonpp] = iono_pierce_point(latR, lonR, azS(e), elS(e));
        ipp_lat(s,e) = latpp*180/pi;
        ipp_lon(s,e) = lonpp*180/pi;
    end
end
