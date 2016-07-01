function [ ipp_lat,ipp_lon,elR] = IPP_satspec( elev_series,azim_series,name_series,commontime,interp_sta,stations_idx,PRN)
%IPP_SATSPEC
% [output]
% ipp_lat: n_staxlength(commontime)
% ipp_lon: n_staxlength(commontime)
% 
% interpolate the elevation,azimuth,diff_L4 for each satellite

n_sta=length(name_series);

%read position of stations
posr = readrpos(name_series);

%define azimuth and elevation
azel=NaN(n_sta,length(commontime),2);

for i=1:n_sta
    %extract azimuth
    azM = azim_series{i}(PRN,stations_idx(i,:));
    pos = find(azM>180);
    azM(pos) = azM(pos) -360;
    
    %extract elevation
    elR = elev_series{i}(PRN,stations_idx(i,:));
    
    azel(i,:,1) = azM;
    azel(i,:,2) = elR;
end


%compute IPP
ipp_lat=NaN(n_sta,length(commontime));
ipp_lon=NaN(n_sta,length(commontime));
for sta = 1:n_sta
    for i = 1 : length(commontime)
        %[stations(sta).posp(i,:),~]=iono_pierce2(azel(sta,i,:).*pi/180,posr(sta,:),350000);
        [latpp, lonpp] = iono_pierce_point(posr(sta,1), posr(sta,2), azel(sta,i,1).*pi/180, azel(sta,i,2).*pi/180);
        stations(sta).posp(i,:) = [latpp, lonpp];
        stations(sta).posp(i,:)=stations(sta).posp(i,:)*180/pi;
    end
    ipp_lat(sta,:) = stations(sta).posp(:,1);
    ipp_lon(sta,:) = stations(sta).posp(:,2);
end
