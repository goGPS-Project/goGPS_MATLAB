function [ ipp_lat,ipp_lon,elR] = IPP_satspec( elev_series,azim_series,name_series,commontime,interp_sta,stations_idx,PRN)
%IPP_SATSPEC この関数の概要をここに記述
% [output]
% ipp_lat: n_staxlength(commontime)
% ipp_lon: n_staxlength(commontime)
% 
% interpolate the elevation,azimuth,diff_L4 for each satellite
%   詳細説明をここに記述

n_sta=length(name_series);

%decide station for interpolation 
for i=1:n_sta
    if strcmp(interp_sta,name_series{i})==1
    target_sta=i;   
    end
end

%read position of stations
 posr = readrpos(name_series);
 
 %should be modified before launch
 
 
 %define azimuth and elevation
 azel=NaN(n_sta,length(commontime),2);
 isnan_idx=zeros(n_sta,length(commontime));
 
 for i=1:n_sta
 %extract azimuth    
 azM = azim_series{i}(PRN,stations_idx(i,:));
 pos = find(azM>180);
 azM(pos) = azM(pos) -360;

 %extract elevation
 elR = elev_series{i}(PRN,stations_idx(i,:));

 azel(i,:,1) = azM;
 azel(i,:,2) = elR;
 
 %find idx which isn't hole
 isnan_idx(i,:)=~isnan(azM);
 end


 %find epoch requires fixing
%  for i=1:n_sta
%     stations(i).interp_id=find((isnan_idx(target_sta,:)-isnan_idx(i,:)==1));
%  end

%  azel_fix=azel;
%  
%  for i=1:n_sta
%     if i~=target_sta;
%         for k=1:length(stations(i).interp_id)
%             
%             if any(isnan(azel_fix(i,stations(i).interp_id(k):stations(i).interp_id(k)+10,1)))==1
%             %determine timeseries for polyfit
%             x_raw=commontime(stations(i).interp_id(k)-10:stations(i).interp_id(k));
%             %determine elevations for polyfit
%             y1_azim=azel_fix(i,stations(i).interp_id(k)-10:stations(i).interp_id(k),1);
%             y1_elev=azel_fix(i,stations(i).interp_id(k)-10:stations(i).interp_id(k),2);
%             
%             else 
%             x_raw=commontime(stations(i).interp_id(k):stations(i).interp_id(k)+10);
%             %determine elevations for polyfit
%             y1_azim=azel_fix(i,stations(i).interp_id(k):stations(i).interp_id(k)+10,1);
%             y1_elev=azel_fix(i,stations(i).interp_id(k):stations(i).interp_id(k)+10,2);
%             
%             end
%             x=x_raw(~isnan(y1_azim));
%             
%             %polyfit azimuth
%             y_azim=y1_azim(~isnan(y1_azim));
%             p_azim=polyfit(x,y_azim',2);
%             azel_fix(i,stations(i).interp_id(k),1)=polyval(p_azim,commontime(stations(i).interp_id(k)));
%             
%             %polyfit elev
%             y_elev=y1_elev(~isnan(y1_azim));
%             p_elev=polyfit(x,y_elev',2);
%             azel_fix(i,stations(i).interp_id(k),1)=polyval(p_elev,commontime(stations(i).interp_id(k)));
%         end
%     end
%  end
 
%compute IPP
ipp_lat=NaN(n_sta,length(commontime));
ipp_lon=NaN(n_sta,length(commontime));
 for sta = 1:n_sta
        for i = 1 : length(commontime)
%          [stations(sta).posp(i,:),~]=iono_pierce2(azel(sta,i,:).*pi/180,posr(sta,:),350000);
             [stations(sta).posp(i,:),~]=iono_pierce2(azel(sta,i,:).*pi/180,posr(sta,:),350000);
             stations(sta).posp(i,:)=stations(sta).posp(i,:)*180/pi;
        end
%      ipp_lat(sta,:) = stations(sta).posp(:,1)*180/pi;
     ipp_lat(sta,:) = stations(sta).posp(:,1);
     ipp_lon(sta,:) = stations(sta).posp(:,2);
%        ipp_lon(sta,:) = stations(sta).posp(:,2)*180/pi;
 end

end

