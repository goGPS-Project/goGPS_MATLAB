function [ til_P2 ] = interpolateP2( P2_series,commontime,ipp_lon,ipp_lat,PRN,target_sta )
%COMPUTE_DIFFP2 この関数の概要をここに記述
%   詳細説明をここに記述
isnan_idx=~isnan(ipp_lon(target_sta,:));

timestamps=1:length(commontime);
timestamps=timestamps(isnan_idx);
til_P2=P2_series(PRN,:,target_sta);
n_sta= length(ipp_lon(:,1));

for i=1:length(timestamps)
if timestamps(i)>length(til_P2);break;end
%verify the stations used for interpolation
%notnan stations
plane_idx=~isnan(ipp_lon(:,timestamps(i)));

%exclude L1 stations
plane_idx(target_sta)=0;

    if any(plane_idx)==1
        %idx for all stations
        raw_idx=1:n_sta;

        %stations idx used for interpolation
        interp_idx=intersect(raw_idx,raw_idx(plane_idx));

        %execute interpolation
        [til_P2(timestamps(i)),~]=plane_fitting(ipp_lon(interp_idx,timestamps(i)),ipp_lat(interp_idx,timestamps(i)),squeeze(P2_series(PRN,timestamps(i),interp_idx)),ipp_lon(target_sta,timestamps(i)),ipp_lat(target_sta,timestamps(i)));

    else
        %apply value of 1 epoch ago    
        if i>2 && isnan(til_P2(timestamps(i-1)))==1
            til_P2(timestamps(i))=til_P2(timestamps(i-1));
        else
            til_P2(timestamps(i))=NaN;
        end
    end
end


end

