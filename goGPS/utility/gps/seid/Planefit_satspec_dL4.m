function [ til_L4 ] = Planefit_satspec_dL4(  diff_L4,commontime,ipp_lon,ipp_lat,PRN,target_sta )
%DL4_PLANEFIT この関数の概要をここに記述
%   詳細説明をここに記述

 %find idx to execute interpolation
isnan_idx=~isnan(ipp_lon(target_sta,:));

timestamps=1:length(commontime);
timestamps=timestamps(isnan_idx);
interp_dL4=diff_L4(PRN,:,target_sta);
n_sta= length(ipp_lon(:,1));

for i=1:length(timestamps)
if timestamps(i)>length(interp_dL4);break;end
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
        [interp_dL4(timestamps(i)),~]=plane_fitting(ipp_lon(interp_idx,timestamps(i)),ipp_lat(interp_idx,timestamps(i)),squeeze(diff_L4(PRN,timestamps(i),interp_idx)),ipp_lon(target_sta,timestamps(i)),ipp_lat(target_sta,timestamps(i)));

    else
        %apply value of 1 epoch ago    
        if i>2 && isnan(interp_dL4(timestamps(i-1)))==1
            interp_dL4(timestamps(i))=interp_dL4(timestamps(i-1));
        else
            interp_dL4(timestamps(i))=NaN;
        end
    end
end

%compute sum of ~dL4
 notnan_id=find(isnan_idx==1);
  %check if isnan_idx exceeds size of interp_dL4
     if any(notnan_id>length(interp_dL4))
         isnan_idx(notnan_id(notnan_id>length(interp_dL4)))=0;
     end 
     
    dL4_for_sum=interp_dL4(isnan_idx);  
%  search NaN and replace to zero 
nan_find=isnan(dL4_for_sum);
 if (any(nan_find))==1
      dL4_for_sum(nan_find==1)=0;
 end
    
til_L4=zeros(length(commontime),1);

sum_dL4=zeros(1,length(dL4_for_sum));
for i=2:length(dL4_for_sum)
    sum_dL4(i)=sum(dL4_for_sum(1:i-1));
end

til_L4=NaN(1,length(commontime));
til_L4(isnan_idx)=sum_dL4;

end

