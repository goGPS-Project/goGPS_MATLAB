function [ pcv1_series,pcv2_series,pcv_idx ] = pcv_compute( pcv1,pcv2, elevation )
%PCV_COMPUTE ‚±‚ÌŠÖ”‚ÌŠT—v‚ð‚±‚±‚É‹Lq
%[func]   : read receiver pcv data interpolated by each 0.1 degree of elevation
%[argin]  : pcv1       = pcv data for L1 (receiver) [1*901]
%           pcv2       = pcv data for L2 (receiver) [don't be used now]
%           elevation  = timeseries of satellite elevation   [1*2880, ordinary] 
%[argout] : pcv_series = timeseries pcv corresponds for elevation variation
%         : pcv_idx    = "not NaN index" in timeseries of elevation 
%[notice] : "elevation" must be a single continuous curve of angles. 

finegrid=0:0.1:90;
finegrid_inv=finegrid(901:-1:1);
pcv1_inv=pcv1(901:-1:1);
pcv2_inv=pcv2(901:-1:1);
%round elevation
round_elev=round(elevation*10)/10;

%find "not NaN index" of elevation
idx_series=1:1:length(elevation);
pcv_idx=idx_series(~isnan(elevation));

%devide elevation for upper curve and downer curve
rise_elev=round_elev(diff(elevation)>=0);
set_elev=round_elev(diff(elevation)<0);

%find degrees in common
[~,rise_idx1,rise_idx2]=intersect(rise_elev,finegrid);
[~,set_idx1,set_idx2]=intersect(set_elev,finegrid_inv);


rise_pcv1=zeros(1,length(rise_elev));rise_pcv2=rise_pcv1;
set_pcv1=zeros(1,length(set_elev)); set_pcv2=set_pcv1;

rise_pcv1(rise_idx1)=pcv1(rise_idx2);
rise_pcv2(rise_idx1)=pcv2(rise_idx2);

set_pcv1(sort(set_idx1,'ascend'))=pcv1_inv(sort(set_idx2,'ascend'));
set_pcv2(sort(set_idx1,'ascend'))=pcv2_inv(sort(set_idx2,'ascend'));
%apply pcvs for timeseries
pcv1_series=zeros(1,length(elevation));pcv2_series=pcv1_series;
pcv1_series(diff(elevation)>=0)=rise_pcv1;
pcv2_series(diff(elevation)>=0)=rise_pcv2;
pcv1_series(diff(elevation)<0)=set_pcv1;
pcv2_series(diff(elevation)<0)=set_pcv2;

end

