function SEID_main(L1_station, L2_stations, nav_file, pcv_file, out_path)

global cutoff snr_threshold weights

cutoff = 15;
snr_threshold = 0;
weights = 1;

% RNX_directory='C:\Users\Iwaki\Desktop\RISH_UBLOX_JAVAD_splitter\UBX\';
% doy='074';
% L1_station= 'UBXB';
% L2_stations ={'0644','0334','0355','0324'}
%  RNX_directory='C:\Data\rnx\5s\before\';
%  doy='121';
%  L1_station= 'RISH';
%  L2_stations ={'0644','0334','0335','0324'}

lambda1 = 0.190293672798365;
lambda2 = 0.244210213424568;

GPS_flag = 1; GLO_flag = 0; GAL_flag = 0; BDS_flag = 0; QZS_flag = 0; SBS_flag = 0;
[constellations] = goGNSS.initConstellation(GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag);

%load L1 RNX files
[series(1).L1,series(1).L2,series(1).P2,series(1).elev,series(1).azim,series(1).time,series(1).name,L1_sta]=SEID_input_specific(L1_station,nav_file,constellations,1);

%load L2 RNX files
for k=1:length(L2_stations)
    [series(k+1).L1,series(k+1).L2,series(k+1).P2,series(k+1).elev,series(k+1).azim,series(k+1).time,series(k+1).name] = SEID_input_specific(L2_stations{k},nav_file,constellations,0);
end
%  
 for k=1:length(L2_stations)+1
     L1_series{k}=series(k).L1;
     L2_series{k}=series(k).L2;
     P2_series{k}=series(k).P2;
     for PRN=1:32
        zero_idx=find(P2_series{k}(PRN,:)==0);
        P2_series{k}(PRN,zero_idx)=NaN;
     end
     
     elev_series{k}=series(k).elev;
     azim_series{k}=series(k).azim;
     time_series{k}=series(k).time;
     name_series{k}=series(k).name;
 end

 %TODO
 name_series{1} = 'GRED';
 name_series{2} = 'CATU';
 name_series{3} = 'COMO';
 name_series{4} = 'MILA';
 name_series{5} = 'MONZ';
 name_series{6} = 'NOVR';
 name_series{7} = 'VARE';
target_sta = 1; 
 
% n_sta=length(name_series);
% for i=1:n_sta
%     if strcmp(L1_station,name_series{i})==1
%     target_sta=i;   
%     end
% end

antenna_name_series={'NONE', 'TPSG3_A1','TPSCR3_GGD','TPSCR3_GGD','LEIAS10','LEIAR25.R3','LEIAR25'};
L1_series_fix=L1_series;L2_series_fix=L2_series;
% compute diff_L4
[diff_L4,P2_time_series,commontime,stations_idx,~,~,L4_series] = compute_diffL4( L1_series_fix,L2_series_fix,P2_series,name_series,time_series,elev_series,antenna_name_series,pcv_file);

til_L2=NaN(size(L2_series{target_sta}));
til_P2=til_L2;
fix_til_L2=til_L2;
%  Interpolate dL4
for PRN=1:32
%compute IPP
[ satel(PRN).ipp_lat,satel(PRN).ipp_lon,satel(PRN).elR] = IPP_satspec( elev_series,azim_series,name_series,commontime,L1_station,stations_idx,PRN);

%interpolate dL4 and compute ~L4
  [ satel(PRN).til_L4 ] = Planefit_satspec_dL4(diff_L4,commontime,satel(PRN).ipp_lon,satel(PRN).ipp_lat,PRN,target_sta );

%interpolate dP2 and compute ~P2
   [ til_P2(PRN,stations_idx(target_sta,:)) ] = interpolateP2(P2_time_series,commontime,satel(PRN).ipp_lon,satel(PRN).ipp_lat,PRN,target_sta );
 %fix til_L4
%   satel(PRN).fix_tilL4= dL4_fix( satel(PRN).til_L4,elev_series{target_sta}(:,stations_idx(target_sta,:)),commontime,PRN  );
 
 %compute ~L2
  til_L2(PRN,stations_idx(target_sta,:)) = (L1_series_fix{target_sta}(PRN,stations_idx(target_sta,:))*lambda1 - satel(PRN).til_L4)/lambda2;
 
 %compute fix ~L2 (remove huge outlier)
 fix_til_L2(PRN,:)=fix_jump_L2(til_L2,PRN,0.6*10e7);

end
% 
% write CS fixed non-interpolated RINEX file

% outputfile_path=strcat(out_path,'/SEID_output_only_clock_fix.obs') ;
% 
% write_RINEX_obs(outputfile_path, 'u-blox', 'NONE', L1_sta.pr1_R, ...
%     L1_sta.pr2_R, L1_series_fix{target_sta}, L2_series_fix{target_sta}, L1_sta.dop1_R, L1_sta.dop2_R, ...
%     L1_sta.snr1_R, L1_sta.snr2_R, L1_sta.time_R, L1_sta.date_R, L1_sta.pos_R, L1_sta.interval,L1_sta.flag_P1);


%write new RINEX file
outputfile_path=strcat(out_path,['/' name_series{target_sta} '_SEID_output.obs']) ;

new_interval = 30;

write_RINEX_obs(outputfile_path, 'u-blox', 'NONE', L1_sta.pr1_R, ...
    til_P2, L1_series_fix{target_sta}, fix_til_L2, L1_sta.dop1_R, L1_sta.dop2_R, ...
    L1_sta.snr1_R, L1_sta.snr2_R, L1_sta.time_R, L1_sta.date_R, L1_sta.pos_R, new_interval,L1_sta.flag_P1);

undersamplingRINEX(outputfile_path, [L1_station '_new'], 0, new_interval, L1_sta.interval);
