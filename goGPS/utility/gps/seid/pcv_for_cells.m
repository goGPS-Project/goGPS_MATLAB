function [ L1_series_pcv,L2_series_pcv ] = pcv_for_cells(L1_series,L2_series,name_series,time_series,elev_series,antenna_name_series,filepath_pcv)
%PCV_FOR_CELLS ‚±‚ÌŠÖ?”‚ÌŠT—v‚ð‚±‚±‚É‹L?q
%   ?Ú?×?à–¾‚ð‚±‚±‚É‹L?q
lambda1 = 0.190293672798365;
lambda2 = 0.244210213424568;

%load 'stations' field

 num_sta=length(name_series);
% station_names={num_sta};
%load pcv file

[ants,~,~,pcv1,pcv2,~]=readantex2(filepath_pcv);
grid = linspace(0,90,19);
finegrid=0:0.1:90;
% rec_idx=find(strcmp('JAV_GRANT-G3T',ants)==1);
% rec_pcv1=interp1(grid,pcv1(rec_idx,:),finegrid);
% rec_pcv2=interp1(grid,pcv2(rec_idx,:),finegrid);

L1_series_pcv={};L2_series_pcv=L1_series_pcv;

for target_sta=1:num_sta
    sta(target_sta).pcv1 = zeros(32,length(time_series{target_sta}));
    sta(target_sta).pcv2=sta(target_sta).pcv1;

    rec_idx=find(strcmp(antenna_name_series(target_sta),ants)==1);
    rec_idx=find(strcmp('TRM29659.00',ants)==1);
    rec_pcv1=interp1(grid,pcv1(rec_idx,:),finegrid);
    rec_pcv2=interp1(grid,pcv2(rec_idx,:),finegrid);

        
for i=1:32
diff_elev=diff(elev_series{target_sta}(i,:));
idx_elev=1:1:length(diff_elev);
isnan_elev=idx_elev(~isnan(diff_elev));
    switch length(find(diff(isnan_elev)>100))
        case 0
            [sta(target_sta).pcv1(i,:),sta(target_sta).pcv2(i,:),~]=pcv_compute(rec_pcv1,rec_pcv2,elev_series{target_sta}(i,:));

        case 1
        fake_elev=elev_series{target_sta}(i,:);
        fake_elev(isnan_elev(diff(isnan_elev)>100)+1:end)=NaN;
        [pcv1_series_1,pcv2_series_1,idx1]=pcv_compute(rec_pcv1,rec_pcv2,fake_elev);
        
        fake_elev=elev_series{target_sta}(i,:);
        fake_elev(1:isnan_elev(find(diff(isnan_elev)>100)+1)-1)=NaN;
        [pcv1_series_2,pcv2_series_2,idx2]=pcv_compute(rec_pcv1,rec_pcv2,fake_elev);
        
        sta(target_sta).pcv1(i,idx1)=pcv1_series_1(idx1);sta(target_sta).pcv1(i,idx2)=pcv1_series_2(idx2);
        sta(target_sta).pcv2(i,idx1)=pcv2_series_1(idx1);sta(target_sta).pcv2(i,idx2)=pcv2_series_2(idx2);

%         pcv1_series(i,idx1)=pcv1_series_1(idx1);pcv1_series(i,idx2)=pcv1_series_2(idx2);
%         pcv2_series(i,idx1)=pcv2_series_1(idx1);pcv2_series(i,idx2)=pcv2_series_2(idx2);
        
        case 2
        cut_point=ones(1,4);cut_point(2:3)=isnan_elev(diff(isnan_elev)>100);cut_point(end)=isnan_elev(end);
        fake_elev=elev_series{target_sta}(i,:);fake_elev(cut_point(2)+3:end)=NaN;
        [pcv1_series_1,pcv2_series_1,idx1]=pcv_compute(rec_pcv1,rec_pcv2,fake_elev);
        
         fake_elev=elev_series{target_sta}(i,:);
         fake_elev(cut_point(1):cut_point(2)+3)=NaN;fake_elev(cut_point(3)+2:cut_point(2)+end)=NaN;
        [pcv1_series_2,pcv2_series_2,idx2]=pcv_compute(rec_pcv1,rec_pcv2,fake_elev);
        
        fake_elev=elev_series{target_sta}(i,:);fake_elev(1:cut_point(3)+3)=NaN;
        [pcv1_series_3,pcv2_series_3,idx3]=pcv_compute(rec_pcv1,rec_pcv2,fake_elev);
        
        sta(target_sta).pcv1(i,idx1)=pcv1_series_1(idx1);sta(target_sta).pcv1(i,idx2)=pcv1_series_2(idx2);sta(target_sta).pcv1(i,idx3)=pcv1_series_3(idx3);
        sta(target_sta).pcv2(i,idx1)=pcv2_series_1(idx1);sta(target_sta).pcv2(i,idx2)=pcv2_series_2(idx2);sta(target_sta).pcv2(i,idx3)=pcv2_series_3(idx3);
    end
    
end


L1_series_pcv{target_sta}=(L1_series{target_sta}*lambda1+sta(target_sta).pcv1)/lambda1;
L2_series_pcv{target_sta}=(L2_series{target_sta}*lambda2+sta(target_sta).pcv2)/lambda2;

end


end