%%%% comaprison file 
close all
project_path = '../../remote_data/project/EPN_test_set'
goGPS([project_path '/config/config.ini'],false);
getResults();
asi = parseMultiStationTropoSinex([project_path '/station/asi19751.tro']);
mode = 'interpolate', % possible values : 'aggregate', 'interpolate'
ztd_diff = [];
tge_diff = [];
tgn_diff = [];
xyz_diff = [];
sta_id = [];
for r = 1 : max(size(rec))
    sta_code = rec(r).getMarkerName4Ch;
    if isfield(asi,sta_code)
        res = asi.(sta_code);
        t1 = res.time.getGpsTime();
        [ztd, t2] = rec(r).getZtd();
        ztd_diff_t = timeSeriesComparison(t1,res.ztd,t2.getGpsTime(),ztd,mode);
        ztd_diff = [ztd_diff; ztd_diff_t];
        [tgn, tge, t2] = rec(r).getGradient();
        tge_diff = [tge_diff; timeSeriesComparison(t1,res.tgn,t2.getGpsTime(),tgn,mode)];
        tgn_diff = [tgn_diff; timeSeriesComparison(t1,res.tge,t2.getGpsTime(),tge,mode)];
        sta_id = [sta_id; r*ones(size(ztd_diff_t))];
        xyz_diff = [xyz_diff; res.xyz - rec(r).xyz];
        figure
        plot(t1,res.ztd)
        title(['ZTD ' sta_code]);
        hold on; plot(t2.getGpsTime,ztd,'.')        
%         figure
%         plot(t1,res.tgn)
%         title(['TGN ' sta_code]);
%         hold on; plot(t2.getGpsTime,tge,'.')
%         figure
%         plot(t1,res.tge)
%         title(['TGE ' sta_code]);
%         hold on; plot(t2.getGpsTime,tge,'.')
%         figure
%         plot(ztd_diff(sta_id == r));
%         title(['ZTD diff' sta_code]);
%         figure
%         plot(tgn_diff(sta_id == r));
%         title(['TGN diff' sta_code]);
%         figure
%         plot(tge_diff(sta_id == r));
%         title(['TGE diff' sta_code]);
    end
end