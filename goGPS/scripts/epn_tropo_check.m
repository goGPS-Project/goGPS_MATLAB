%%%% comaprison file 
%close all
project_path = '../../remote_data/project/EPN_test_set';
%goGPS([project_path '/config/config.ini'],false);
getResults();
asi = parseMultiStationTropoSinex([project_path '/station/eur19757.tro']);
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
        t2 = t2.getGpsTime();
        idx_el = t1 < min(t2) | t1 > max(t2);
        t1(idx_el) = [];
        res.TROTOT(idx_el) = [];
        
        
        ztd_diff_t = timeSeriesComparison(t1,res.TROTOT/1000,t2,ztd,mode);
        ztd_diff = [ztd_diff; ztd_diff_t];
        if isfield(res,'TGNTOT')
            res.TGNTOT(idx_el) = [];
            res.TGETOT(idx_el) = [];
            [tgn, tge, ~] = rec(r).getGradient();
            tge_diff = [tge_diff; timeSeriesComparison(t1,res.TGNTOT/1000,t2,tgn,mode)];
            tgn_diff = [tgn_diff; timeSeriesComparison(t1,res.TGETOT/1000,t2,tge,mode)];
        end
        sta_id = [sta_id; r*ones(size(ztd_diff_t))];
        xyz_diff = [xyz_diff; res.xyz - rec(r).xyz];
        figure
        plot(t1,res.TROTOT/1000)
        title(['ZTD ' sta_code]);
        hold on; plot(t2,ztd,'.')        
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
figure; hist(ztd_diff(abs(ztd_diff)<0.05),30)