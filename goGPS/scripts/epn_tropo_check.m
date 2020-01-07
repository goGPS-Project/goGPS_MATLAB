%%%% comaprison file 
project_path = '../../remote_data/project/EPN_test_set/station/';
%goGPS([project_path '/config/config.ini'],false);
getResults();
solutions = {'asi19752.tro'  'bek19751.tro'  'bkg19751.tro' 'eur19757.tro'};
ext_sol = parseMultiStationTropoSinex([project_path solutions{1}]);
mode = 'interpolate', % possible values : 'aggregate', 'interpolate'
ztd_diff = [];
tge_diff = [];
tgn_diff = [];
xyz_diff = [];
sta_id = [];
for r = 1 : max(size(rec))
    sta_code = rec(r).getMarkerName4Ch;
    if isfield(ext_sol,sta_code)
        res = ext_sol.(sta_code);
        t1 = res.time.getMatlabTime();
        [ztd] = rec(r).out.getZtd();
        t2 = rec(r).out.getTime().getMatlabTime();
        
        idx_el = t1 < min(t2) | t1 > max(t2);
        t1(idx_el) = [];
        res.TROTOT(idx_el) = [];
        
        
        ztd_diff_t = timeSeriesComparison(t1,res.TROTOT/1000,t2,ztd,mode);
        ztd_diff = [ztd_diff; ztd_diff_t];
        if isfield(res,'TGNTOT')
            res.TGNTOT(idx_el) = [];
            res.TGETOT(idx_el) = [];
            [tgn, tge] = rec(r).out.getGradient();
            tge_diff = [tge_diff; timeSeriesComparison(t1,res.TGNTOT/1000,t2,tgn,mode)];
            tgn_diff = [tgn_diff; timeSeriesComparison(t1,res.TGETOT/1000,t2,tge,mode)];
        end
        sta_id = [sta_id; r*ones(size(ztd_diff_t))];
        xyz_diff = [xyz_diff; res.xyz - rec(r).out.xyz(1)];
        figure
        plot(t1,res.TROTOT/1000)
        title(['ZTD ' sta_code]);
        hold on; plot(t2,ztd,'.')        
        figure
        plot(t1,res.TGNTOT/1000)
        title(['TGN ' sta_code]);
        hold on; plot(t2,tgn,'.')
         setTimeTicks(6,'dd/mm/yyyy HH:MMPM')
        figure
        plot(t1,res.TGETOT/1000)
        title(['TGE ' sta_code]);
        hold on; plot(t2,tge,'.')
         setTimeTicks(6,'dd/mm/yyyy HH:MMPM')
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
figure; hist(ztd_diff(abs(ztd_diff)<0.05),30); title('ZTD diff')
figure; hist(tge_diff(abs(tge_diff)<0.05),30); title('TGE diff')
figure; hist(tgn_diff(abs(tgn_diff)<0.05),30); title('TGN diff')