
%load('/Volumes/Data/ArchiveGNSS/Results/RealTime/core_20191023_141610.mat')

%close all;

sta = {'CTWN', 'SUTM'};

for r = 1:numel(sta)
    y_lim = round(minMax(core.rec(r).getZtd * 1e2) + [-1 1]);
    
    base_dir1 = {};
    base_dir2 = {};
    for doy = 268:270
        doy
        if (ispc)
            base_dir1 = [base_dir1; {sprintf('V:/ArchiveGNSS/Results/RealTime_Final_Fixed/2019/%03d/', doy)}];
        else
            base_dir1 = [base_dir1; {sprintf('/Volumes/Data/ArchiveGNSS/Results/RealTime_FinalOL/2019/%03d/', doy)}];
        end
%         if (ispc)
%             base_dir1 = [base_dir1; {sprintf('V:/ArchiveGNSS/Results/RealTime/2019/%03d/', doy)}];
%         else
%             base_dir1 = [base_dir1; {sprintf('/Volumes/Data/ArchiveGNSS/Results/RealTime/2019/%03d/', doy)}];
%         end
        %base_dir2 = [base_dir2; {sprintf('/Volumes/Data/goGPS_data/project/Lampo/OUT_SID/2019/%03d/', doy)}];
    end
    
    tm = Tropo_Mat_Reader(base_dir1, sta{r});
    
    n_sample = min(numel(tm.data_set), 1e100);
    diff_ztd = nan(12*60*2, n_sample);
    
    % load sample 1
    for s = 1:n_sample
        %Per ogni processing...
        [ztd_final, t_final] = core.rec(r).getZtd();
        t_final.toUtc();
        t_final = t_final.getMatlabTime();
        t_final = round(t_final * 86400);
        t_nrt = round(tm.data_set(s).utc_time * 86400);
        
        [t_common, ia, ib] = intersect(t_final, t_nrt);
        dt = (t_common(end) - t_common);
        dt_epoch = -(dt / median(diff(dt), 'omitnan')) + 1;
        %per ogni ora..
        
      
%         figure(1); clf; plot(t_nrt(ib(1):ib(end)), tm.data_set(s).ztd(ib(1):ib(end)));
%         hold on; plot(t_final(ia(1):ia(end)), ztd_final(ia(1):ia(end)));


% export figureZTD
%         figure(1000); hold off; plot(t_nrt(ib(1):ib(end))/86400, tm.data_set(s).ztd(ib(1):ib(end))*1e2, 'Color', Core_UI.getColor(3, 7));
%         hold on; plot(t_final(ia(1):ia(end))/86400, ztd_final(ia(1):ia(end))*1e2, 'Color', Core_UI.getColor(4, 7));
%         axis tight;
%         setAllLinesWidth(2);
%         setTimeTicks(4);
%         ylim(y_lim);
%         title([sta{r} ' ZTD']);
%         legend('NRT', 'FINAL');
%         Core_UI.beautifyFig(gcf,'light');
%         file_name = sprintf('ZTD_comparison_%s_%s.png', sta{r}, datestr(t_final(ia(1)), 'yyyymmddHHMM'));
%         Core_Utils.exportCurFig(file_name);
%       
%         drawnow;
        diff_ztd(dt_epoch, s) = tm.data_set(s).ztd(ib) - ztd_final(ia);
    end
    
     ztd_std = sqrt(mean(diff_ztd.^2, 2, 'omitnan')) * 1e2;
     ztd_mean = mean(diff_ztd, 2, 'omitnan') * 1e2;
     ztd_median = median(diff_ztd, 2, 'omitnan') * 1e2;
     ztd_rmse = sqrt(ztd_mean.^2 + ztd_std.^2);
%     ztd_nmad = mean(abs(diff_ztd * 1e2 - ztd_mean), 2, 'omitnan') / max(diff_ztd(:) * 1e2);
%     
     ztd_max = max(diff_ztd') * 1e2;
     ztd_min = min(diff_ztd') * 1e2;
%     
     ztd_amax = max(abs(diff_ztd)') * 1e2;
     ztd_amin = min(abs(diff_ztd)') * 1e2;
      
%       i = 1;
%       i = i+1; figure(i); plot(t_nrt(ib(1):ib(end))/86400, tm.data_set(s).ztd(ib(1):ib(end)));
%       title('ZTD');
%       Core_UI.beautifyFig(gcf,'light');
%       hold on; plot(t_final(ia(1):ia(end))/86400, ztd_final(ia(1):ia(end)));
%       title('ZTD');
%       setAllLinesWidth(2);
%       Core_UI.beautifyFig(gcf,'light');    
%       setTimeTicks(4);
%     
    i = 1;
    i = i + 1; figure(i); plot((1:size(diff_ztd,1)) / 120, ztd_std); ylabel('std [cm]'); xlabel('hours');
    %title(core.rec(r).getMarkerName4Ch());
    title('STANDARD DEVIATION');
    %ylim([0 2]);
    setAllLinesWidth(2);
    ax(i) = gca;
    Core_UI.beautifyFig(gcf,'light');
    hold on;
    
    i = i + 1; figure(i); plot((1:size(diff_ztd,1)) / 120, ztd_mean); ylabel('mean [cm]'); xlabel('hours');
    %title(core.rec(r).getMarkerName4Ch());
    title('MEAN');
    ylim([-0.4 0.3]);
    setAllLinesWidth(2);
    ax(i) = gca;
    Core_UI.beautifyFig(gcf,'light');
    hold on;
%     
    i = i + 1; figure(i); plot((1:size(diff_ztd,1)) / 120, ztd_median); ylabel('median [cm]'); xlabel('hours');
    %title(core.rec(r).getMarkerName4Ch());
    title('MEDIAN');
    ylim([-0.3 0.4]);
    setAllLinesWidth(2);
    ax(i) = gca;
    Core_UI.beautifyFig(gcf,'light');
    hold on;
%     
    i = i + 1; figure(i); plot((1:size(diff_ztd,1)) / 120, ztd_rmse); ylabel('rmse [cm]');
    %title(core.rec(r).getMarkerName4Ch());
    title('RMSE');
    setAllLinesWidth(2);
    ax(i) = gca;
    Core_UI.beautifyFig(gcf,'light');
    hold on;
%     
% %     i = i + 1; figure(i); plot((1:size(diff_ztd,1)) / 120, ztd_nmad); ylabel('nmad [cm]');
% %     title(core.rec(r).getMarkerName4Ch());
% %     setAllLinesWidth(2);
% %     ax(i) = gca;
% %     Core_UI.beautifyFig(gcf,'light');
% %     hold on;
%     
%     % i = i + 1; figure(i); plot((1:size(diff_ztd,1)) / 120, ztd_max); ylabel('max [cm]');
%     % title(core.rec(r).getMarkerName4Ch());
%     % setAllLinesWidth(2);
%     % ax(i) = gca;
%     % Core_UI.beautifyFig(gcf,'light');
%     %
%     % i = i + 1; figure(i); plot((1:size(diff_ztd,1)) / 120, ztd_min); ylabel('min [cm]');
%     % title(core.rec(r).getMarkerName4Ch());
%     % setAllLinesWidth(2);
%     % ax(i) = gca;
%     % Core_UI.beautifyFig(gcf,'light');
%     %
    i = i + 1; figure(i); plot((1:size(diff_ztd,1)) / 120, ztd_amax); ylabel('abs max [cm]'); xlabel('hours');
    %title(core.rec(r).getMarkerName4Ch());
    title('MAXIMUM');
    ylim([0 7]);
    setAllLinesWidth(2);
    ax(i) = gca;
    Core_UI.beautifyFig(gcf,'light');
    hold on;
%     %
%     i = i + 1; figure(i); plot((1:size(diff_ztd,1)) / 120, ztd_amin); ylabel('abs min [cm]'); xlabel('hours');
%     %title(core.rec(r).getMarkerName4Ch());
%     title('MINIMUM');
%     setAllLinesWidth(2);
%     ax(i) = gca;
%     Core_UI.beautifyFig(gcf,'light');
%     hold on;
    
    linkaxes(ax, 'x');
end
for f = 1:i
    figure(f);
    xlim([0 6]);
    legend(sta, 'Location', 'NorthEast')
end
dockAllFigures;
