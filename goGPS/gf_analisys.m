for r=1:17; gf(r) = core.rec(r).work.getGeometryFree('L1', 'L2', 'G'); end
%%

id_ref = 1;
close all
for id_cmp = [2:17]
    [~, time_ref, time_cmp] = intersect(round(gf(id_ref).time.getRefTime(gf(id_ref).time.first.getMatlabTime())), round(gf(id_cmp).time.getRefTime(gf(id_ref).time.first.getMatlabTime())));
    [~, prn_ref, prn_cmp] = intersect(gf(id_ref).prn, gf(id_cmp).prn);
    
    % geometry-free diff
    tmp = zero2nan(gf(id_ref).obs(time_ref, prn_ref)) - zero2nan(gf(id_cmp).obs(time_cmp,prn_cmp));
    
    % for each satellite
    for s = 1 : size(tmp, 2)
        lim = getOutliers(~isnan(tmp(:,s)));
        for l = 1 : size(lim,1)
            % get limit of the arc
            id_ok = lim(l, 1) : lim(l, 2);
            tmp(id_ok, s) = tmp(id_ok, s) - median(tmp(id_ok, s));
        end
    end
    fh = figure; plot(1e2 * tmp)
    fprintf('Rec %s: %f\n', rec(id_cmp).getMarkerName4Ch, std(1e2 * tmp(:), 'omitnan'))
    ylim(25 * [-1 1])
    iono_diff{id_cmp} = tmp;
    fh.Name = rec(id_cmp).getMarkerName4Ch;
    
end
dockAllFigures();

%% std ZTD
ztd = rec.getZtd_mr;

ztd_diff = ztd(:,2:end) - ztd(:, 1);

for id_cmp = [2:17]
    fprintf('Rec %s: %f\n', rec(id_cmp).getMarkerName4Ch, std(1e2 * ztd_diff(:,id_cmp-1), 'omitnan'))
end

%% plot ZTD @ 25Km
rec([1, 4, 14, 17]).showZtd

%% plot based on elevation
id_ref = 1;
close all
for id_cmp = [2:17]
    [~, time_ref, time_cmp] = intersect(round(gf(id_ref).time.getRefTime(gf(id_ref).time.first.getMatlabTime())), round(gf(id_cmp).time.getRefTime(gf(id_ref).time.first.getMatlabTime())));
    [~, prn_ref, prn_cmp] = intersect(gf(id_ref).prn, gf(id_cmp).prn);
    
    % geometry-free diff
    tmp = zero2nan(gf(id_ref).obs(time_ref, prn_ref)) - zero2nan(gf(id_cmp).obs(time_cmp,prn_cmp));
    
    % for each satellite
    for s = 1 : size(tmp, 2)
        lim = getOutliers(~isnan(tmp(:,s)));
        for l = 1 : size(lim,1)
            % get limit of the arc
            id_ok = lim(l, 1) : lim(l, 2);
            tmp(id_ok, s) = tmp(id_ok, s) - median(tmp(id_ok, s));
            %tmp(id_ok, s) = median(tmp(id_ok, s));
        end
    end
    
    el = gf(id_cmp).el(time_cmp,prn_cmp);
    
    fh = figure; 
    for s = 1 : size(tmp, 2)
        id_ok = ~isnan(tmp(:,s));
        plot(el(id_ok,s), 1e2 * abs(tmp(id_ok,s)), '.'); hold on;
    end
    fh.Name = rec(id_cmp).getMarkerName4Ch;
    ylim([0 15]);
end
dockAllFigures();


% %% plot based on elevation
% id_ref = 1;
% close all
% fh = figure;
% for id_cmp = [3, 13]
%     [~, time_ref, time_cmp] = intersect(round(gf(id_ref).time.getRefTime(gf(id_ref).time.first.getMatlabTime())), round(gf(id_cmp).time.getRefTime(gf(id_ref).time.first.getMatlabTime())));
%     [~, prn_ref, prn_cmp] = intersect(gf(id_ref).prn, gf(id_cmp).prn);
%     
%     % geometry-free diff
%     tmp = zero2nan(gf(id_ref).obs(time_ref, prn_ref)) - zero2nan(gf(id_cmp).obs(time_cmp,prn_cmp));
%     
%     % for each satellite
%     for s = 1 : size(tmp, 2)
%         lim = getOutliers(~isnan(tmp(:,s)));
%         for l = 1 : size(lim,1)
%             % get limit of the arc
%             id_ok = lim(l, 1) : lim(l, 2);
%             tmp(id_ok, s) = tmp(id_ok, s) - 0.2*l - median(tmp(id_ok, s));
%         end
%     end
%     
%     el = gf(id_cmp).el(time_cmp,prn_cmp);
%     
%     for s = 1 : size(tmp, 2)
%         id_ok = ~isnan(tmp(:,s));
%         %plot(el(id_ok,s), 1e2 * (tmp(id_ok,s)), '.'); hold on;
%         plot(1e2 * (tmp(id_ok,s)), '.'); hold on;
%     end
%     fh.Name = rec(id_cmp).getMarkerName4Ch;
%     %ylim([0 15]);
% end
% dockAllFigures();
% 
