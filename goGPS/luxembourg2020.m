%[l, m] = Core_Utils.getAllZdegree(3, 3);
%Core_Utils.showZernike3(l, m, ones(size(l)), 0); colormap(Cmap.get('RdBu', 1024));
%Core_Utils.showAllZernike(l, m, ones(size(l)), 0, true); colormap([1 1 1] * linspace(1,0,1024));


std_gr = [];
[zwd, p_time, id_sync, tge, tgn] = rmap.core.rec.getZtd_mr();
std_gr(1, :) = [mean(std(zwd - mean(zwd,2,'omitnan'),'omitnan')) mean(std(tgn - mean(tgn,2,'omitnan'),'omitnan')) mean(std(tge - mean(tge,2,'omitnan'),'omitnan'))];
[zwd, p_time, id_sync, tge, tgn] = core.rec.getZtd_mr();
std_gr(2, :) = [mean(std(zwd - mean(zwd,2,'omitnan'),'omitnan')) mean(std(tgn - mean(tgn,2,'omitnan'),'omitnan')) mean(std(tge - mean(tge,2,'omitnan'),'omitnan'))];
figure; plot(bsxfun(@rdivide, std_gr, std_gr(1,:)) * 100, 'o-', 'LineWidth', 2);

% Coordinates
for i=1
    rec_set = (1 : 6)
    std_enu = [];
    enu = [];
    for r = rec_set
        enu(:,:,r) = rmap.core.rec(r).out.add_coo(2).coo.getENU;
        enu(:,:,r) = enu(:,:,r) - median(enu(:,:,r), 1, 'omitnan');
    end
    enu = permute(enu, [1 3 2]);
    std_enu(1, :) = [mean(std(enu(:,:,1), 1, 'omitnan')) mean(std(enu(:,:,2), 1, 'omitnan')) mean(std(enu(:,:,3), 1, 'omitnan'))];
    
    enu = [];
    for r = rec_set
        enu(:,:,r) = core.rec(r).out.add_coo(2).coo.getENU;
        enu(:,:,r) = enu(:,:,r) - median(enu(:,:,r), 1, 'omitnan');
    end
    enu = permute(enu, [1 3 2]);
    std_enu(2, :) = [mean(std(enu(:,:,1), 1, 'omitnan')) mean(std(enu(:,:,2), 1, 'omitnan')) mean(std(enu(:,:,3), 1, 'omitnan'))];
    
    figure; plot(bsxfun(@rdivide, std_enu, std_enu(1,:)) * 100, 'o-', 'LineWidth', 2);
    if numel(rec_set) == 1 
        title(sprintf('Coordinates stability improvement of %s', rec(rec_set).getMarkerName4Ch));
    else
        title(sprintf('Coordinates stability improvement'));
    end
end

%% Load cores
none = load('/Users/Andrea/Repositories/goGPS_MATLAB/data/project/Antenna_Calibration/out/core_IF_test_none.mat');
zmap = load('/Users/Andrea/Repositories/goGPS_MATLAB/data/project/Antenna_Calibration/out/core_IF_test_zmap.mat');
rmap = load('/Users/Andrea/Repositories/goGPS_MATLAB/data/project/Antenna_Calibration/out/core_IF_test_rmap.mat');
gmap = load('/Users/Andrea/Repositories/goGPS_MATLAB/data/project/Antenna_Calibration/out/core_IF_test_gmap.mat');
cmap = load('/Users/Andrea/Repositories/goGPS_MATLAB/data/project/Antenna_Calibration/out/core_IF_test_cmap.mat');
g1map = load('/Users/Andrea/Repositories/goGPS_MATLAB/data/project/Antenna_Calibration/out/core_IF_test_g1map.mat');
c1map = load('/Users/Andrea/Repositories/goGPS_MATLAB/data/project/Antenna_Calibration/out/core_IF_test_c1map.mat');
Core_Utils.playAlert
%% Load cores
none = load('/Users/Andrea/Repositories/goGPS_MATLAB/data/project/Antenna_Calibration/out/test_uff/core_IF_test_none.mat');
zmap = load('/Users/Andrea/Repositories/goGPS_MATLAB/data/project/Antenna_Calibration/out/test_uff/core_IF_test_zmap.mat');
rmap = load('/Users/Andrea/Repositories/goGPS_MATLAB/data/project/Antenna_Calibration/out/test_uff/core_IF_test_rmap.mat');
gmap = load('/Users/Andrea/Repositories/goGPS_MATLAB/data/project/Antenna_Calibration/out/test_uff/core_IF_test_gmap.mat');
cmap = load('/Users/Andrea/Repositories/goGPS_MATLAB/data/project/Antenna_Calibration/out/test_uff/core_IF_test_cmap.mat');
g1map = load('/Users/Andrea/Repositories/goGPS_MATLAB/data/project/Antenna_Calibration/out/test_uff/core_IF_test_g1map.mat');
c1map = load('/Users/Andrea/Repositories/goGPS_MATLAB/data/project/Antenna_Calibration/out/test_uff/core_IF_test_c1map.mat');

%% Troposphere
std_gr = [];
[zwd, p_time, id_sync, tge, tgn] = none.core.rec.getZtd_mr();
std_gr(1, :) = [mean(std(zwd - mean(zwd,2,'omitnan'),'omitnan')) mean(std(tgn - mean(tgn,2,'omitnan'),'omitnan')) mean(std(tge - mean(tge,2,'omitnan'),'omitnan'))];

[zwd, p_time, id_sync, tge, tgn] = zmap.core.rec.getZtd_mr();
std_gr(2, :) = [mean(std(zwd - mean(zwd,2,'omitnan'),'omitnan')) mean(std(tgn - mean(tgn,2,'omitnan'),'omitnan')) mean(std(tge - mean(tge,2,'omitnan'),'omitnan'))];

[zwd, p_time, id_sync, tge, tgn] = rmap.core.rec.getZtd_mr();
std_gr(3, :) = [mean(std(zwd - mean(zwd,2,'omitnan'),'omitnan')) mean(std(tgn - mean(tgn,2,'omitnan'),'omitnan')) mean(std(tge - mean(tge,2,'omitnan'),'omitnan'))];

[zwd, p_time, id_sync, tge, tgn] = gmap.core.rec.getZtd_mr();
std_gr(4, :) = [mean(std(zwd - mean(zwd,2,'omitnan'),'omitnan')) mean(std(tgn - mean(tgn,2,'omitnan'),'omitnan')) mean(std(tge - mean(tge,2,'omitnan'),'omitnan'))];

[zwd, p_time, id_sync, tge, tgn] = cmap.core.rec.getZtd_mr();
std_gr(5, :) = [mean(std(zwd - mean(zwd,2,'omitnan'),'omitnan')) mean(std(tgn - mean(tgn,2,'omitnan'),'omitnan')) mean(std(tge - mean(tge,2,'omitnan'),'omitnan'))];

[zwd, p_time, id_sync, tge, tgn] = g1map.core.rec.getZtd_mr();
std_gr(6, :) = [mean(std(zwd - mean(zwd,2,'omitnan'),'omitnan')) mean(std(tgn - mean(tgn,2,'omitnan'),'omitnan')) mean(std(tge - mean(tge,2,'omitnan'),'omitnan'))];

%[zwd, p_time, id_sync, tge, tgn] = c1map.core.rec.getZtd_mr();
%std_gr(7, :) = [mean(std(zwd - mean(zwd,2,'omitnan'),'omitnan')) mean(std(tgn - mean(tgn,2,'omitnan'),'omitnan')) mean(std(tge - mean(tge,2,'omitnan'),'omitnan'))];

figure; plot(bsxfun(@rdivide, std_gr, std_gr(1,:)) * 100, 'o-', 'LineWidth', 2);

%% stds
figure;

core_list = [none.core zmap.core rmap.core gmap.core cmap.core g1map.core c1map.core];
std_core = [];
for c = 2:numel(core_list)
    tmp_core = core_list(c);
    std_rec = [];
    for r = 1:6
        std_rec(r) = mean(std(zero2nan(tmp_core.rec(r).out.sat.res.value), 'omitnan'), 'omitnan');
    end
    std_core(c) = mean(std_rec(r));
end
figure; plot(bsxfun(@rdivide, std_core, std_core(1)) * 100, 'o-', 'LineWidth', 2);

%% Coordinates
for i=1
    rec_set = 1:5;
    std_enu = [];
    enu = [];
    for r = rec_set
        enu(:,:,r) = none.core.rec(r).out.add_coo(1).coo.getENU;
        enu(:,:,r) = enu(:,:,r) - median(enu(:,:,r), 1, 'omitnan');
    end
    enu = permute(enu, [1 3 2]);
    std_enu(1, :) = [mean(std(enu(:,:,1), 1, 'omitnan')) mean(std(enu(:,:,2), 1, 'omitnan')) mean(std(enu(:,:,3), 1, 'omitnan'))];
    
    enu = [];
    for r = rec_set
        enu(:,:,r) = zmap.core.rec(r).out.add_coo(1).coo.getENU;
        enu(:,:,r) = enu(:,:,r) - median(enu(:,:,r), 1, 'omitnan');
    end
    enu = permute(enu, [1 3 2]);
    std_enu(2, :) = [mean(std(enu(:,:,1), 1, 'omitnan')) mean(std(enu(:,:,2), 1, 'omitnan')) mean(std(enu(:,:,3), 1, 'omitnan'))];
    
    enu = [];
    for r = rec_set
        enu(:,:,r) = rmap.core.rec(r).out.add_coo(1).coo.getENU;
        enu(:,:,r) = enu(:,:,r) - median(enu(:,:,r), 1, 'omitnan');
    end
    enu = permute(enu, [1 3 2]);
    std_enu(3, :) = [mean(std(enu(:,:,1), 1, 'omitnan')) mean(std(enu(:,:,2), 1, 'omitnan')) mean(std(enu(:,:,3), 1, 'omitnan'))];
    
    enu = [];
    for r = rec_set
        enu(:,:,r) = gmap.core.rec(r).out.add_coo(1).coo.getENU;
        enu(:,:,r) = enu(:,:,r) - median(enu(:,:,r), 1, 'omitnan');
    end
    enu = permute(enu, [1 3 2]);
    std_enu(4, :) = [mean(std(enu(:,:,1), 1, 'omitnan')) mean(std(enu(:,:,2), 1, 'omitnan')) mean(std(enu(:,:,3), 1, 'omitnan'))];
    
    enu = [];
    for r = rec_set
        enu(:,:,r) = cmap.core.rec(r).out.add_coo(1).coo.getENU;
        enu(:,:,r) = enu(:,:,r) - median(enu(:,:,r), 1, 'omitnan');
    end
    enu = permute(enu, [1 3 2]);
    std_enu(5, :) = [mean(std(enu(:,:,1), 1, 'omitnan')) mean(std(enu(:,:,2), 1, 'omitnan')) mean(std(enu(:,:,3), 1, 'omitnan'))];
    
    enu = [];
    for r = rec_set
        enu(:,:,r) = g1map.core.rec(r).out.add_coo(1).coo.getENU;
        enu(:,:,r) = enu(:,:,r) - median(enu(:,:,r), 1, 'omitnan');
    end
    enu = permute(enu, [1 3 2]);
    std_enu(6, :) = [mean(std(enu(:,:,1), 1, 'omitnan')) mean(std(enu(:,:,2), 1, 'omitnan')) mean(std(enu(:,:,3), 1, 'omitnan'))];
    
    enu = [];
    for r = rec_set
        enu(:,:,r) = c1map.core.rec(r).out.add_coo(1).coo.getENU;
        enu(:,:,r) = enu(:,:,r) - median(enu(:,:,r), 1, 'omitnan');
    end
    enu = permute(enu, [1 3 2]);
    std_enu(7, :) = [mean(std(enu(:,:,1), 1, 'omitnan')) mean(std(enu(:,:,2), 1, 'omitnan')) mean(std(enu(:,:,3), 1, 'omitnan'))];
    
    figure; plot(bsxfun(@rdivide, std_enu, std_enu(1,:)) * 100, 'o-', 'LineWidth', 2);
    if numel(rec_set) == 1 
        title(sprintf('Coordinates stability improvement of %s', rec(rec_set).getMarkerName4Ch));
    else
        title(sprintf('Coordinates stability improvement'));
    end
end
%%
coo = [none.core.rec.get('WTZA').out.add_coo(1).coo ...
    zmap.core.rec.get('WTZA').out.add_coo(1).coo ...
    rmap.core.rec.get('WTZA').out.add_coo(1).coo ...
    g1map.core.rec.get('WTZA').out.add_coo(1).coo ...
    c1map.core.rec.get('WTZA').out.add_coo(1).coo];
    coo.showPositionENU
