dt_archive = [];

%% BRDC ================================================================================================

ini_file = '/Volumes/Data/goGPS_data/project/JRC_timing_test/config/config_1s.ini';
ini = Main_Settings(ini_file);

% igs_glo, igs_gps, code, code_mgex, gfz, jaxa
ini.selected_center = 'igs_gps';
% 'final' 'rapid' 'ultra' 'broadcast'
ini.preferred_eph = 'broadcast';
%'final' 'predicted1' 'predicted2' 'broadcast'
ini.preferred_iono = 'predicted2';

goGPS(ini, false);

work = rec.work;

time = work.time.getEpoch(work.getIdSync).getMatlabTime();
nans = zero2nan(double(~work.getMissingEpochs()));
prepro_dt = work.getDtPrePro .* nans(work.getIdSync);
ppp_dt = work.getTotalDt .* nans(work.getIdSync);

dt_archive.igs_brdc = struct('time', time, 'prepro_dt', prepro_dt, 'ppp_dt', ppp_dt);

f = figure; plot(dt_archive.igs_brdc.time, [0; diff(dt_archive.igs_brdc.prepro_dt)] * Core_Utils.V_LIGHT,'.')
hold on; plot(dt_archive.igs_brdc.time, [0; diff(dt_archive.igs_brdc.ppp_dt)] * Core_Utils.V_LIGHT, '.');
axis tight; ylim([-0.5 0.5]); setTimeTicks(); grid on; ylabel('Clock error [m]')
setAllLinesWidth(2)
Core_UI.beautifyFig(f);
title('BRDC (DIFF)');

%% REALTIME ============================================================================================

ini_file = '/Volumes/Data/goGPS_data/project/JRC_timing_test/config/config_1s.ini';
ini = Main_Settings(ini_file);

% igs_gps, code, code_mgex, gfz, jaxa
ini.selected_center = 'bnc';
ini.eph_dir = '/Volumes/Data/goGPS_data/satellite/EPH/BNC_IGS03';
ini.clk_dir = '/Volumes/Data/goGPS_data/satellite/CLK/BNC_IGS03';
% 'final' 'rapid' 'ultra' 'broadcast'
ini.preferred_eph = 'final';
% 'final' 'predicted1' 'predicted2' 'broadcast'
ini.preferred_iono = 'predicted2';

goGPS(ini, false);

work = rec.work;

time = work.time.getEpoch(work.getIdSync).getMatlabTime();
nans = zero2nan(double(~work.getMissingEpochs()));
prepro_dt = work.getDtPrePro .* nans(work.getIdSync);
ppp_dt = work.getTotalDt .* nans(work.getIdSync);

dt_archive.igs_realtime = struct('time', time, 'prepro_dt', prepro_dt, 'ppp_dt', ppp_dt);

f = figure; plot(dt_archive.igs_realtime.time, [0; diff(dt_archive.igs_realtime.prepro_dt)] * Core_Utils.V_LIGHT,'.')
hold on; plot(dt_archive.igs_realtime.time, [0; diff(dt_archive.igs_realtime.ppp_dt)] * Core_Utils.V_LIGHT, '.');
axis tight; ylim([-0.5 0.5]); setTimeTicks(); grid on; ylabel('Clock error [m]')
setAllLinesWidth(2)
Core_UI.beautifyFig(f);
title('REALTIME (DIFF)');

%% ULTRA ===============================================================================================

ini_file = '/Volumes/Data/goGPS_data/project/JRC_timing_test/config/config_1s.ini';
ini = Main_Settings(ini_file);

% igs_glo, igs_gps, code, code_mgex, gfz, jaxa
ini.selected_center = 'igs_gps';
% 'final' 'rapid' 'ultra' 'broadcast'
ini.preferred_eph = 'ultra';
% 'final' 'predicted1' 'predicted2' 'broadcast'
ini.preferred_iono = 'final';

goGPS(ini, false);

work = rec.work;

time = work.time.getEpoch(work.getIdSync).getMatlabTime();
nans = zero2nan(double(~work.getMissingEpochs()));
prepro_dt = work.getDtPrePro .* nans(work.getIdSync);
ppp_dt = work.getTotalDt .* nans(work.getIdSync);

dt_archive.igs_ultra = struct('time', time, 'prepro_dt', prepro_dt, 'ppp_dt', ppp_dt);

f = figure; plot(dt_archive.igs_ultra.time, [0; diff(dt_archive.igs_ultra.prepro_dt)] * Core_Utils.V_LIGHT,'.')
hold on; plot(dt_archive.igs_ultra.time, [0; diff(dt_archive.igs_ultra.ppp_dt)] * Core_Utils.V_LIGHT, '.');
axis tight; ylim([-0.5 0.5]); setTimeTicks(); grid on; ylabel('Clock error [m]')
setAllLinesWidth(2)
Core_UI.beautifyFig(f);
title('ULTRA (DIFF)');

%% RAPID ===============================================================================================

ini_file = '/Volumes/Data/goGPS_data/project/JRC_timing_test/config/config_1s.ini';
ini = Main_Settings(ini_file);

% igs_glo, igs_gps, code, code_mgex, gfz, jaxa
ini.selected_center = 'igs_gps';
% 'final' 'rapid' 'ultra' 'broadcast'
ini.preferred_eph = 'rapid';
% 'final' 'predicted1' 'predicted2' 'broadcast'
ini.preferred_iono = 'final';

goGPS(ini, false);

work = rec.work;

time = work.time.getEpoch(work.getIdSync).getMatlabTime();
nans = zero2nan(double(~work.getMissingEpochs()));
prepro_dt = work.getDtPrePro .* nans(work.getIdSync);
ppp_dt = work.getTotalDt .* nans(work.getIdSync);

dt_archive.igs_rapid = struct('time', time, 'prepro_dt', prepro_dt, 'ppp_dt', ppp_dt);

f = figure; plot(dt_archive.igs_rapid.time, [0; diff(dt_archive.igs_rapid.prepro_dt)] * Core_Utils.V_LIGHT,'.')
hold on; plot(dt_archive.igs_rapid.time, [0; diff(dt_archive.igs_rapid.ppp_dt)] * Core_Utils.V_LIGHT, '.');
axis tight; ylim([-0.5 0.5]); setTimeTicks(); grid on; ylabel('Clock error [m]')
setAllLinesWidth(2)
Core_UI.beautifyFig(f);
title('RAPID (DIFF)');

%% FINAL ===============================================================================================

ini_file = '/Volumes/Data/goGPS_data/project/JRC_timing_test/config/config_1s.ini';
ini = Main_Settings(ini_file);

% igs_gps, code, code_mgex, gfz, jaxa
ini.selected_center = 'igs_gps';
% 'final' 'rapid' 'ultra' 'broadcast'
ini.preferred_eph = 'final';
% 'final' 'predicted1' 'predicted2' 'broadcast'
ini.preferred_iono = 'final';

goGPS(ini, false);

work = rec.work;

time = work.time.getEpoch(work.getIdSync).getMatlabTime();
nans = zero2nan(double(~work.getMissingEpochs()));
prepro_dt = work.getDtPrePro .* nans(work.getIdSync);
ppp_dt = work.getTotalDt .* nans(work.getIdSync);

dt_archive.igs_final = struct('time', time, 'prepro_dt', prepro_dt, 'ppp_dt', ppp_dt);

f = figure; plot(dt_archive.igs_final.time, [0; diff(dt_archive.igs_final.prepro_dt)] * Core_Utils.V_LIGHT,'.')
hold on; plot(dt_archive.igs_final.time, [0; diff(dt_archive.igs_final.ppp_dt)] * Core_Utils.V_LIGHT, '.');
axis tight; ylim([-0.5 0.5]); setTimeTicks(); grid on; ylabel('Clock error [m]')
setAllLinesWidth(2)
Core_UI.beautifyFig(f);
title('FINAL (DIFF)');

%% ALL DIFF PLOT =======================================================================================

f = figure; plot(dt_archive.igs_brdc.time, [0; diff(dt_archive.igs_brdc.prepro_dt)] * Core_Utils.V_LIGHT,'.')
hold on; plot(dt_archive.igs_brdc.time, [0; diff(dt_archive.igs_brdc.ppp_dt)] * Core_Utils.V_LIGHT, '.');
axis tight; ylim([-0.5 0.5]); axis tight; ylim([-0.5 0.5]); setTimeTicks(); grid on; ylabel('Clock error [m]')
Core_UI.beautifyFig(f);
drawnow;
title('BRDC (DIFF)');

f = figure; plot(dt_archive.igs_realtime.time, [0; diff(dt_archive.igs_realtime.prepro_dt)] * Core_Utils.V_LIGHT,'.')
hold on; plot(dt_archive.igs_realtime.time, [0; diff(dt_archive.igs_realtime.ppp_dt)] * Core_Utils.V_LIGHT, '.');
axis tight; ylim([-0.5 0.5]); setTimeTicks(); grid on; ylabel('Clock error [m]')
Core_UI.beautifyFig(f);
drawnow;
title('REALTIME (DIFF)');

f = figure; plot(dt_archive.igs_ultra.time, [0; diff(dt_archive.igs_ultra.prepro_dt)] * Core_Utils.V_LIGHT,'.')
hold on; plot(dt_archive.igs_ultra.time, [0; diff(dt_archive.igs_ultra.ppp_dt)] * Core_Utils.V_LIGHT, '.');
axis tight; ylim([-0.5 0.5]); setTimeTicks(); grid on; ylabel('Clock error [m]')
Core_UI.beautifyFig(f);
drawnow;
title('ULTRA (DIFF)');

f = figure; plot(dt_archive.igs_rapid.time, [0; diff(dt_archive.igs_rapid.prepro_dt)] * Core_Utils.V_LIGHT,'.')
hold on; plot(dt_archive.igs_rapid.time, [0; diff(dt_archive.igs_rapid.ppp_dt)] * Core_Utils.V_LIGHT, '.');
axis tight; ylim([-0.5 0.5]); setTimeTicks(); grid on; ylabel('Clock error [m]')
Core_UI.beautifyFig(f);
drawnow;
title('RAPID (DIFF)');

f = figure; plot(dt_archive.igs_final.time, [0; diff(dt_archive.igs_final.prepro_dt)] * Core_Utils.V_LIGHT,'.')
hold on; plot(dt_archive.igs_final.time, [0; diff(dt_archive.igs_final.ppp_dt)] * Core_Utils.V_LIGHT, '.');
axis tight; ylim([-0.5 0.5]); setTimeTicks(); grid on; ylabel('Clock error [m]')
Core_UI.beautifyFig(f);
drawnow;
title('FINAL (DIFF)');

dockAllFigures
%%

%% ALL DIFF PLOT =======================================================================================

f = figure; plot(dt_archive.igs_brdc.time, dt_archive.igs_brdc.prepro_dt * Core_Utils.V_LIGHT,'.')
hold on; plot(dt_archive.igs_brdc.time, dt_archive.igs_brdc.ppp_dt * Core_Utils.V_LIGHT, '.');
axis tight; setTimeTicks(); grid on; ylabel('Clock error [m]')
Core_UI.beautifyFig(f);
drawnow;
title('BRDC');


f = figure; plot(dt_archive.igs_realtime.time, dt_archive.igs_realtime.prepro_dt * Core_Utils.V_LIGHT,'.')
hold on; plot(dt_archive.igs_realtime.time, dt_archive.igs_realtime.ppp_dt * Core_Utils.V_LIGHT, '.');
axis tight; setTimeTicks(); grid on; ylabel('Clock error [m]')
Core_UI.beautifyFig(f);
drawnow;
title('REALTIME');

f = figure; plot(dt_archive.igs_ultra.time, dt_archive.igs_ultra.prepro_dt * Core_Utils.V_LIGHT,'.')
hold on; plot(dt_archive.igs_ultra.time, dt_archive.igs_ultra.ppp_dt * Core_Utils.V_LIGHT, '.');
axis tight; setTimeTicks(); grid on; ylabel('Clock error [m]')
Core_UI.beautifyFig(f);
drawnow;
title('ULTRA');

f = figure; plot(dt_archive.igs_rapid.time, dt_archive.igs_rapid.prepro_dt * Core_Utils.V_LIGHT,'.')
hold on; plot(dt_archive.igs_rapid.time, dt_archive.igs_rapid.ppp_dt * Core_Utils.V_LIGHT, '.');
axis tight; setTimeTicks(); grid on; ylabel('Clock error [m]')
Core_UI.beautifyFig(f);
drawnow;
title('RAPID');

f = figure; plot(dt_archive.igs_final.time, dt_archive.igs_final.prepro_dt * Core_Utils.V_LIGHT,'.')
hold on; plot(dt_archive.igs_final.time, dt_archive.igs_final.ppp_dt * Core_Utils.V_LIGHT, '.');
axis tight; setTimeTicks(); grid on; ylabel('Clock error [m]')
Core_UI.beautifyFig(f);
drawnow;
title('FINAL');

dockAllFigures
