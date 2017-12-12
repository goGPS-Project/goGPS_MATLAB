[ph, wl, id_ph1] = rec(1).getPhases();
rec(1).updateErrTropo('all',1)
rec(1).updateSolidEarthCorr();
sensor1 = Core_Pre_Processing.diffAndPred(zero2nan( zero2nan(rec(1).getSyntPhObs)' )); %sensor1 = bsxfun(@minus, sensor1, median(sensor1, 2, 'omitnan'));%- zero2nan(rec(1).getSyntPhObs)'
rec(1).updateAzimuthElevation()
sensor1 = Core_Pre_Processing.diffAndPred(sensor1);
sensor1 = Core_Pre_Processing.diffAndPred(sensor1);
% [ph, wl, id_ph2] = rec(2).getPhases;
% sensor2 = Core_Pre_Processing.diffAndPred(ph - zero2nan(rec(2).getSyntPhObs)'); sensor2 = bsxfun(@minus, sensor2, median(sensor2, 2, ‘omitnan’)); figure; plot(sensor2)
% rec(2).updateAzimuthElevation()

% [ph, wl, id_ph3] = rec(3).getPhases;
% sensor3 = Core_Pre_Processing.diffAndPred(ph - zero2nan(rec(3).getSyntPhObs)'); sensor3 = bsxfun(@minus, sensor3, median(sensor3, 2, ‘omitnan’)); figure; plot(sensor3)
% rec(3).updateAzimuthElevation()

id_ok = (~isnan(sensor1));
az = rec(1).sat.az(:,rec(1).go_id(id_ph1));
el = rec(1).sat.el(:,rec(1).go_id(id_ph1));
flag = abs(sensor1(id_ok)) > 1;%flagExpand(, 3);
%figure; clf; scatplot(serialize(az(id_ok))/180*pi,serialize(90-el(id_ok))/180*pi, 50, flag, 'filled');
az = az(id_ok);
el = el(id_ok);
figure; clf; scatter(serialize(az),serialize(el), 30, flag, 'filled');
figure; clf; scatter(serialize(az(flag)),serialize(el(flag)), 10, 'k', 'filled');