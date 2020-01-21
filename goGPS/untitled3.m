res_tmp_s = Receiver_Commons.smoothMat(res_tmp, 'spline', 10);
step = 0.5;
[gData, wGrid] = simpleGridder.go(this.sat.el(id_ok).*deg2rad, (this.sat.az(id_ok)).*deg2rad, scale * res_tmp_s(id_ok), step); 

[el_grid, ~] = getGrid(step(end));
[~, az_grid] = getGrid(step(1));
[az, el] = meshgrid(az_grid, el_grid);
id_okg = find(~isnan(gData));
[~, id_sort] = sort(abs(gData(id_okg)));
id_okg = id_okg(id_sort);
figure; polarScatter(az(id_okg)*deg2rad, (90-el(id_okg)).*deg2rad, 20, gData(id_okg), 'filled');
dockAllFigures
%figure; imagesc(nan2zero(gData(1:size(gData,1)/2,:)));  
colormap((Cmap.get('PuOr', 2^11))); colorbar; caxis([-4  4])
fh = gcf; Core_UI.addBeautifyMenu(fh); Core_UI.beautifyFig(fh, 'dark');
%%
figure; 
gain = [];
for r=1:6
    std_noz = std(serialize(noz.core.rec(r).work.getResPhases([],'1')), 'omitnan');
    std_z = std(serialize(core.rec(r).work.getResPhases([],'1')), 'omitnan');
    plot(r, std_noz, '.', 'Color', [0.4 0.4 0.4], 'MarkerSize', 20); hold on
    plot(r, std_z, '*', 'Color', Core_UI.getColor(r,6), 'MarkerSize', 10);
    gain(r) = 100 - (std_z * 100) / std_noz;
end
%%

sat_id = 1;
res_tmp = noz.core.rec(1).work.getResPhases('G', '1');
figure; plotSep(mod((1:size(res_tmp,1))*30, 86150), Receiver_Commons.smoothMat(res_tmp(:,sat_id) * 1e3, 'spline', 10),  'LineWidth', 2)
%figure; plotSep(mod((1:size(res_tmp,1))*30, 86150), res_tmp(:,sat_id) * 1e3,  'LineWidth', 1)
res_tmp = core.rec(1).work.getResPhases('G', '1');
hold on; plotSep(mod((1:size(res_tmp,1))*30, 86150), Receiver_Commons.smoothMat(res_tmp(:,sat_id) * 1e3, 'spline', 10), 'LineWidth', 2)
%hold on; plotSep(mod((1:size(res_tmp,1))*30, 86150), res_tmp(:,sat_id) * 1e3, 'LineWidth', 1)
