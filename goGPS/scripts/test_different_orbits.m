tic;
CENTERS = {'igs_gps'; 'emr'; 'esoc'; 'code'; 'code1d'; 'code5s'; 'gfz'; 'cnes'; 'wuhan'};

state = Main_Settings('last_settings.ini');
sta_co = GNSS_Station(); sta_co(1) = [];
sta_u2 = GNSS_Station(); sta_u2(1) = [];
for i = 1 : numel(CENTERS)
    state.selected_center = CENTERS{i};
    goGPS(state, 0);
    [res_zwd{i}, time{i}] = rec(1).getZwd();
    sta_u2(i) = rec(1);
    sta_co(i) = rec(2);
end
toc;

sta = sta_co;
fh = sta.showZtd;
ax = fh.Children(end);
axes(ax);
ylim([212.5 220]);
title(sprintf('%s of %s station \nusing different GPS orbits providers - combined engine\\fontsize{5} \n', ax.Title.String, rec(1).getMarkerName4Ch));
legend(CENTERS, 'Interpreter', 'none', 'Location', 'NorthWest');
fh = gcf; Core_UI.addBeautifyMenu(fh); Core_UI.beautifyFig(fh, 'light');

sta = sta_u2;
fh = sta.showZtd;
ax = fh.Children(end);
axes(ax);
ylim([212.5 220]);
title(sprintf('%s of %s station \nusing different GPS orbits providers - uncombined engine\\fontsize{5} \n', ax.Title.String, rec(1).getMarkerName4Ch));
legend(CENTERS, 'Interpreter', 'none', 'Location', 'NorthWest');
fh = gcf; Core_UI.addBeautifyMenu(fh); Core_UI.beautifyFig(fh, 'light');

%%
%sta_co_ca = sta_co;
%sta_u2_ca = sta_u2;
