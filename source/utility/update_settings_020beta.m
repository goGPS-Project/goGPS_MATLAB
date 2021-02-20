function update_settings_020beta(settings_folder)

if (nargin == 0)
    settings_folder = '../data/settings';
end

out = zeros(21,1);

out(1)  = update_settings(settings_folder, 'u_com_detect', -1);
out(2)  = update_settings(settings_folder, 'rinex_out', -1);
out(3)  = update_settings(settings_folder, 'data_streams', -1);
out(4)  = update_settings(settings_folder, 'server_delay', -1);
out(5)  = update_settings(settings_folder, 'plotproc', 1);
out(6)  = update_settings(settings_folder, 'stopGOstop', 0);
out(7)  = update_settings(settings_folder, 'diff_amb_estim', -1);
out(8)  = update_settings(settings_folder, 'LS_amb_estim', -1);
out(9)  = update_settings(settings_folder, 'RINEX_rover_nav', -1);
out(10) = update_settings(settings_folder, 'RINEX_master_nav', 'RINEX_nav');
out(11) = update_settings(settings_folder, 'flag_doppler', 0);
out(12) = update_settings(settings_folder, 'amb_select', 2);
out(13) = update_settings(settings_folder, 'num_receivers', 1);
out(14) = update_settings(settings_folder, 'com_select_0', 1);
out(15) = update_settings(settings_folder, 'com_select_1', 1);
out(16) = update_settings(settings_folder, 'com_select_2', 1);
out(17) = update_settings(settings_folder, 'com_select_3', 1);
out(18) = update_settings(settings_folder, 'protocol_select_0', 1);
out(19) = update_settings(settings_folder, 'protocol_select_1', 1);
out(20) = update_settings(settings_folder, 'protocol_select_2', 1);
out(21) = update_settings(settings_folder, 'protocol_select_3', 1);

if (sum(out) == length(out))
    disp('Done!');
end
