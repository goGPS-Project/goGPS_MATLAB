function createEmptyProject(base_dir, prj_name)
% create empty config file
fnp = File_Name_Processor();
[status, msg, msgID] = mkdir(fnp.checkPath([base_dir '/' prj_name]));
[status, msg, msgID] = mkdir(fnp.checkPath([base_dir '/' prj_name '/config']));
[status, msg, msgID] = mkdir(fnp.checkPath([base_dir '/' prj_name '/out']));
[status, msg, msgID] = mkdir(fnp.checkPath([base_dir '/' prj_name '/RINEX']));
[status, msg, msgID] = mkdir(fnp.checkPath([base_dir '/' prj_name '/station']));
[status, msg, msgID] = mkdir(fnp.checkPath([base_dir '/' prj_name '/station/CRD']));
[status, msg, msgID] = mkdir(fnp.checkPath([base_dir '/' prj_name '/station/MET']));
state = Main_Settings();
state.setPrjHome(fnp.checkPath([base_dir '/' prj_name]));
state.prj_name = prj_name;
fid = fopen(fnp.checkPath([base_dir '/' prj_name '/config/config.ini']),'w');
fprintf(fid,'%s',state.toString());
fclose(fid);
end