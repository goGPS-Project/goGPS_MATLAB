Core.load('/Volumes/Data/goGPS_data/project/GIMS_Potoska/out/core_20200630_161239.mat')
time_r = rec(1).out.time.first.getCopy();
time_1 = round( rec(1).out.time.getRefTime(time_r.getMatlabTime()) / rec(1).out.time.getRate);
time_2 = round( rec(7).out.time.getRefTime(time_r.getMatlabTime()) / rec(7).out.time.getRate);
[time_common,idx_a,idx_b] = intersect(time_1,time_2);
time_common_gps_1f = time_r.getCopy();
time_common_gps_1f.addSeconds(time_common*rec(7).out.time.getRate);
dt_1 = rec(1).out.dt + rec(1).out.dt_ip;
dt_2 = rec(7).out.dt + rec(7).out.dt_ip;
dt_1f = (dt_1(idx_a) - dt_2(idx_b));
Core.load('/Volumes/Data/goGPS_data/project/GIMS_Potoska/out/core_20200630_161408.mat')
time_r = rec(1).out.time.first.getCopy();
time_1 = round( rec(1).out.time.getRefTime(time_r.getMatlabTime()) / rec(1).out.time.getRate);
time_2 = round( rec(7).out.time.getRefTime(time_r.getMatlabTime()) / rec(7).out.time.getRate);
[time_common,idx_a,idx_b] = intersect(time_1,time_2);
time_common_gps_2f = time_r.getCopy();
time_common_gps_2f.addSeconds(time_common*rec(7).out.time.getRate);
dt_1 = rec(1).out.dt + rec(1).out.dt_ip;
dt_2 = rec(7).out.dt + rec(7).out.dt_ip;
dt_2f = (dt_1(idx_a) - dt_2(idx_b));



time_r = time_common_gps_1f.first.getCopy();
time_1 = round( time_common_gps_1f.getRefTime(time_r.getMatlabTime()) / rec(1).out.time.getRate);
time_2 = round( time_common_gps_2f.getRefTime(time_r.getMatlabTime()) / rec(7).out.time.getRate);
[time_common,idx_a,idx_b] = intersect(time_1,time_2);
time_common_diff = time_r.getCopy();
time_common_diff.addSeconds(time_common*rec(7).out.time.getRate);
figure; plot(time_common_diff.getMatlabTime,(dt_1f(idx_a)-dt_2f(idx_b))*Core_Utils.V_LIGHT)
setTimeTicks