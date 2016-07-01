function [diff_L4, diff_P4, commontime, stations_idx, L1_series, L2_series, L4_series, P4_series] = compute_diffL4(L1_series, L2_series, P1_series, P2_series, name_series, time_series)

lambda1 = 0.190293672798365;
lambda2 = 0.244210213424568;

%compute dL4
n_sta = length(name_series);

%choose station for interpolation 
commontime = time_series{1}; 
for i = 2 : n_sta
    commontime = intersect(commontime,time_series{i});
end

%matrices same length as L4
stations_idx = zeros(n_sta, length(commontime));
L4_series = zeros(32, length(commontime), n_sta);
P4_series = zeros(32, length(commontime), n_sta);

%matrices same size as diff_L4
diff_L4 = zeros(32, length(commontime) - 1, n_sta); 
diff_P4 = zeros(32, length(commontime) - 1, n_sta); 

%compute L4 = L1 - L2, dL4, P4 = P1 - P2, dP4
for sta = 1 : n_sta
    [~, ~, stations_idx(sta,:)] = intersect(commontime, time_series{sta});
    
    for PRN = 1 : 32
        L4_series(PRN,:,sta)=L1_series{sta}(PRN,stations_idx(sta,:))*lambda1 - L2_series{sta}(PRN,stations_idx(sta,:))*lambda2;
        P4_series(PRN,:,sta)=P2_series{sta}(PRN,stations_idx(sta,:))         - P1_series{sta}(PRN,stations_idx(sta,:));
        diff_L4(PRN,:,sta)=diff(L4_series(PRN,:,sta));
        diff_P4(PRN,:,sta)=diff(P4_series(PRN,:,sta));
    end
end
