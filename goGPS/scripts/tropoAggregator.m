%%% aggregator
clear
% min_num_epochs = 7900;
min_num_epochs = 0;
GNSS_ZTD = struct();
in_dir = 'V:\goGPS_data\project\South_Africa\out';
files = dir(in_dir);
% estimated relationheoght ztd
x = [];
y = [];
for i = 1:length(files)
     fname = files(i).name;
    if strfind(fname,'.mat')
        tropo_data = load([ in_dir '/' fname]);
        x = [x; tropo_data.h_ortho*ones(size(tropo_data.ztd(~isnan(tropo_data.ztd))))];
        y = [y; tropo_data.ztd(~isnan(tropo_data.ztd))];
    end
end
A = [ones(size(x)) x x.^2];
p = A\y;
idx_rem = y < 0 & y > 3;
A(idx_rem,:) = [];
y(idx_rem) = [];
% remove outlier
res = y - A*p;
res  = res / mean(abs(res));
idx_rem = abs(res) > 2.5;
A_no_out= A(~idx_rem,:);
y_no_out = y(~idx_rem);
p = A_no_out\y_no_out;
% remove outlier
res = y - A*p;
res  = res / mean(abs(res));
idx_rem = abs(res) > 2.5;
A_no_out= A(~idx_rem,:);
y_no_out = y(~idx_rem);
p = A_no_out\y_no_out;
 figure; 
 scatter(x,y);
 hold on
 plot(x,A*p);
for i = 1:length(files)
    fname = files(i).name;
    if strfind(fname,'.mat') & (length(fname) > 4)
        tropo_data = load([ in_dir '/' fname]);
        ztd_ok = p(1) + p(2) * tropo_data.h_ortho + p(3) * tropo_data.h_ortho.^2;
        d_ztd = diff(tropo_data.ztd)./diff(tropo_data.utc_time)*(30/86400);
        mv_avg = ones(1,20)/20;
        md_ztd = filter(mv_avg,1,d_ztd);
        if sum(abs(tropo_data.ztd - ztd_ok) > 0.3) > 1 || sum(isnan(tropo_data.ztd)) > 1 || length(tropo_data.ztd) < min_num_epochs || max(abs(md_ztd)) > 0.005
            %disp([' Outlier!! ' fname(1:4)])
            if max(abs(md_ztd)) > 0.005 & sum(abs(tropo_data.ztd - ztd_ok) > 0.3) == 0 || sum(abs(tropo_data.ztd - ztd_ok) > 0.3) > 1
                disp([' Outlier!! ' fname(1:4)])
            else
                disp([' Empty or partly empty!! ' fname(1:4)])
                disp( [sum(isnan(tropo_data.ztd)) length(tropo_data.ztd)])
            end
        else
            if isfield(GNSS_ZTD,fname(1:4))
                if abs(GNSS_ZTD.(fname(1:4)).h_ellips - tropo_data.h_ellips) > 0.1
                    sprintf('Warning dh = %f',(GNSS_ZTD.(fname(1:4)).h_ellips - tropo_data.h_ellips) );
                end
                GNSS_ZTD.(fname(1:4)).ztd = [GNSS_ZTD.(fname(1:4)).ztd ; tropo_data.ztd];
                GNSS_ZTD.(fname(1:4)).utc_time = [GNSS_ZTD.(fname(1:4)).utc_time ; tropo_data.utc_time];
            else
                GNSS_ZTD.(fname(1:4)) = tropo_data;
            end
        end
    end
end
