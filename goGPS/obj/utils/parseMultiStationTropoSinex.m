function [results] = parseMultiStationTropoSinex(filename)
% parse a sinex fil containing tropopsheric products (e.g. the ones produced by epn)
fid = fopen(filename);
txt = fread(fid,'*char')';
fclose(fid);

% get new line separators
nl = regexp(txt, '\n')';
if nl(end) <  numel(txt)
    nl = [nl; numel(txt)];
end
lim = [[1; nl(1 : end - 1) + 1] (nl - 1)];
lim = [lim lim(:,2) - lim(:,1)];

% get starting block lines
st_idxes = find(txt(lim(:,1)) == '+');
end_idxes = find(txt(lim(:,1)) == '-');
for i = 1:length(st_idxes)
    st_idx = st_idxes(i);
    end_idx = end_idxes(i);
    if strfind(txt(lim(st_idx,1): lim(st_idx,2)),'TROP/STA_COORDINATES')
        n_lin = end_idx - st_idx -2;
        sta_4char = txt(repmat(lim(st_idx+2:end_idx-1,1),1,4) + repmat([1:4],n_lin ,1));
        xyz = reshape(sscanf( (txt(repmat(lim(st_idx+2:end_idx-1,1),1,38) + repmat([17:54], n_lin,1)))','%f %f %f\n'),3,n_lin)';
    elseif strfind(txt(lim(st_idx,1): lim(st_idx,2)),'TROP/SOLUTION')
        n_lin = end_idx - st_idx -2;
        sta_4char_trp = txt(repmat(lim(st_idx+2:end_idx-1,1),1,4) + repmat([1:4],n_lin,1));
        data = reshape(sscanf( (txt(repmat(lim(st_idx+2:end_idx-1,1),1,57) + repmat([6:62],n_lin,1)))','%f:%f:%f %f %f %f %f %f %f\n'),9,n_lin)';
        year = data(:,1);
        idx_70 = year<70;
        year(idx_70) = year(idx_70) + 2000;
        year(~idx_70) = year(~idx_70) + 1900;
        doy = data(:,2);
        sod = data(:,3);
        ztd = data(:,4)/1e3;
        ztd_std = data(:,5)/1e3;
        tgn = data(:,6)/1e3;
        tgn_std = data(:,7)/1e3;
        tge = data(:,8)/1e3;
        tge_std = data(:,9)/1e3;
    end
           
end
for s = 1: length(sta_4char)
    c_sta_4char = sta_4char(s,:);
    idx_sta  = idxCharLines(sta_4char_trp,c_sta_4char);
    if sum(idx_sta) > 0
        results.(c_sta_4char).ztd = ztd(idx_sta);
        results.(c_sta_4char).ztd_std = ztd_std(idx_sta);
        results.(c_sta_4char).tgn = tgn(idx_sta);
        results.(c_sta_4char).tgn_std = tgn_std(idx_sta);
        results.(c_sta_4char).tge = tge(idx_sta);
        results.(c_sta_4char).tge_std = tge_std(idx_sta);
        results.(c_sta_4char).time = GPS_Time.fromDoySod(year(idx_sta),doy(idx_sta),sod(idx_sta));
        results.(c_sta_4char).xyz = xyz(s,:);
    end
end
end
