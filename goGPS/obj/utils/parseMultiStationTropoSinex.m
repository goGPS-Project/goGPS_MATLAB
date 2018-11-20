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
    if strfind(txt(lim(st_idx,1): lim(st_idx,2)),'TROP/DESCRIPTION')
        for l = st_idx : end_idx
            if strfind(txt(lim(l,1): lim(l,2)),' SOLUTION_FIELDS_1')
                pars = strsplit(strtrim(txt(lim(l,1): lim(l,2))));
                pars(1) = [];
                n_par = length(pars);
                for i = 1: n_par
                    pars{i} = strrep(pars{i},'#','NUM_');
                    all_res.(pars{i}) = [];
                end
            end
        end
        
    elseif strfind(txt(lim(st_idx,1): lim(st_idx,2)),'TROP/STA_COORDINATES')
        n_lin = end_idx - st_idx -2;
        sta_4char = txt(repmat(lim(st_idx+2:end_idx-1,1),1,4) + repmat([1:4],n_lin ,1));
        xyz = reshape(sscanf( (txt(repmat(lim(st_idx+2:end_idx-1,1),1,38) + repmat([17:54], n_lin,1)))','%f %f %f\n'),3,n_lin)';
    elseif strfind(txt(lim(st_idx,1): lim(st_idx,2)),'TROP/SOLUTION')
        n_lin = end_idx - st_idx -2;
        sta_4char_trp = txt(repmat(lim(st_idx+2:end_idx-1,1),1,4) + repmat([1:4],n_lin,1));
        obs_lines = st_idx+2:end_idx-1;
        idx_no_ast = txt(lim(obs_lines,1)) ~= '*';
        idx_no_ast = obs_lines(idx_no_ast);
         n_lin = length(idx_no_ast);
        sta_4char_trp = txt(repmat(lim(idx_no_ast,1),1,4) + repmat([1:4],n_lin,1));
        end_col = min(lim(idx_no_ast,3));
       
        data = reshape(sscanf( (txt(repmat(lim(idx_no_ast,1),1,end_col) + repmat([0 :(end_col-1)],n_lin,1)))',['%*s %f:%f:%f' repmat(' %f',1,n_par)] ),n_par+3,n_lin)';
        year = data(:,1);
        idx_70 = year<70;
        year(idx_70) = year(idx_70) + 2000;
        year(~idx_70) = year(~idx_70) + 1900;
        doy = data(:,2);
        sod = data(:,3);
        for i = 1 : n_par
            all_res.(pars{i}) = [all_res.(pars{i}) ;  data(:,3+i)];
        end
    end
           
end
for s = 1: size(sta_4char,1)
    c_sta_4char = sta_4char(s,:);
    idx_sta  = strLineMatch(sta_4char_trp,c_sta_4char);
    if sum(idx_sta) > 0
        for i = 1 : n_par
        results.(c_sta_4char).(pars{i}) = all_res.(pars{i})(idx_sta);
        end
        results.(c_sta_4char).time = GPS_Time.fromDoySod(year(idx_sta),doy(idx_sta),sod(idx_sta));
        results.(c_sta_4char).xyz = xyz(s,:);
    end
end
end
