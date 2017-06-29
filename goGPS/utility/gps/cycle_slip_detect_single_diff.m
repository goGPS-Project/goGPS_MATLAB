function [ph_R, ph_M] = cycle_slip_detect_single_diff(ph_R, ph_M, interval)

for s = 1 : size(ph_R,1)
    if (any(ph_R(s,:)) && any(ph_M(s,:)))
        idxR = find(ph_R(s,:) ~= 0);
        idxM = find(ph_M(s,:) ~= 0);
        idx = intersect(idxR,idxM);
        ph_R_tmp = ph_R(s,idx);
        ph_M_tmp = ph_M(s,idx);
        ph_R_tmp(1,ph_R_tmp(1,:)==0) = NaN;
        ph_M_tmp(1,ph_R_tmp(1,:)==0) = NaN;
        delta_sd = diff(ph_R_tmp(1,:) - ph_M_tmp(1,:))';
        delta_sd(delta_sd == 0) = NaN;
        not_zero = find(delta_sd ~= 0);
        not_nan  = find(~isnan(delta_sd));
        avail_sd = intersect(not_zero, not_nan);
        
        if (~isempty(delta_sd(avail_sd)))
            outliers = batch_outlier_detection(delta_sd(avail_sd),median(round(interval)));
            [~,jmp_sd] = intersect(delta_sd,outliers);
        end
        jmp_sd = sort(jmp_sd);
        ph_R(s, idx(jmp_sd)) = 0;
        ph_M(s, idx(jmp_sd)) = 0;
    end
end
