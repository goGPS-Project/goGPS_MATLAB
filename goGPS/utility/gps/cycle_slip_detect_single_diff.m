function [ph_R, ph_M] = cycle_slip_detect_single_diff(ph_R, ph_M, interval)

ph_R_tmp = nan(size(ph_R));
ph_M_tmp = nan(size(ph_M));
delta_sd = nan(size(ph_R,1),size(ph_R,2)-1);
avail_sd = nan(size(delta_sd));
emp_mean = nan(size(ph_R,1),1);
emp_std = nan(size(ph_R,1),1);
avail_ep_num = nan(size(ph_R,1),1);

for s = 1 : size(ph_R,1)
    if (any(ph_R(s,:)) && any(ph_M(s,:)))
        ph_R_tmp(s,:) = ph_R(s,:);
        ph_M_tmp(s,:) = ph_M(s,:);
        ph_R_tmp(s,ph_R_tmp(s,:)==0) = NaN;
        ph_M_tmp(s,ph_M_tmp(s,:)==0) = NaN;
        delta_sd(s,:) = diff(ph_R_tmp(s,:) - ph_M_tmp(s,:))/interval';
        delta_sd(s,delta_sd(s,:) == 0) = NaN;
        not_zero = find(delta_sd(s,:) ~= 0);
        not_nan  = find(~isnan(delta_sd(s,:)));
        avail_ep = intersect(not_zero, not_nan);
        avail_ep_num(s,1) = length(avail_ep);
        avail_sd(s,1:avail_ep_num(s,1)) = avail_ep;

        if (~isempty(delta_sd(s,avail_sd(s,1:avail_ep_num(s,1)))))
            outliers = batch_outlier_detection(delta_sd(s,avail_sd(s,1:avail_ep_num(s,1))),median(round(interval)));
            [~,jmp_sd] = intersect(delta_sd(s,:),outliers);
        end

        %figure; plot(avail_sd(s,1:avail_ep_num(s,1)), delta_sd(s,avail_sd(s,1:avail_ep_num(s,1))),'.')

        jmp_sd = sort(jmp_sd);
        delta_sd(s,jmp_sd) = NaN;
        p = setdiff(avail_sd(s,1:avail_ep_num(s,1)), jmp_sd);
        avail_ep_num(s,1) = avail_ep_num(s,1) - length(jmp_sd);
        avail_sd(s,1:avail_ep_num(s,1)) = p;
        ph_R(s, jmp_sd+1) = 0;
        ph_M(s, jmp_sd+1) = 0;

        %mean(delta_sd(s,avail_sd(s,1:avail_ep_num(s,1))))
        %std(delta_sd(s,avail_sd(s,1:avail_ep_num(s,1))))
        emp_mean(s,1) = mean(delta_sd(s,avail_sd(s,1:avail_ep_num(s,1))));
        emp_std(s,1) = std(delta_sd(s,avail_sd(s,1:avail_ep_num(s,1))));
    end
end

emp_std(isnan(emp_std)) = [];
emp_std = median(emp_std);

for s = 1 : size(ph_R,1)
    if (any(ph_R(s,:)) && any(ph_M(s,:)))
        jmp_thres = find(abs(delta_sd(s,:)-emp_mean(s,1)) > 3*emp_std)';
%         ph_R(s, jmp_thres) = 0;
%         ph_M(s, jmp_thres) = 0;
        ph_R(s, jmp_thres+1) = 0;
        ph_M(s, jmp_thres+1) = 0;
    end
end
