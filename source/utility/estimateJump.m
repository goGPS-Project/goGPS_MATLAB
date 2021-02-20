function jmp = estimateJump(series,jmp_idx,same_slope,mode)
    %Description: estimate jmp in a series doing a dtrend using medians
    %INPUT:
    % series : the series
    % jmp_idx : where the jmp occurs
    % mode : 1 joint detrend , 2: detrend separately
    if nargin < 3
        same_slope = true;
        mode = 'median';
    end
    if nargin < 4
        mode = 'median';
    end
    if strcmp(mode, 'median')
        diff_pre = diff(series(1:jmp_idx-1));
        diff_post = diff(series(jmp_idx:end));
        if same_slope
            trnd = median([diff_pre ;diff_post],'omitnan');
            series = series - trnd * [1:length(series)]';
        else
            trend1 = median([diff_pre],'omitnan');
            trend2 = median([diff_post],'omitnan');
            %series(1:jmp_idx-1) = series(1:jmp_idx-1) - trend1 * [1 : jmp_idx-1]';
            %series(jmp_idx:end) = series(jmp_idx:end) - trend2 * [1 : (length(series) - jmp_idx + 1)]';
        end
        jmp = median(series(jmp_idx:end),'omitnan') - median(series(1:jmp_idx-1),'omitnan')  ;
    elseif strcmp(mode, 'ls')
        valid_idx_l = ~isnan(series);
        valid_idx = find(valid_idx_l);
        idx_bf = valid_idx < jmp_idx;
        idx_aft = valid_idx >= jmp_idx;
        n_valid = sum(valid_idx_l);
        if same_slope
            A= zeros(n_valid,3);
            A(:,1) = find(valid_idx_l);
            A(idx_bf,2) = 1;
            A(idx_aft,3) = 1;
            x = A\series(valid_idx_l);
            jmp = x(3) -x(2);
        else
            A= zeros(n_valid,4);
            A(idx_bf, 1) = valid_idx(idx_bf);
            A(idx_aft, 2) = valid_idx(idx_aft);
            A(idx_bf, 3) = 1;
            A(idx_aft, 4) = 1;
            x = A\series(valid_idx_l);
            jmp = x(4) -x(3);
        end
    end
end
    
    