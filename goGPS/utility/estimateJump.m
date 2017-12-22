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
        idx_bf = find(valid_idx_l(1:jmp_idx - 1));
        idx_aft = jmp_idx -1  + find(valid_idx_l(jmp_idx:end));
        n_valid = sum(valid_idx);
        if same_slope
            
        else
        end
    end
end
    
    