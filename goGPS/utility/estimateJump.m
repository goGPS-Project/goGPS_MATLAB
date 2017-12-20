function jmp = estimateJump(series,jmp_idx,mode)
    %Description: estimate jmp in a series doing a dtrend using medians 
    %INPUT:
    % series : the series
    % jmp_idx : where the jmp occurs
    % mode : 1 joint detrend , 2: detrend separately
    diff_pre = diff(series(1:jmp_idx-1));
    diff_post = diff(series(jmp_idx:end));
    if mode == 1
       trnd = median([diff_pre ;diff_post],'omitnan');
      series = series - trnd * [1:length(series)]';
    else
        trend1 = median([diff_pre],'omitnan');
        trend2 = median([diff_post],'omitnan');
        %series(1:jmp_idx-1) = series(1:jmp_idx-1) - trend1 * [1 : jmp_idx-1]';
        %series(jmp_idx:end) = series(jmp_idx:end) - trend2 * [1 : (length(series) - jmp_idx + 1)]';
    end
    jmp = median(series(jmp_idx:end),'omitnan') - median(series(1:jmp_idx-1),'omitnan')  ;
end
    
    