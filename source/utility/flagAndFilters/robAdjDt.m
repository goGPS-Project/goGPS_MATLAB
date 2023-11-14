function [dts] = robAdjDt(dph)
    % wrapper of fallback in case no mex
    % it is not the same code, instead of huber it uses the simple median
    
    % Fast simple approach: dt = median(sensor, 2, 'omitnan');
    
    [n_ep, n_stream] = size(dph);
    thrs = 0.02;
    
    dts = zeros(n_ep, 1);

    % for all the epoch
    parfor k = 1:n_ep
        % get the phases of the row
        dph_tmp = dph(k, :);
        id_ok = isfinite(dph_tmp);
        dph_tmp2 = dph_tmp(id_ok);
        if ~isempty(dph_tmp2) % if we have phases
            w = ones(size(dph_tmp2));
            j = 0;
            dt = 1e9;
            dt_prev = -1e9;
            while (j < 20 && abs(dt - dt_prev) > 0.005) % limit the reweight to 30 or less than 0.005 improvement
                dt_prev = dt;
                tmp = dph_tmp2 .* w;
                dt = sum(tmp) / sum(w); % weighted mean
                
                ares_n = abs(dph_tmp2 - dt) / thrs; % absolute residuals
                w = ones(size(dph_tmp2));
                idx_rw = find(ares_n > 1); % residual to be reweighted
                if ~isempty(idx_rw)
                    w(idx_rw) = 1 ./ ares_n(idx_rw).^2; % compute the weight
                end
                j = j + 1;  
            end
            dts(k) = dt; % put in the results
        end
    end
end

end
