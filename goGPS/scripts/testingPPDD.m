% %%
% 
% A_idx_s = ls.A_idx;
% A_s = ls.A_ep;
% sat_s = ls.sat;
% rec_id_s = ls.receiver_id;
% epoch_s = ls.epoch;
% y_s = ls.y;
% %%
% A_s = A_rec;
% A_idx_s = Aidx_rec;
% sat_s = sat_rec;
% rec_id_s = ones(size(y_rec))*i;
% epoch_s = ep_rec;
% y_s = y_rec;
% %%
% 
% A_s = A;
% A_idx_s = Aidx;
% sat_s = sat;
% rec_id_s = r;
% epoch_s = ep;
% y_s = y;
% 
% %% Aidx to A
% tic;
% A_full = sparse(size(A_idx_s,1), max(A_idx_s(:)));
% A_serial_id = serialize(repmat((1 : size(A_idx_s,1))', 1, size(A_idx_s,2)) + (A_idx_s - 1) * size(A_idx_s, 1));
% type = nan(max(A_idx_s(:)), 1);
% for t = 1 : size(A_idx_s, 2)
%     type(A_idx_s(:, t)) = t;
% end
% A_full(A_serial_id) = A_s(:);
% toc;
% %%
% tic
% max_ep = 1000;
% max_sat = max(sat_s(:));
% sat_err = nan(max_ep, max(sat_s), max(rec_id_s));
% sat_obs = nan(max_ep, max(sat_s), max(rec_id_s));
% id_ok = epoch_s < max_ep;
% col_ok = sum(abs(A_full(id_ok, :)) ~= 0);
% type_tmp = p_class_rec(type(col_ok > 0));
% A_tmp = A_full(id_ok, col_ok > 0);
% A_tmp(:, type_tmp == 6) = []; % remove clock estimation
% type_tmp(type_tmp == 6) = []; % remove clock estimation
% 
% toc
% %%
% x1 = A_tmp \ y_s(id_ok);
% 
% y1 = (A_tmp * x1); 
% 
% sat_obs(epoch_s(id_ok) + (sat_s(id_ok) - 1) * max_ep + (rec_id_s(id_ok) - 1) * (max_sat * max_ep)) = y_s(id_ok);
% %figure; plot(sat_obs(:,:,1))
% %figure; plot(sat_obs(:,:,2))
% 
% 
% sat_err(epoch_s(id_ok) + (sat_s(id_ok) - 1) * max_ep + (rec_id_s(id_ok) - 1) * (max_sat * max_ep)) = y_s(id_ok) - y1;
% figure; plot(sat_err(:,:,1))
% figure; plot(sat_err(:,:,2))
% figure; plot(sum(sat_err(:,:,:),3));
% toc
% 
% %%
%             % reference
%             [synt_obs_ref, xs_loc_ref] = rec(1).getSyntTwin(obs_set(2));
%             xs_loc_ref = zero2nan(xs_loc_ref);
%             diff_obs_ref = nan2zero(zero2nan(obs_set(1).obs) - zero2nan(synt_obs_ref));
% 
%             % target
%             [synt_obs_trg, xs_loc_trg] = rec(2).getSyntTwin(obs_set(2));
%             xs_loc_trg = zero2nan(xs_loc_trg);
%             diff_obs_trg = nan2zero(zero2nan(obs_set(2).obs) - zero2nan(synt_obs_trg));
%             
%             xs_loc = xs_loc_ref - xs_loc_trg;
%             
%             % Removing common clock error (between trg and ref receivers)
%             dd = zero2nan(diff_obs_ref) - zero2nan(diff_obs_trg);
%             sensor = Core_Utils.diffAndPred(dd);
%             sensor = median(sensor, 2, 'omitnan'); % using median to diregard cycle-slips
%             % rough estimation of clock, the median is not a good estimator
%             % but for now it could stay like this
%             clock_diff = detrend(cumsum(sensor));
%             % Estimate clock_diff (and remove it, like if I was doing DD)
%             dd = dd - clock_diff;   
%         
%         
% %%
% 
% % Aidx to A
% tic;
% A_full = sparse(size(Aidx,1), max(Aidx(:)));
% A_serial_id = serialize(repmat((1 : size(Aidx,1))', 1, size(Aidx,2)) + (Aidx - 1) * size(Aidx, 1));
% A_full(A_serial_id) = A(:);
% toc;
% %%
% tic;
% max_ep = 300;
% max_sat = max(sat(:));
% sat_err = nan(max_ep, max(sat), max(r));
% id_ok = ep < max_ep;
% A_tmp = A_full(id_ok, sum(abs(A_full(id_ok, :))) > 0);
% x = A_tmp \ y(id_ok);
% 
% y1 = (A_tmp*x);
% 
% sat_err(ep(id_ok) + (sat(id_ok) - 1) * max_ep + (r(id_ok)-1) * (max_sat * max_ep)) = y(id_ok) - y1;
% figure; plot(sat_err(:,:,1))
% figure; plot(sat_err(:,:,2))
% toc
% 

%% TESTING OUT DET
% run this only after a good pre-processing or PPP
tic
n_rec = numel(rec);
sys_list = rec(1).work.getAvailableSS();
work_list = [rec.work];
[s_time, id_sync] = Receiver_Commons.getSyncTimeExpanded(work_list);
for sys_c = sys_list
    % for each constellation
    freq = {};
    
    n_sat = 0;
    codes = {};
    for r = 1 : n_rec
        work = rec(r).work;
        work.remShortArc(max(2, core.state.getMinArc));
        % supposing only one observation set per band
        % at this point I aspect all the "worse" frequencies to be purged from the work object
        codes{r} = work.getAvailableObsCode('L??', sys_c);
        n_sat = max([n_sat; work.go_id]);
    end
    n_freq = numel(unique(Core_Utils.code3Char2Num(cell2mat(codes'))));
    
    obs_diff_mr = nan(size(id_sync, 1), n_sat * n_freq, n_rec);
    
    for r = 1 : n_rec
        work = rec(r).work;
        freq_offset{r} = [];
        
        % Prepare clock free data
        [ph{r}, wl{r}, id_ph{r}] = work.getPhases();
        go_id{r} = work.go_id(id_ph{r});
        obs_code = work.obs_code(id_ph{r}, :);
        for f = 1 : size(obs_code, 1)
            [~, freq_offset{r}(f)] = intersect(Core_Utils.code3Char2Num(codes{r}), Core_Utils.code3Char2Num(obs_code(f,:)));
        end
        freq_offset{r} = freq_offset{r}(:) - 1;
        synt_obs{r} = work.getSyntPhases();
        
        obs_diff = (ph{r} - synt_obs{r});
        
        % First approach
        %if any(work.dt)
        %    est_clock{r} = work.dt .* Global_Configuration.V_LIGHT;
        %else
        est_clock{r} = detrend(cumsum(mean(Core_Utils.diffAndPred(obs_diff), 2, 'omitnan')));
        %end
        ph_clean1{r} = bsxfun(@minus, ph{r}, est_clock{r});
        obs_diff1 = (ph_clean1{r} - synt_obs{r});
        
        obs_diff_mr(~isnan(id_sync(:, r)), go_id{r} + freq_offset{r} * n_sat, r) = obs_diff1(noNaN(id_sync(:,r)), :);
    end
    
    % Using multiple receivers I remove the common part among all the observations of the same satellite
    cmn = nan(size(id_sync, 1), n_sat);
    inan = all(isnan(obs_diff_mr), 3);
    % Removing satellite common part
    for s = 1 : n_sat
        % Computing the mean among frequencies
        % this can be improved using the error of the frequency
        sat_diff = zeros(size(id_sync, 1), n_rec);
        n_freq_sat = zeros(size(id_sync, 1), n_rec);
        
        % Use multi frequencies as more observations for a single receiver
        % Compute the mean among all the observed frequencies
        %figure;
        for f = unique(freq_offset{r}')
            id_ok = squeeze(~isnan(obs_diff_mr(:,s + f * n_sat, :)));
            n_freq_sat = n_freq_sat + id_ok;
            for r = 1 : n_rec
                if any(id_ok(:,r))
                    tmp = Core_Utils.diffAndPred(squeeze(obs_diff_mr(:, s + f * n_sat, r)));
                    sat_diff(id_ok(:,r), r) = sat_diff(id_ok(:,r), r) + tmp(id_ok(:,r));
                    %hold on; plot(tmp); dockAllFigures
                end
            end
        end
        sat_diff = sat_diff ./ n_freq_sat;
        cmn(:, s) = detrend(cumsum(nan2zero(median(sat_diff, 2, 'omitnan'))));
    end
    cmn = repmat(cmn, 1, n_freq);
    cmn(inan) = nan;
    cmn = bsxfun(@minus, cmn, detrend(cumsum(median(Core_Utils.diffAndPred(cmn), 2, 'omitnan'))));
    
    %     for r = 1 : n_rec
    %         work = rec(r).work;
    %         work.keepEpochs(id_sync(:,r));
    %         ph_new{r} = ph_clean1{r}(noNaN(id_sync(:, r)), :) - cmn(~isnan(id_sync(:, r)), go_id{r} + freq_offset{r} * n_sat);
    %         work.setPhases(ph_new{r}, wl{r}, id_ph{r});
    %     end
end
toc

for r = 1 : n_rec
    fprintf('Receiver %d\n', r);
    data_diff0 = ph{r}(noNaN(id_sync(:, r)), :) - synt_obs{r}(noNaN(id_sync(:, r)), :);
    data_diff1 = ph_clean1{r}(noNaN(id_sync(:, r)), :) - synt_obs{r}(noNaN(id_sync(:, r)), :);
    data_diff2 = ph_clean1{r}(noNaN(id_sync(:, r)), :) - cmn(~isnan(id_sync(:, r)), go_id{r} + freq_offset{r} * n_sat) - synt_obs{r}(noNaN(id_sync(:, r)), :);
    
    fprintf(' - std %.2f mm \n', [mean(std(diff(data_diff0, 2), 'omitnan')), ...
        mean(std(diff(data_diff1, 2), 'omitnan')), ...
        mean(std(diff(data_diff2, 2), 'omitnan'))] * 1e3);
    
    ax = [];
    fh = figure; fh.NumberTitle = 'off'; fh.Name = sprintf('%d) R%d stage0', fh.Number, r); plot(data_diff0, '.-'); ax(1) = gca();
    fh = figure; fh.NumberTitle = 'off'; fh.Name = sprintf('%d) R%d stage1', fh.Number, r); plot(data_diff1, '.-'); ax(2) = gca();
    fh = figure; fh.NumberTitle = 'off'; fh.Name = sprintf('%d) R%d stage2', fh.Number, r); plot(data_diff2, '.-'); ax(3) = gca();
    linkaxes(ax);
    fh = figure; fh.NumberTitle = 'off'; fh.Name = sprintf('%d) R%d diff stage0', fh.Number, r); plot(diff(data_diff0), '.-'); ylim([-1 1] * 2);
    fh = figure; fh.NumberTitle = 'off'; fh.Name = sprintf('%d) R%d diff stage1', fh.Number, r); plot(diff(data_diff1), '.-'); ylim([-1 1] * 0.1);
    fh = figure; fh.NumberTitle = 'off'; fh.Name = sprintf('%d) R%d diff stage2', fh.Number, r); plot(diff(data_diff2), '.-'); ylim([-1 1] * 0.1);
    dockAllFigures
end

