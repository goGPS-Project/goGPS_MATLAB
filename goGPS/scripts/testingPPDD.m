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
        [ph{r}, wl, id_ph] = work.getPhases();
        go_id{r} = work.go_id(id_ph);
        obs_code = work.obs_code(id_ph, :);
        for f = 1 : size(obs_code, 1)
            [~, freq_offset{r}(f)] = intersect(Core_Utils.code3Char2Num(codes{r}), Core_Utils.code3Char2Num(obs_code(f,:)));
        end
        freq_offset{r} = freq_offset{r}(:) - 1;
        synt_obs{r} = work.getSyntPhases();
        
        obs_diff = (ph{r} - synt_obs{r});
        
        % First approach
        est_clock = detrend(cumsum(mean(Core_Utils.diffAndPred(obs_diff), 2, 'omitnan')));
        ph_clean1{r} = bsxfun(@minus, ph{r}, est_clock);
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
        sat_diff = zeros(size(id_sync, 1), 3);
        n_freq_sat = zeros(size(id_sync, 1), 3);
        
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

