%% CODE residual weights

%% % remS observations with an SNR smaller than snr_thr
%
% SYNTAX
%   this.remUnderSnrThr(abs_snr_thr, scaled_snr_thr)

nargout = 2;
profile off;
profile on;
flag_debug = true;

l_max = 11;
m_max = 11;
std_thr_offset = 0.5; % Every data with variance > 0.75m above the weight mask is considered an outlier!

pr_grid_step = 1;
pr_grid = (1 : pr_grid_step : 70);
cc = Core.getConstellationCollector();
fh_id = 0;

tic;
% for each satellite system

z_pol_err = []; 
clear z_pol_err;
for sys_c = this.getAvailableSS
    cur_ss = cc.(lower(cc.getSysName(sys_c)));
    ss_obs_code = this.getAvailableObsCode('C', sys_c);
    
    % Read all pr and put them in a 3-D matrix
    [pr, id_pr] = this.getObs('C??', sys_c);
    go_id_list = unique(this.go_id(id_pr));    
    prs = this.getSyntObs(go_id_list);
    
    obs_code_num = Core_Utils.code2Char2Num(this.obs_code(id_pr,2:3));
    obs_code_list = unique(obs_code_num);
    
    n_obs = size(pr, 2);
    n_sat = numel(go_id_list);
    n_cod = numel(obs_code_list);
    
    noise_pr = zeros(n_obs, n_sat, n_cod);
    id2noise = zeros(size(pr, 1), 2);
    for  i = 1 : size(pr, 1)
        o = find(obs_code_list == obs_code_num(i));
        s = find(go_id_list == this.go_id(id_pr(i)));
        id2noise(i, 1) = s;
        id2noise(i, 2) = o;
        
        % remove synthetic
        tmp = (pr(i, :) - prs(s, :))';
        % remove DCB
        tmp = tmp - mean(tmp, 1, 'omitnan');
        noise_pr(:, s, o) = tmp;                
    end
    
    % remove the common part (tropo + iono)
    for s = 1 : n_sat
        tmp = bsxfun(@minus, zero2nan(squeeze(noise_pr(:, s, :))), mean(squeeze(noise_pr(:, s, :)), 2, 'omitnan'));         
        % remove 10 minutes splines
        tmp = Receiver_Commons.smoothSatData([],[],zero2nan(mean(noise_pr(:, s, :), 3, 'omitnan')), [], 'spline', 600 / this.getRate, 0);
        noise_pr(:, s, :) = bsxfun(@minus, noise_pr(:, s, :), tmp);
    end    
    noise_pr = zero2nan(movstd(noise_pr, 900 / this.getRate, 1, 'omitnan'));

    az = this.sat.az(:, go_id_list);
    el = this.sat.el(:, go_id_list);
    % Array of data with variance too high
    id_ko = false(size(noise_pr));
    for o = 1 : numel(obs_code_list)        
        data = noise_pr(:, :, o);
        pr_ok = ~isnan(data); % Consider that stdmov have a border effect, ignore the first and last 5 epochs
        pr_int = flagShrink(pr_ok, 5); % Consider that stdmov have a border effect, ignore the first and last 5 epochs
        % if flag_debug
        %    fh = figure; clf; hold on; fh = Core_Utils.polarZerMapQuad(l_max, m_max, az(pr_int)/180*pi, el(pr_int)/180*pi, data(pr_int)); colormap([[1 1 1]; flipud(Cmap.get('plasma')); [0.6 0.6 0.6]]);
        %    lim = [0 1]; subplot(2,2,1); caxis(lim); subplot(2,2,2); caxis(lim); subplot(2,2,3); ylim(lim); subplot(2,2,4); ylim(lim);
        %    fh.Name = sprintf('%d) PRE %c C%s', fh.Number, sys_c, Core_Utils.num2Code2Char(obs_code_list(o)));
        %    fh.NumberTitle = 'off';
        % end
        
        [z_par, l, m] = Core_Utils.zAnalisysAll(l_max, m_max, az(pr_int)/180*pi, el(pr_int)/180*pi, data(pr_int));
        % all the data with a variance bigger than thr_offset w.r.t. std_map is considered outlier
        id_ko(find(pr_ok) + (o-1) * n_obs * n_sat) = (data(pr_ok) - (std_thr_offset/3*2) - Core_Utils.zSinthesys(l, m, az(pr_ok)/180*pi, el(pr_ok)/180*pi, z_par)) > 0;
        
        if any(id_ko(:))
            pr_int = pr_int & ~id_ko;
            %  reestimate zernike without outliers
            z_par = Core_Utils.zAnalisys(l, m, az(pr_ok)/180*pi, el(pr_ok)/180*pi, data(pr_ok));
            % all the data with a variance bigger than thr_offset w.r.t. std_map is considered outlier
            if (nargout >= 2)
                id_ko(find(pr_ok) + (o-1) * n_obs * n_sat) = (data(pr_ok) - (std_thr_offset/3*2) - Core_Utils.zSinthesys(l, m, az(pr_ok)/180*pi, el(pr_ok)/180*pi, z_par)) > 0;
            end
        end
        
        if flag_debug && any(pr_ok(:))
            pr_ok = pr_ok & ~id_ko(:, :, o);                
            fh = figure; clf; hold on; fh = Core_Utils.polarZerMapQuad(l_max, m_max, az(pr_ok)/180*pi, el(pr_ok)/180*pi, data(pr_ok));  colormap([[1 1 1]; flipud(Cmap.get('plasma')); [0.6 0.6 0.6]]);
            lim = [0 1]; subplot(2,2,1); caxis(lim); subplot(2,2,2); caxis(lim); subplot(2,2,3); ylim(lim); subplot(2,2,4); ylim(lim);
            fh.Name = sprintf('%d) POST %c C%s', fh.Number, sys_c, Core_Utils.num2Code2Char(obs_code_list(o)));
            fh.NumberTitle = 'off';
        end
        z_pol_err.(sprintf('%cC%s', sys_c, Core_Utils.num2Code2Char(obs_code_list(o)))) = struct('z_pol', z_par, 'l', l, 'm', m);
    end
end
% Export id_ok to the same format of all pr
if (nargout >= 2)
    id_ko = id_ko(:, id2noise(:,1) + n_sat * (id2noise(:,2) - 1));
end
toc
%profile off;
%profile viewer;






%% Display the Zernike function Z(n=5,m=1)
x = -1:0.01:1;
[X,Y] = meshgrid(x,x);
[theta,r] = cart2pol(X,Y);
idx = r<=1;
z = nan(size(X));
z(idx) = zernfun(4,2,r(idx),theta(idx));
figure(100)
mesh(x,x,z), shading interp
axis square, colorbar
title('Zernike function Z_5^1(r,\theta)')

%% Example 2:

% Display the first 10 Zernike functions
x = -1:0.01:1;
[X,Y] = meshgrid(x,x);
[theta,r] = cart2pol(X,Y);
idx = r<=1;
z = nan(size(X));
n = [0  1  1  2  2  2  3  3  3  3];
m = [0 -1  1 -2  0  2 -3 -1  1  3];
Nplot = [4 10 12 16 18 20 22 24 26 28];
y = zernfun(n,m,r(idx),theta(idx));
figure('Units','normalized')
for k = 1:10
    z(idx) = y(:,k);
    subplot(4,7,Nplot(k))
    pcolor(x,x,z), shading interp
    set(gca,'XTick',[],'YTick',[])
    axis square
    title(['Z_{' num2str(n(k)) '}^{' num2str(m(k)) '}'])
end

%% Analysis

l_max = 21;
m_max = 21;

this = rec(1).work;

sys_c_list = unique(this.system);

cc = Core.getState.getConstellationCollector;

for sys_c = sys_c_list
    for b = 1 : 9 % try all the bands
        [snr_freq, snr_id_freq] = this.getSNR(sys_c, num2str(b));
        snr_id_freq = find(snr_id_freq);
        
        if any(snr_id_freq) && any(snr_freq(:))
            % Get all the trackings for this SNR
            obs_code =  unique(this.obs_code(snr_id_freq, 3));
            obs_code = [repmat(this.obs_code(snr_id_freq(1), 1:2), size(obs_code, 1), 1) obs_code];
            
            for trk = obs_code(:,3)'
                id_ok = this.obs_code(snr_id_freq, 3) == trk;
                snr_id = snr_id_freq(id_ok);
                snr = snr_freq(:, id_ok);
                
                id_ok = (~isnan(snr));
                res = this.sat.res(:,this.go_id(snr_id));
                az = this.sat.az(:,this.go_id(snr_id));
                el = this.sat.el(:,this.go_id(snr_id));
            end
        end
    end
end
%id_ok(4:end,:) = 0;
data = snr(id_ok);
az = az(id_ok)/180*pi;
el = el(id_ok)/180*pi;

% Filter
[z, z_par, l, m] = Core_Utils.zFilter(l_max, m_max, az, el, data, 1);

% Synthesis

f = figure(100);
%scatter3(x, y, z, 5, z, 'filled'), shading interp
polarScatter(az, pi/2 - el, 45, z, 'filled');
colormap(jet);  cax = caxis();
caxis([min(cax(1), 4), max(cax(2), 60)]);
setColorMap('jet', [10 55], 0.9);
colorbar();
h = title(sprintf('S%d%c - receiver %s - %s', b, trk, this.parent.marker_name, cc.getSysExtName(sys_c)),'interpreter', 'none');
h.FontWeight = 'bold';
Core_UI.beautifyFig(f);
Core_UI.addBeautifyMenu(f);
f.Visible = 'on';
title(sprintf('Zernike function up to degree %d)', l_max))

% Synthesis
tic
x = -1 : 0.005 : 1;
y = x;
[X,Y] = meshgrid(x,x);
[theta, r_synt] = cart2pol(X,Y);
idx = r_synt <= 1;
z = nan(size(X));
z(idx) = zernfun(l, m, r_synt(idx), theta(idx)) * z_par;
toc
figure(101)
mesh(x, y, z, 'LineWidth', 2); shading interp

axis square, colorbar
caxis([4 60]);
colormap(jet);
zlim([4 60]);
title(sprintf('Zernike function up to degree %d)', l_max))

% Whaaaaaat
Core_Utils.showZerniche(l, m, z_par); caxis([4 60]);

%%
[z, z_par, l, m] = Core_Utils.zFilter(l_max, m_max, az, el, data, 1);

%[az1, el1] = meshgrid(linspace(0, 2*pi, 21), linspace(0, pi/2, 19));
[az1, el1] = meshgrid([pi/2 -pi/2], linspace(0, pi/2, 11));

z1 = Core_Utils.zAnalisysAll(l_max, m_max, az1, el1, z_par);

%figure; polarScatter(az, pi/2 - el, 45, z, 'filled');
figure; polarScatter(az1(:), pi/2 - el1(:), 160, z1(:), 'filled'); colormap(jet); colorbar; caxis([4 60]);
fh = gcf; Core_UI.addBeautifyMenu(fh); Core_UI.beautifyFig(fh, 'dark');

%% SNR Zerniche test

for a = all_attrib
    
    [z_par{id}, l, m] = Core_Utils.zAnalisysAll(l_max, m_max, az(snr_ok)/180*pi, el(snr_ok)/180*pi, snr(snr_ok), 1);
    if id > 1
        snr_ref = Core_Utils.zSinthesysAll(l_max, m_max, az(snr_ok)/180*pi, el(snr_ok)/180*pi, z_par{1});
        snr_cur = Core_Utils.zSinthesysAll(l_max, m_max, az(snr_ok)/180*pi, el(snr_ok)/180*pi, z_par{id});
        scale = snr_ref./snr_cur;
        fh = Core_Utils.polarZerMapQuad(21, 21, az(snr_ok)/180*pi, el(snr_ok)/180*pi, snr_ref-snr_cur); %colormap(flipud(Cmap.get('plasma')));
        fh = Core_Utils.polarZerMapQuad(21, 21, az(snr_ok)/180*pi, el(snr_ok)/180*pi, snr(snr_ok)); %colormap(flipud(Cmap.get('plasma')));
        snr(snr_ok) = snr(snr_ok) .* scale;
    end
    fh = Core_Utils.polarZerMapQuad(21, 21, az(snr_ok)/180*pi, el(snr_ok)/180*pi, snr(snr_ok)); %colormap(flipud(Cmap.get('plasma')));
    
    
    [snr, id_snr] = this.getObs(['S' f a], sys_c);
    snr = snr';
    az = this.sat.az(:, this.go_id(id_snr));
    el = this.sat.el(:, this.go_id(id_snr));
    
    %snr_f = Receiver_Commons.smoothSatData([],[],zero2nan(snr), [], 'spline', 60 / this.getRate, 10); % smoothing SNR => to be improved
    
    id_ok = ~isnan(snr);
    %snr_f = snr;
    %snr_f(id_ok) = Core_Utils.zFilter(11, 11, az(id_ok)/180*pi, el(id_ok)/180*pi, snr(id_ok));
    
    
    %snr_f = Receiver_Commons.smoothSatData([],[],zero2nan(snr), [], 'spline', 150 / this.getRate, 10); % smoothing SNR => to be improved
    %fh = Core_Utils.polarZerMapDual(11, 11, az(id_ok)/180*pi, el(id_ok)/180*pi, abs(snr_f(id_ok) - snr(id_ok)));
    fh = Core_Utils.polarZerMapQuad(11, 11, az(id_ok)/180*pi, el(id_ok)/180*pi, snr(id_ok)); %colormap(flipud(Cmap.get('plasma')));
    
    [pr, id_pr] = this.getObs(['C' f a], sys_c);
    prs = this.getSyntObs(this.go_id(id_pr));
    pr = pr' - prs';
    pr = abs(pr - median(pr(:), 'omitnan'));
    
    %pr = movstd(Core_Utils.diffAndPred(pr,2) - median(Core_Utils.diffAndPred(pr,2),2,'omitnan'), 11, 'omitnan');
    id_ok = ~isnan(pr);
    fh = Core_Utils.polarZerMapQuad(11, 11, az(id_ok)/180*pi, el(id_ok)/180*pi, pr(id_ok)); %colormap(flipud(Cmap.get('plasma')));
    
end
