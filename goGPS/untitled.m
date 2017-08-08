[time, xyz, enu, std_enu, std_3d] = steCrdImporter ('/Users/Andrea/Repositories/goGPS_MATLAB/data/project/default_DD_batch/station/CRD/CAC2-(CAC3).CRD');
time.addIntSeconds(86400);
%%
ids = time.getMatlabTime()>datenum('28 Nov 2016') & time.getMatlabTime()<=datenum('28 Dec 2016');
Fig_Lab.plotExtractionPos('../data/project/default_DD_batch/outFloat_FB/CAC2_CAC3_2016333_2016363_extraction.txt')
Fig_Lab.plotENU(time.getId(ids), enu(ids,:), 14, -1)
Fig_Lab.plotExtractionPos('/Users/Andrea/Repositories/goGPS_MATLAB/data/project/default_DD_batch/outBlock/CAC2_CAC3_2016333_2016363_extraction.txt', 1);

%%

Fig_Lab.plotExtractionPos('../data/project/PET1-CAC1/outLong/CAC1_PET1_2016055_2016161_extraction.txt')
Fig_Lab.plotExtractionPos('../data/project/PET1-CAC1/outGatTest/CAC1_PET1_2016070_2016073_extraction.txt',1)


%%
figure;
Fig_Lab.plotENU(time.getId(ids), enu(ids,:), 14, 0);
Fig_Lab.plotExtractionPos('/Users/Andrea/Repositories/goGPS_MATLAB/data/project/default_DD_batch/outBlock/CAC2_CAC3_2016333_2016363_extraction.txt', 1);

%%
profile off
profile on
[pr1_RM, ph1_RM, pr2_RM, ph2_RM, dop1_RM, dop2_RM, snr1_RM, snr2_RM, ...
                    zero_time, time_GPS_diff, time_RM_diff, week_RM, date_RM, pos_RM, interval, antoff_RM, antmod_RM, codeC1_RM, marker_RM] = ...
                    load_RINEX_obs(filename_obs, state.getConstellationCollector(), processing_interval);
profile off
profile viewer

%%
[pr1_RM, ph1_RM, pr2_RM, ph2_RM, dop1_RM, dop2_RM, snr1_RM, snr2_RM, ...
                    zero_time, time_GPS_diff, time_RM_diff, week_RM, date_RM, pos_RM, interval, antoff_RM, antmod_RM, codeC1_RM, marker_RM] = ...
                    load_RINEX_obs(filename_obs, state.getConstellationCollector(), processing_interval);

tic;
[pr1_RM_o, ph1_RM_o, pr2_RM_o, ph2_RM_o, dop1_RM_o, dop2_RM_o, snr1_RM_o, snr2_RM_o, ...
                    zero_time, time_GPS_diff, time_RM_diff, week_RM, date_RM, pos_RM, interval, antoff_RM, antmod_RM, codeC1_RM, marker_RM] = ...
                    load_RINEX_obs_old(filename_obs, state.getConstellationCollector(), processing_interval);
toc;

%%
    tic
    fid = fopen(filename_obs{1}, 'r');
    txt = fread(fid,'*char')';
    fclose(fid);
    txt(txt == 13) = [];
    eol = [0 find(txt == 10)];
    line = cell(numel(eol),1);
    for l = 1 : numel(eol)-1
        line{l} = txt(eol(l)+1:eol(l+1)-1);
    end
    toc
    tic
    for l = 1 : numel(eol)-1
        tmp = txt(eol(l)+1:eol(l+1)-1);
    end
    toc

    tic
    for l = 1 : numel(eol)-1
        tmp = line{l};
    end
    toc

    tic
    fid = fopen(filename_obs{1}, 'r');
    line2 = textscan(fid, '%s', 'Delimiter', '\n', 'whitespace', ''); line2 = line2{1};
    fclose(fid);
    toc

%%
Fig_Lab.plotExtractionPos('/Users/Andrea/Repositories/goGPS_MATLAB/data/project/CAC1-CAC3/out/BlockYear/CAC1_CAC3_2016215_2017202_extraction.txt');
subplot(3,1,1); ylim([-10 10]); subplot(3,1,2); ylim([-10 10]); subplot(3,1,3); ylim([-10 10]);
subplot(3,1,1); ylim([-4 4]); subplot(3,1,2); ylim([-7 7]); subplot(3,1,3); ylim([-9 9]);
Fig_Lab.plotExtractionPos('/Users/Andrea/Repositories/goGPS_MATLAB/data/project/CAC1-CAC3/out/BlockYear/CAC1_CAC3_2016215_2017202_extraction_old.txt', -1);

%%

time_span = 3600;
rate = 15;
time_step = [(time_span / rate) : (time_span / rate) : sat_track(end,1)  sat_track(end,1)]' ;
block_id_lim = zeros(numel(time_step),2);
s = 1;
i = 1;
last_i = 1;
while i < size(sat_track, 1)
    while (sat_track(i,1) > time_step(s))
        s = s + 1;
    end
    while  i < size(sat_track, 1) && (sat_track(i,1) < time_step(s)+1) 
        i = i + 1;
    end
    block_id_lim(s, 1) = last_i;
    block_id_lim(s, 2) = i - 1;
    last_i = i;
    i = i + 1;
end

% Building the new Design matrix
valid_block_id_lim = block_id_lim(block_id_lim(:,1) > 0, :);
num_est = size(valid_block_id_lim,1);

A2 = sparse(size(A, 1), size(A, 2) + 3 * (num_est - 1));
for i = 1 : num_est
    id = valid_block_id_lim(i,1) : valid_block_id_lim(i,2);
    A2(id, (3 * (i - 1) + 1 : 3 * i)) = A(id,1:3);
end
i = i + 1;
A2(:,(3 * (i - 1) + 1) : end) = A(:, 4 : end);

[x_hr, Cxx, sigma02_hat, v_hat] = fast_least_squares_solver(y0, b, A2, Q);
x_hr = reshape(x_hr(1 : 3 * num_est), 3, num_est)';
x_m = 0.5 * median(x_hr) + 0.5 * mean(x_hr);

figure;
subplot(3,1,1); plot(time_step(block_id_lim(:,1) > 0) * rate, 1*(x_hr(:,1) - x(1)),'.r', 'MarkerSize', 20); hold on;
subplot(3,1,2); plot(time_step(block_id_lim(:,1) > 0) * rate, 1*(x_hr(:,2) - x(2)),'.g', 'MarkerSize', 20); hold on;
subplot(3,1,3); plot(time_step(block_id_lim(:,1) > 0) * rate, 1*(x_hr(:,3) - x(3)),'.b', 'MarkerSize', 20); hold on;

% Post fix:

%[deltaX, estim_amb, sigma_amb, sigma_pos] = lambdafix(x(1:3), x(4:end), cov_X, cov_N, cov_XN);
estim_amb = estim_amb_bk;
num_amb = numel(estim_amb);
num_obs = size(A2,1);
A_fix = [A2; sparse(num_amb, size(A2, 2))];
y0_fix = [y0; estim_amb];
b_fix = [b; zeros(num_amb, 1)];
Q_fix = [[Q sparse(size(Q, 1), num_amb)]; sparse(num_amb, size(Q, 1) + num_amb)]; 
for i = 1 : numel(estim_amb)
    A_fix(num_obs + i, 3 * num_est + i) = 1;
    Q_fix(num_obs + i,num_obs + i) = 0.001;
end

[x_hr, Cxx, sigma02_hat, v_hat] = fast_least_squares_solver(y0_fix, b_fix, A_fix, Q_fix);
x_hr = reshape(x_hr(1 : 3 * num_est), 3, num_est)';
x_m = 0.5 * median(x_hr) + 0.5 * mean(x_hr);

subplot(3,1,1); plot(time_step(block_id_lim(:,1) > 0) * rate, 1*(x_hr(:,1) - x(1)),'.k', 'MarkerSize', 20); hold on;
subplot(3,1,2); plot(time_step(block_id_lim(:,1) > 0) * rate, 1*(x_hr(:,2) - x(2)),'.k', 'MarkerSize', 20); hold on;
subplot(3,1,3); plot(time_step(block_id_lim(:,1) > 0) * rate, 1*(x_hr(:,3) - x(3)),'.k', 'MarkerSize', 20); hold on;


%%
%switch from SD to DD
                        D = zeros(amb_num-1,amb_num);
                        D(:,1) = 1;
                        D(:,2:end) = -eye(amb_num-1);
                        G = zeros(3*num_est+amb_num-1,3*num_est+amb_num);
                        G(1:3*num_est,1:3*num_est) = eye(num_est*3);
                        G((3*num_est+1):end,(3*num_est+1):end) = D;
                        x = G*x;
                        Cxx = G*Cxx*G';

                        cov_X  = Cxx(1:3*num_est,1:3*num_est);     %position covariance block
                        cov_N  = Cxx((3*num_est+1):end,(3*num_est+1):end); %ambiguity covariance block
                        cov_XN = Cxx(1:3*num_est,(3*num_est+1):end);   %position-ambiguity covariance block

                        fixed_amb = 0;
                        try
                            [U] = chol(cov_N);
                            cov_N = U'*U;
                        catch ex
                            logger.addWarning(sprintf('Covariance matrix unstable %s', ex.message));
                        end

                        %integer phase ambiguity solving by LAMBDA
                        [deltaX, estim_amb, sigma_amb, sigma_pos] = lambdafix(x(1:3*num_est), x((3*num_est+1):end), cov_X, cov_N, cov_XN);
                        %%
%%
x_hr = zeros(ceil(86400 / time_span)-1,3);
for i = 1 : size(block_id_lim, 1)
    id = block_id_lim(i, 1) : block_id_lim(i, 2);
    if id(1) > 0
        id2 = full(sum(A(id,:)~=0, 1)) > 0;
        tmp = (A(id,id2)'/Q(id,id)*A(id,id2)) \ (A(id,id2)'/Q(id,id)) * (y0(id)-b(id));
        x_hr(i,:) = tmp(1:3);
    else
        x_hr(i,:) = NaN;
    end 
end
subplot(3,1,1); plot(time_step * rate, 1*(x_hr(:,1) - x(1)),'.', 'MarkerSize', 15, 'Color', [0.5 0.5 0.5]); hold on;
subplot(3,1,2); plot(time_step * rate, 1*(x_hr(:,2) - x(2)),'.', 'MarkerSize', 15, 'Color', [0.5 0.5 0.5]); hold on;
subplot(3,1,3); plot(time_step * rate, 1*(x_hr(:,3) - x(3)),'.', 'MarkerSize', 15, 'Color', [0.5 0.5 0.5]); hold on;
