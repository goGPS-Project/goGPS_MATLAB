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
Fig_Lab.plotExtractionPos('/Users/Andrea/Repositories/goGPS_MATLAB/data/project/CAC1-CAC3/out/BlockYear/CAC1_CAC3_2016215_2017202_extraction_merge.txt');
subplot(3,1,1); ylim([-10 10]); subplot(3,1,2); ylim([-10 10]); subplot(3,1,3); ylim([-10 10]);
subplot(3,1,1); ylim([-5 5]); subplot(3,1,2); ylim([-5 5]); subplot(3,1,3); ylim([-10 10]);
