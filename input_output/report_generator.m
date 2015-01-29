function report_generator(report)

% SYNTAX:
%   report_generator(report);
%
% INPUT:
%   report      = goGPS report struct 


%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Code by Stefano Caldera
%----------------------------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%----------------------------------------------------------------------------------------------


if is_subfield(report,'opt.write') && report.opt.write == 1
    
    fout_report=fopen([report.opt.outfolder,'_report.txt'],'wt');
    fprintf(fout_report,'----------------------------------------------------------------------------------------------\n');
    fprintf(fout_report,'                           goGPS v0.4.3\n\n');
    fprintf(fout_report,' Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini\n');
    fprintf(fout_report,'----------------------------------------------------------------------------------------------\n\n');
    
    fprintf(fout_report,'INPUT OPTIONS\n');
    fprintf(fout_report,'*************\n');
    
    if report.opt.is_batch==1
        option_i = 'YES';
    else
        option_i = 'NO';
    end
    
    fprintf(fout_report,'BATCH processing                      : %s\n', option_i);
    switch report.opt.mode
        case goGNSS.MODE_RT_NAV
            option_i='Real Time Navigation (Kalman Filter on Code and Phase Double Differences (with/without a constraint)';
        case goGNSS.MODE_RT_R_MON
            option_i='% Real Time Rover Monitor';
        case goGNSS.MODE_RT_M_MON
            option_i='Real Time Master Monitor';
        case goGNSS.MODE_RT_RM_MON
            option_i='Real Time Master + Rover Monitor';
        case goGNSS.MODE_PP_LS_C_SA
            option_i='Post Proc Least Squares on Code Stand Alone';
        case goGNSS.MODE_PP_LS_CP_SA
            option_i='Post Proc Least Squares on Code and Phase Stand Alone';
        case goGNSS.MODE_PP_LS_CP_VEL
            option_i='Post Proc Least Squares on Code and Phase for Velocity estimation';
        case goGNSS.MODE_PP_LS_C_DD
            option_i='Post Proc Least Squares on Code Double Differences';
        case goGNSS.MODE_PP_LS_CP_DD_L
            option_i='Post Proc Least Squares on Code Double Differences with LAMBDA';
        case goGNSS.MODE_PP_LS_CP_DD_MR
            option_i='Post Proc Least Squares on Code and Phase Double Differences, Multiple Receivers';
        case goGNSS.MODE_PP_LS_C_SA_MR
            option_i='Post Proc Least Squares on Code Stand Alone, Multiple Receivers';
        case goGNSS.MODE_PP_KF_C_SA
            option_i='Post Proc Kalman Filter on Code Stand Alone';
        case goGNSS.MODE_PP_KF_C_DD
            option_i='Post Proc Kalman Filter on Code Double Differences';
        case goGNSS.MODE_PP_KF_CP_SA
            option_i='Post Proc Kalman Filter on Code and Phase Stand Alone';
        case goGNSS.MODE_PP_KF_CP_DD
            option_i='Post Proc Kalman Filter on Code and Phase Double Differences';
        case goGNSS.MODE_PP_KF_CP_DD_MR
            option_i='Post Proc Kalman Filter on Code and Phase Double Differences, Multiple Receivers';
    end
    fprintf(fout_report,'Processing MODE                       : %s\n', option_i);
    fprintf(fout_report,'Constellations enabled                : %s\n', report.opt.constellations);
    if report.opt.flag_SP3 == 1
        option_i = 'PRECISE';
    else
        option_i = 'BROADCAST';
    end
    fprintf(fout_report,'ORBITS type                           : %s\n', option_i);
    
    if goGNSS.isDD(report.opt.mode)
        if report.opt.flag_ms_pos == 1
            option_i = 'YES';
        else
            option_i = 'NO';
        end
        fprintf(fout_report,'Get MASTER STATION position from RINEX: %s\n', option_i);
    end
    
    if report.opt.flag_SBAS == 1
        option_i = 'YES';
    else
        option_i = 'NO';
    end
    fprintf(fout_report,'Apply SBAS corrections                : %s\n', option_i);
    
    if goGNSS.isPH(report.opt.mode)
        if report.opt.flag_IAR == 1
            option_i = 'YES';
        else
            option_i = 'NO';
        end
        fprintf(fout_report,'Solve ambiguities with LAMBDA         : %s\n', option_i);
    end
    fprintf(fout_report,'\n\n');
    
    
    
    fprintf(fout_report,'INPUT FILENAMES\n');
    fprintf(fout_report,'***************\n');
    if is_subfield(report,'inp.iniFile')
        fprintf(fout_report,'INI file                      : %s\n', report.inp.iniFile);
    end
    if goGNSS.isDD(report.opt.mode)
        fprintf(fout_report,'MASTER STATION filename       : %s\n', report.inp.filename_M_obs);
    end
    if ~iscell(report.inp.filename_R_obs)
        fprintf(fout_report,'ROVER STATION filename        : %s\n', report.inp.filename_R_obs);
    else
        fprintf(fout_report,'ROVER STATIONS filename       : ');
        for i=2:length(report.inp.filename_R_obs)
            fprintf(fout_report, '%s%s\n', char(report.inp.filename_R_obs(1)),char(report.inp.filename_R_obs(i)));
            if i < length(report.inp.filename_R_obs)
                fprintf(fout_report, '                                ');
            end
        end
    end
    
    fprintf(fout_report,'NAVIGATIONAL filename         : %s\n', report.inp.filename_nav);
    fprintf(fout_report,'PHASE CENTER OFFSET filename  : %s\n', report.inp.filename_pco);
    if is_subfield(report,'inp.sta_coord_file')
        fprintf(fout_report,'STATION COORDINATES filename  : %s\n', report.inp.sta_coord_file);
    end
    fprintf(fout_report,'\n\n');
    
    fprintf(fout_report,'PROCESSING OPTIONS\n');
    fprintf(fout_report,'******************\n');
    fprintf(fout_report,'Minimum number of observed epochs    : %d \n', report.opt.min_epoch);    
    fprintf(fout_report,'Cut off elevation                    : %.0f degrees\n', report.opt.cutoff);
    if goGNSS.isKM(report.opt.mode)
        fprintf(fout_report,'Min number of sat. for Kalman Filter : %2d\n', report.opt.min_nsat);
    end
    fprintf(fout_report,'Variance of CODE 1 observation       : %7.4f m\n', sqrt(report.opt.sigmaq_cod1));
    fprintf(fout_report,'Variance of CODE 2 observation       : %7.4f m\n', sqrt(report.opt.sigmaq_cod2));
    if goGNSS.isPH(report.opt.mode)
        fprintf(fout_report,'Variance of PHASE observations       : %7.4f m\n', sqrt(report.opt.sigmaq_ph));
        fprintf(fout_report,'Variance of ambiguity combinations   : %d cycles \n', report.opt.sigmaq0_N);
        fprintf(fout_report,'Cycle slip threshold                 : %.2f cycles\n', report.opt.cs_threshold);
    end
    fprintf(fout_report,'Outlier test significance            :  0.995\n');
    switch report.opt.weights
        case 0
            option_i='same weight for all the observations';
        case 1
            option_i='weight based on satellite elevation (sin)';
        case 2
            option_i='weight based on signal-to-noise ratio';
        case 3
            option_i='weight based on combined elevation and signal-to-noise ratio';
        case 4
            option_i='weight based on satellite elevation (exp)';
    end
    fprintf(fout_report,'Weight of observations               : %s\n', option_i);
    
    if goGNSS.isPH(report.opt.mode) && report.opt.flag_IAR == 1
        fprintf(fout_report,'LAMBDA options:\n');
        switch report.opt.IAR_method
            case 0
                option_i='ILS method with numeration in search (LAMBDA2)';
            case 1
                option_i='ILS method with shrinking ellipsoid during search (LAMBDA3)';
            case 2
                option_i='ILS method with numeration in search (LAMBDA3)';
            case 3
                option_i='Integer rounding method (LAMBDA3)';
            case 4
                option_i='Integer bootstrapping method (LAMBDA3)';
            case 5
                option_i='Partial Ambiguity Resolution (PAR) (LAMBDA3)';
        end
        
        fprintf(fout_report,'   -Integer Least Squares estimator  : %s\n', option_i);
        if report.opt.flag_default_P0 == 1
            option_i = 'YES';
        else
            option_i = 'NO';
        end
        fprintf(fout_report,'   -Automatic determination of P0    : %s\n', option_i);
        if report.opt.flag_default_P0 == 0
            if report.opt.flag_IAR==1 || report.opt.flag_IAR==2
                fprintf(fout_report,'   -Fixed failure rate               : %.5d\n', report.opt.P0);
            elseif report.opt.flag_IAR==5
                fprintf(fout_report,'   -Minimum required success rate    : %.2d\n', report.opt.P0);
            end
        end
        fprintf(fout_report,'   -Minimum required success rate    : %.2d\n', report.opt.P0);
        
        if report.opt.flag_auto_mu == 1
            option_i = 'YES';
        else
            option_i = 'NO';
        end
        fprintf(fout_report,'   -Automatic determination of MU    : %s\n', option_i);
        if report.opt.flag_auto_mu == 0
            fprintf(fout_report,'   -Threshold for ratio test (MU)    : %.4f\n', report.opt.mu);
        end
    end
    
    
    % table header
    fprintf(fout_report,'\n\nSTATION INFORMATION\n');
    fprintf(fout_report,'*******************\n\n');
    fprintf(fout_report,'Observations (RAW)              Start time               End time                 Rate  #Sat   #Epoch    #Frq   #C1/P1  #C2/P2     #L1     #L2   #DOP1   #DOP2  %%Epoch %%L2/L1\n');
    fprintf(fout_report,'-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
    for f = 1 : length(report.obs.filename)
        fprintf(fout_report,'%-30s  %23s  %23s  %4.1f    %2d   %6d  %6d  %7d %7d %7d %7d %7d %7d  %6.1f %-6s\n', char(report.obs.filename(f)), char(report.obs_raw.time_start(f)), ...
            char(report.obs_raw.time_end(f)), report.obs_raw.interval(f), report.obs_raw.n_sat(f), report.obs_raw.n_epoch(f), report.obs_raw.nfreq(f), report.obs_raw.n_pr1(f), report.obs_raw.n_pr2(f), ...
            report.obs_raw.n_ph1(f), report.obs_raw.n_ph2(f), report.obs_raw.n_dop1(f), report.obs_raw.n_dop2(f), char(report.obs_raw.epoch_completeness(f)), char(report.obs_raw.L1L2_completeness(f)));
    end
    
    if goGNSS.isDD(report.opt.mode)
        fprintf(fout_report,'\nObservations (after sync)       Start time               End time                 Rate  #Sat   #Epoch    #Frq   #C1/P1  #C2/P2     #L1     #L2   #DOP1   #DOP2  %%Epoch %%L2/L1\n');
        fprintf(fout_report,'-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
        for f = 1 : length(report.obs.filename)
            fprintf(fout_report,'%-30s  %23s  %23s  %4.1f    %2d   %6d  %6d  %7d %7d %7d %7d %7d %7d  %6.1f %-6s\n', char(report.obs.filename(f)), char(report.obs_sync.time_start(f)), ...
                char(report.obs_sync.time_end(f)), report.obs_sync.interval(f), report.obs_sync.n_sat(f), report.obs_sync.n_epoch(f), report.obs_sync.nfreq(f), report.obs_sync.n_pr1(f), report.obs_sync.n_pr2(f), ...
                report.obs_sync.n_ph1(f), report.obs_sync.n_ph2(f), report.obs_sync.n_dop1(f), report.obs_sync.n_dop2(f), char(report.obs_sync.epoch_completeness(f)), char(report.obs_sync.L1L2_completeness(f)));
        end
    end
    
    fprintf(fout_report,'\n\nAntenna information             Antenna model            EAST   NORTH      UP     PCO\n');
    fprintf(fout_report,'-------------------------------------------------------------------------------------\n');
    %if size(report.inp.filename_R_obs,2) == 1 % incase of MR
    fprintf(fout_report,'%-30s  %-20s  %7.4f %7.4f %7.4f %7s\n', char(report.obs.filename(1)), char(report.obs.antname_R), report.obs.antoff_R(1), report.obs.antoff_R(2), report.obs.antoff_R(3), char(report.obs.pcv_yn(1)));
    %else   % this block may work incase of MR
    %    for f = 1 : size(report.inp.filename_R_obs,2)-1
    %        fprintf(fout_report,'%-30s  %-20s  %7.4f %7.4f %7.4f %7s\n', char(report.obs.filename(f)), char(report.obs.antname_R(f)), report.obs.antoff_R(f,1), report.obs.antoff_R(f,2), report.obs.antoff_R(f,3), char(report.obs.pcv_yn(f)));
    %    end
    %end
    if  goGNSS.isDD(report.opt.mode)
        fprintf(fout_report,'%-30s  %-20s  %7.4f %7.4f %7.4f %7s\n', char(report.obs.filename(end)), char(report.obs.antname_M), report.obs.antoff_M(1), report.obs.antoff_M(2), report.obs.antoff_M(3), char(report.obs.pcv_yn(end)));
    end
    
    if is_subfield(report,'obs.coord_R')
        fprintf(fout_report,'\n\nA-priori coordinates                    X             Y             Z      Method                        \n');
        fprintf(fout_report,'-------------------------------------------------------------------------------------------------------------\n');
        fprintf(fout_report,'%s\n',char(report.obs.coord_R));
        if goGNSS.isDD(report.opt.mode)
            fprintf(fout_report,'%s\n',char(report.obs.coord_M));
        end
    end


    fprintf(fout_report,'\n\n\nOBSERVATION PREPROCESSING (RECEIVER CLOCKS SYNCHRONIZATION WITH CODE)\n');
    fprintf(fout_report,'*********************************************************************');
    
    if  report.errors.few_epochs == 1                                               
        fprintf(fout_report,'\n\n     ******* ERROR! THE NUMBER OF AVAILABLE EPOCHS IS LOWER THAN THE THRESHOLD ******\n');
        fprintf(fout_report,'             Number of available epochs: %d\n', report.obs_sync.n_epoch(1));
        fprintf(fout_report,'             Minimum number            : %d\n', report.opt.min_epoch);
        fclose(fout_report);
        return        
    end
    
    fprintf(fout_report,'\n\nRMS threshold for SPP single epoch solution (m) = %.1f\n',report.prep.spp_threshold);
    
    fprintf(fout_report,'\n\nReceiver clock estimation           #Epoch    A-posteriori RMS (m)      #Bad epochs     RMS of worst epoch    Coord.estimated\n');
    fprintf(fout_report,'-----------------------------------------------------------------------------------------------------------------------------\n');
    if report.prep.flag_R==0
        option_i = 'YES';
    else
        option_i = ' NO';
    end
    fprintf(fout_report,'%-30s  %10d %16.2f %14d (%5.1f%%) %16.2f               %s\n', char(report.obs.filename(1)), report.prep.proc_epoch_R(1),  report.prep.varSPP_R(1), ...
        report.prep.bad_epoch_R(1), report.prep.bad_epoch_R(1)/report.prep.tot_epoch_R(1)*100,  report.prep.max_varSPP_R(1), option_i);
    if goGNSS.isDD(report.opt.mode)
        fprintf(fout_report,'%-30s  %10d %16.2f %14d (%5.1f%%) %16.2f               %s\n', char(report.obs.filename(end)), report.prep.proc_epoch_M(1),  report.prep.varSPP_M(1), ...
            report.prep.bad_epoch_M(1), report.prep.bad_epoch_M(1)/report.prep.tot_epoch_M(1)*100,  report.prep.max_varSPP_M(1), ' NO');
    end
    
    
    
    fprintf(fout_report,'\nC1/P1 Observation distribution        #Sat     #Total          #Used              #Cutoff            #Outlier    \n');
    fprintf(fout_report,'-----------------------------------------------------------------------------------------------------------------\n');
    fprintf(fout_report,'%-30s  %10d %10d %10d (%5.1f%%) %10d (%5.1f%%) %10d (%5.1f%%)\n', char(report.obs.filename(1)), sum(report.prep.obs_stat_R(:,2,1)~=0), ...
        report.prep.tot_obs_R(1), report.prep.obs_used_R(1), report.prep.obs_used_R(1)/report.prep.tot_obs_R(1)*100, report.prep.obs_undercutoff_R(1), ...
        report.prep.obs_undercutoff_R(1)/report.prep.tot_obs_R(1)*100, report.prep.obs_outlier_R(1), report.prep.obs_outlier_R(1)/(report.prep.obs_outlier_R(1)+report.prep.obs_used_R(1))*100);
    if goGNSS.isDD(report.opt.mode)
        fprintf(fout_report,'%-30s  %10d %10d %10d (%5.1f%%) %10d (%5.1f%%) %10d (%5.1f%%)\n', char(report.obs.filename(end)), sum(report.prep.obs_stat_M(:,2)~=0), ...
            report.prep.tot_obs_M, report.prep.obs_used_M, report.prep.obs_used_M/report.prep.tot_obs_M*100, report.prep.obs_undercutoff_M, ...
            report.prep.obs_undercutoff_M/report.prep.tot_obs_M*100, report.prep.obs_outlier_M, report.prep.obs_outlier_M/(report.prep.obs_outlier_M+report.prep.obs_used_M)*100);
    end
    

    fprintf(fout_report,'\nC1/P1 Observation by satellite              ');
    for i=1:size(report.prep.obs_stat_R,1)
        fprintf(fout_report,'  %4s',char(report.opt.sat_id(i)));        
    end
    fprintf(fout_report,'\n---------------------------------------------');
    for i=1:size(report.prep.obs_stat_R,1)
        fprintf(fout_report,'------');        
    end
    fprintf(fout_report,'\n%-30s        # Obs ', char(report.obs.filename(1)));
    for i=1:size(report.prep.obs_stat_R,1)
        fprintf(fout_report,'%6d',report.prep.obs_stat_R(i,1,1)+report.prep.obs_stat_R(i,2,1)+report.prep.obs_stat_R(i,3,1));  
    end
        
    fprintf(fout_report,'\n                                      %%Used ');
    for i=1:size(report.prep.obs_stat_R,1)
        if (report.prep.obs_stat_R(i,1,1)+report.prep.obs_stat_R(i,2,1)+report.prep.obs_stat_R(i,3,1)) > 0
            result_i = report.prep.obs_stat_R(i,2,1)/(report.prep.obs_stat_R(i,1,1)+report.prep.obs_stat_R(i,2,1)+report.prep.obs_stat_R(i,3,1))*100;
        else
            result_i = 0;
        end
        fprintf(fout_report,'%6.1f',result_i);  
    end
    fprintf(fout_report,'\n                                      %% Out ');
    for i=1:size(report.prep.obs_stat_R,1)
        if (report.prep.obs_stat_R(i,1,1)+report.prep.obs_stat_R(i,2,1)+report.prep.obs_stat_R(i,3,1)) > 0
            result_i = report.prep.obs_stat_R(i,3,1)/(report.prep.obs_stat_R(i,1,1)+report.prep.obs_stat_R(i,2,1)+report.prep.obs_stat_R(i,3,1))*100;
        else
            result_i = 0;
        end
        fprintf(fout_report,'%6.1f',result_i);  
    end    
    fprintf(fout_report,'\n                                      %% Cut ');
    for i=1:size(report.prep.obs_stat_R,1)
        if report.prep.obs_stat_R(i,1,1)+report.prep.obs_stat_R(i,2,1)+report.prep.obs_stat_R(i,3,1) > 1
            result_i = report.prep.obs_stat_R(i,1,1)/(report.prep.obs_stat_R(i,1,1)+report.prep.obs_stat_R(i,2,1)+report.prep.obs_stat_R(i,3,1))*100;
        else
            result_i = 0;
        end
        fprintf(fout_report,'%6.1f',result_i);
    end
    fprintf(fout_report,'\n');
    
    if goGNSS.isDD(report.opt.mode)
        fprintf(fout_report,'\n%-30s        # Obs ', char(report.obs.filename(end)));
        for i=1:size(report.prep.obs_stat_M,1)
            fprintf(fout_report,'%6d',report.prep.obs_stat_M(i,1)+report.prep.obs_stat_M(i,2)+report.prep.obs_stat_M(i,3));
        end
        
        fprintf(fout_report,'\n                                      %%Used ');
        for i=1:size(report.prep.obs_stat_R,1)
            if report.prep.obs_stat_M(i,1)+report.prep.obs_stat_M(i,2)+report.prep.obs_stat_M(i,3) > 0
                result_i = report.prep.obs_stat_M(i,2)/(report.prep.obs_stat_M(i,1)+report.prep.obs_stat_M(i,2)+report.prep.obs_stat_M(i,3))*100;
            else
                result_i = 0;
            end
            fprintf(fout_report,'%6.1f',result_i);
        end
        fprintf(fout_report,'\n                                      %% Out ');
        for i=1:size(report.prep.obs_stat_R,1)
            if report.prep.obs_stat_M(i,1)+report.prep.obs_stat_M(i,2)+report.prep.obs_stat_M(i,3) > 0
                result_i = report.prep.obs_stat_M(i,3)/(report.prep.obs_stat_M(i,1)+report.prep.obs_stat_M(i,2)+report.prep.obs_stat_M(i,3))*100;
            else
                result_i = 0;
            end
            fprintf(fout_report,'%6.1f',result_i);
        end
        fprintf(fout_report,'\n                                      %% Cut ');
        for i=1:size(report.prep.obs_stat_R,1)
            if report.prep.obs_stat_M(i,1)+report.prep.obs_stat_M(i,2)+report.prep.obs_stat_M(i,3) > 0
                result_i = report.prep.obs_stat_M(i,1)/(report.prep.obs_stat_M(i,1)+report.prep.obs_stat_M(i,2)+report.prep.obs_stat_M(i,3))*100;
            else
                result_i = 0;
            end
            fprintf(fout_report,'%6.1f',result_i);
        end
        fprintf(fout_report,'\n');
    end
    
      
    
    fprintf(fout_report,'\nCycle slips fixing                      #CS    Sat    Frq      Epoch      Float cycles     Fixed?       Fixed cycles\n');
    fprintf(fout_report,'--------------------------------------------------------------------------------------------------------------------\n');
    CS_R = report.prep.CS_R{1};
    fprintf(fout_report,'%-30s %12d   ', char(report.obs.filename(1)), size(CS_R,1));
    if size(CS_R,1) == 0
        fprintf(fout_report,'\n');
    else
        for i = 1 : size(CS_R,1)
            if CS_R(i,6) == 1
                result_i='YES';
            else
                result_i='NO';
            end
            if i>1
                fprintf(fout_report,'                                              ');
            end
            fprintf(fout_report,'%4s      %1d  %9d   %15.3f       %3s     %15.3f\n', char(report.opt.sat_id(CS_R(i,1))), CS_R(i,2), CS_R(i,3), CS_R(i,5), result_i, CS_R(i,4));
        end
    end
    if goGNSS.isDD(report.opt.mode)
        CS_M = report.prep.CS_M;
        fprintf(fout_report,'%-30s %12d   ', char(report.obs.filename(end)), size(CS_M,1));
        if size(CS_M,1) == 0
            fprintf(fout_report,'\n');
        else
            for i = 1 : size(CS_M,1)
                if CS_M(i,6) == 1
                    result_i='YES';
                else
                    result_i='NO';
                end
                if i>1
                    fprintf(fout_report,'                                              ');
                end
                fprintf(fout_report,'%4s      %1d  %9d   %15.3f       %3s     %15.3f\n', char(report.opt.sat_id(CS_M(i,1))), CS_M(i,2), CS_M(i,3), CS_M(i,5), result_i, CS_M(i,4));
            end
        end
        
    end
    
    

    fclose(fout_report);
    
    
end
end



function issbf=is_subfield(s,field)
issbf = 1;
argv_com = ['s','.',field,';'];
try eval(argv_com);
catch
    issbf = 0;
end
end


