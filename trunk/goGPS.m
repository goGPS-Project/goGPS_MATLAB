%--------------------------*--. --- --. .--. ...*---------------------------------------------
%
%                    %%%%%  %%%%%   %%%%%  %%%%% %%%%%
%                    %      %   %   %      %   % %
%                    %  %%% %   %   %  %%% %%%%% %%%%%
%                    %   %  %   %   %   %  %         %
%                    %%%%%  %%%%%   %%%%%  %     %%%%%
%
%
%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Media Center, Osaka City University, Japan
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
%---------------------------------------------------------------------------------------------

% clear all the variables in the workspace
clear all

% close all windows
close all

% clear the command prompt
%clc

% close all the opened files
fclose('all');

% disable warnings
warning off;

% start evaluating computation time
tic

%----------------------------------------------------------------------------------------------
% INTERFACE TYPE DEFINITION
%----------------------------------------------------------------------------------------------

mode_user = 1;  % user interface type
                % mode_user=0 --> use text interface
                % mode_user=1 --> use GUI

%----------------------------------------------------------------------------------------------
% INTERFACE STARTUP
%----------------------------------------------------------------------------------------------

%initialization of global variables/constants
global_init;

global o1 o2 o3 h_antenna

if (mode_user == 1)

    [mode, mode_vinc, mode_data, mode_ref, flag_ms_rtcm, flag_ms, flag_ge, flag_cov, flag_COM, flag_NTRIP, flag_amb, ...
        flag_cpu, filerootIN, filerootOUT, filename_R_obs, filename_R_nav, filename_M_obs, filename_M_nav, filename_ref, ...
        pos_M] = goGPS_gui;

    if (isempty(mode))
        return
    end
else
    
    %-------------------------------------------------------------------------------------------
    % DEFINITION OF THE FUNCTIONING MODE (TEXTUAL INTERFACE)
    %-------------------------------------------------------------------------------------------
    
    mode =   11;     % functioning mode
                     % POST-PROCESSING
                     % mode=1  --> KALMAN FILTER ON PHASE AND CODE DOUBLE DIFFERENCES WITH/WITHOUT A CONSTRAINT
                     % mode=2  --> KALMAN FILTER ON PHASE AND CODE, WITHOUT INTERNET CONNECTION AND WITHOUT A CONSTRAINT (to be implemented)
                     % mode=3  --> LEAST SQUARES ADJ. ON CODE DOUBLE DIFFERENCES, NO CONSTRAINT
                     % mode=4  --> LEAST SQUARES ADJ. ON CODE, NO CONSTRAINT
                     % mode=5  --> KALMAN FILTER ON CODE DOUBLE DIFFERENCES, NO CONSTRAINT
                     % mode=6  --> KALMAN FILTER ON CODE, NO CONSTRAINT
                     % mode=7  --> ....
                     % mode=8  --> ....
                     % mode=9  --> ....
                     % REAL-TIME
                     % mode=11 --> KALMAN FILTER ON PHASE AND CODE DOUBLE DIFFERENCES WITH/WITHOUT A CONSTRAINT
                     % mode=12 --> U-BLOX MONITORING
                     % mode=13 --> MASTER MONITORING
    
    mode_vinc = 0;   % navigation mode
                     % mode_vinc=0 --> without linear constraint
                     % mode_vinc=1 --> with linear constraint
    
    mode_data = 1;   % data loading mode
                     % mode_data=0 --> RINEX data
                     % mode_data=1 --> goGPS data (saved during a real-time session)
                     % mode_data=2 --> stream data (saved during a real-time session)
    
    mode_ref = 0;    % reference path mode
                     % mode_ref=0 --> do not use a reference path
                     % mode_ref=1 --> use a reference path (plot it and use it for statistics)
    
    flag_ms_rtcm = 1;% read master station position from RTCM
    
    flag_ms = 0;     % plot master station position --> no=0, yes=1
    
    flag_ge = 0;     % use google earth --> no=0, yes=1
    
    flag_cov = 1;    % plot error ellipse --> no=0, yes=1
    
    flag_COM = 0;    % u-blox COM automatic detection --> no=0, yes=1
    
    flag_NTRIP = 1;  % use NTRIP --> no=0, yes=1
    
    flag_amb = 0;    % plot ambiguities (only in post-processing)
    
    flag_cpu = 0;    % save CPU (don't draw skyplot and SNR graph)
    
    %----------------------------------------------------------------------------------------------
    % USER-DEFINED SETTINGS
    %----------------------------------------------------------------------------------------------
    
    %User-defined global settings
    global_settings;
    
end

%-------------------------------------------------------------------------------------------
% REFERENCE PATH LOAD
%-------------------------------------------------------------------------------------------

if (mode_ref == 1)
    
    d = dir(filename_ref);
    
    if ~isempty(d)
        load(filename_ref, 'ref_path', 'mat_path');
        
        %adjust the reference path according to antenna height
        [ref_phi, ref_lam, ref_h] = cart2geod(ref_path(:,1),ref_path(:,2),ref_path(:,3));
        ref_h = ref_h + h_antenna;
        [ref_X, ref_Y, ref_Z] = geod2cart(ref_phi, ref_lam, ref_h, a, f);
        ref_path = [ref_X , ref_Y , ref_Z];
        
    else
        ref_path = [];
        mat_path = [];
    end
    
else
    ref_path = [];
    mat_path = [];
end

%----------------------------------------------------------------------------------------------
% FILE READING
%----------------------------------------------------------------------------------------------

if (mode < 10) %post-processing
    
    if (mode_data == 0)
        
        %read data from RINEX files
        [pr1_R, pr1_M, ph1_R, ph1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
            Eph_R, Eph_M, iono_R, iono_M, snr_R, snr_M, ...
            pr1_RR, pr1_MR, ph1_RR, ph1_MR, pr2_RR, pr2_MR, ph2_RR, ph2_MR, ...
            Eph_RR, Eph_MR, snr_RR, snr_MR, ...
            time_GPS, date] = ...
            load_RINEX(filename_R_obs, filename_R_nav, filename_M_obs, filename_M_nav);

        %select ephemerides source
        %Eph = Eph_R;
        Eph = Eph_M;
        %Eph_GLO = Eph_MR;
        
        %convert signal-to-noise ratio
        snr_R = 6 * snr_R;
        snr_M = 6 * snr_M;
        %snr_RR = 6 * snr_RR;
        %snr_MR = 6 * snr_MR;
        
        %remove satellites without ephemerides (GPS)
        delsat = setdiff(1:32,unique(Eph(1,:)));
        pr1_R(delsat,:) = 0;
        pr1_M(delsat,:) = 0;
        pr2_R(delsat,:) = 0;
        pr2_M(delsat,:) = 0;
        ph1_R(delsat,:) = 0;
        ph1_M(delsat,:) = 0;
        ph2_R(delsat,:) = 0;
        ph2_M(delsat,:) = 0;
        snr_R(delsat,:) = 0;
        snr_M(delsat,:) = 0;
        
        %%remove satellites without ephemerides (GLONASS)
        %delsat = setdiff(1:32,unique(Eph_GLO(1,:)));
        %pr1_RR(delsat,:) = 0;
        %pr1_MR(delsat,:) = 0;
        %pr2_RR(delsat,:) = 0;
        %pr2_MR(delsat,:) = 0;
        %ph1_RR(delsat,:) = 0;
        %ph1_MR(delsat,:) = 0;
        %ph2_RR(delsat,:) = 0;
        %ph2_MR(delsat,:) = 0;
        %snr_RR(delsat,:) = 0;
        %snr_MR(delsat,:) = 0;
        
        %%reverse the path (GPS)
        %pr1_R = pr1_R(:,end:-1:1);
        %pr1_M = pr1_M(:,end:-1:1);
        %ph1_R = ph1_R(:,end:-1:1);
        %ph1_M = ph1_M(:,end:-1:1);
        %pr2_R = pr2_R(:,end:-1:1);
        %pr2_M = pr2_M(:,end:-1:1);
        %ph2_R = ph2_R(:,end:-1:1);
        %ph2_M = ph2_M(:,end:-1:1);
        %snr_R = snr_R(:,end:-1:1);
        %snr_M = snr_M(:,end:-1:1);
        
        %%reverse the path (GLONASS)
        %pr1_RR = pr1_RR(:,end:-1:1);
        %pr1_MR = pr1_MR(:,end:-1:1);
        %ph1_RR = ph1_RR(:,end:-1:1);
        %ph1_MR = ph1_MR(:,end:-1:1);
        %pr2_RR = pr2_RR(:,end:-1:1);
        %pr2_MR = pr2_MR(:,end:-1:1);
        %ph2_RR = ph2_RR(:,end:-1:1);
        %ph2_MR = ph2_MR(:,end:-1:1);
        %snr_RR = snr_RR(:,end:-1:1);
        %snr_MR = snr_MR(:,end:-1:1);
        
        %time_GPS = time_GPS(end:-1:1);
        %date = date(end:-1:1,:);
        
    elseif (mode_data == 1)
        
        %read data from goGPS saved files
        [time_GPS, time_R, time_M, pr1_R, pr1_M, ph1_R, ph1_M, snr_R, snr_M, ...
            pos_M, Eph, delay, loss_R, loss_M] = load_goGPSinput (filerootIN);

        %reference GPS time
        time_GPS = time_GPS(1) + (0 : 1 : length(time_M)-1)';
        
        %nominal date for Google Earth
        date0 = datenum(2009,1,1,0,0,0);
        date = datevec(date0 + time_GPS/86400);
        
        %other variables
        iono_R = zeros(8,1);
        iono_M = zeros(8,1);
        pr2_M = zeros(size(pr1_M));
        pr2_R = zeros(size(pr1_R));
        ph2_M = zeros(size(ph1_M));
        ph2_R = zeros(size(ph1_R));
        
        %complete/partial path
        %tMin = 689;
        %tMax = 1466;
        tMin = 1;
        tMax = 1e30;
        tMin = max(tMin,1);
        tMax = min(tMax,length(time_GPS));
        time_GPS = time_GPS(tMin:tMax);
        time_R = time_R(tMin:tMax);
        time_M = time_M(tMin:tMax);
        pr1_R = pr1_R(:,tMin:tMax);
        pr1_M = pr1_M(:,tMin:tMax);
        ph1_R = ph1_R(:,tMin:tMax);
        ph1_M = ph1_M(:,tMin:tMax);
        pr2_R = pr2_R(:,tMin:tMax);
        pr2_M = pr2_M(:,tMin:tMax);
        ph2_R = ph2_R(:,tMin:tMax);
        ph2_M = ph2_M(:,tMin:tMax);
        snr_R = snr_R(:,tMin:tMax);
        snr_M = snr_M(:,tMin:tMax);
        pos_M = pos_M(:,tMin:tMax);
        Eph = Eph(:,:,tMin:tMax);
        delay = delay(tMin:tMax);
        loss_R = loss_R(tMin:tMax);
        loss_M = loss_M(tMin:tMax);
        date = date(tMin:tMax,:);
        
    else %if (mode_data == 2)
        
        %read data from saved stream files
        [time_GPS, time_R, time_M, pr1_R, pr1_M, ph1_R, ph1_M, snr_R, snr_M, pos_M, ...
            Eph, loss_R, loss_M, stream_R, stream_M] = load_stream (filerootIN);
        satObs_R = find( (pr1_R(:,1) ~= 0) & (ph1_R(:,1) ~= 0) );
        satObs_M = find( (pr1_M(:,1) ~= 0) & (ph1_M(:,1) ~= 0) );
        satEph = find(sum(abs(Eph(:,:,1)))~=0);
        satObs = find( (pr1_R(:,1) ~= 0) & (pr1_M(:,1) ~= 0));
        while (length(satEph) < length(satObs_M)) | (length(satObs) < 4)
            
            time_GPS(1) = [];
            time_R(1) = [];
            time_M(1) = [];
            pr1_R(:,1) = [];
            pr1_M(:,1) = [];
            ph1_R(:,1) = [];
            ph1_M(:,1) = [];
            snr_R(:,1) = [];
            snr_M(:,1) = [];
            pos_M(:,1) = [];
            Eph(:,:,1) = [];
            loss_R(1) = [];
            loss_M(1) = [];
            
            satObs_R = find( (pr1_R(:,1) ~= 0) & (ph1_R(:,1) ~= 0) );
            satObs_M = find( (pr1_M(:,1) ~= 0) & (ph1_M(:,1) ~= 0) );
            satObs = find( (pr1_R(:,1) ~= 0) & (pr1_M(:,1) ~= 0));
            satEph = find(sum(abs(Eph(:,:,1)))~=0);
        end

        %remove observations without ephemerides
        for i = 1 : length(time_GPS)
            satEph = find(sum(abs(Eph(:,:,i)))~=0);
            delsat = setdiff(1:32,satEph);
            pr1_R(delsat,i) = 0;
            pr1_M(delsat,i) = 0;
            ph1_R(delsat,i) = 0;
            ph1_M(delsat,i) = 0;
            snr_R(delsat,i) = 0;
            snr_M(delsat,i) = 0;
        end
        
        %nominal date for Google Earth
        date0 = datenum(2009,1,1,0,0,0);
        date = datevec(date0 + time_GPS/86400);
        
        %other variables
        iono_R = zeros(8,1);
        iono_M = zeros(8,1);
        pr2_M = zeros(size(pr1_M));
        pr2_R = zeros(size(pr1_R));
        ph2_M = zeros(size(ph1_M));
        ph2_R = zeros(size(ph1_R));
        
        %complete/partial path
        tMin = 1;
        tMax = 190e30;
        tMin = max(tMin,1);
        tMax = min(tMax,length(time_GPS));
        time_GPS = time_GPS(tMin:tMax);
        time_R = time_R(tMin:tMax);
        time_M = time_M(tMin:tMax);
        pr1_R = pr1_R(:,tMin:tMax);
        pr1_M = pr1_M(:,tMin:tMax);
        ph1_R = ph1_R(:,tMin:tMax);
        ph1_M = ph1_M(:,tMin:tMax);
        pr2_R = pr2_R(:,tMin:tMax);
        pr2_M = pr2_M(:,tMin:tMax);
        ph2_R = ph2_R(:,tMin:tMax);
        ph2_M = ph2_M(:,tMin:tMax);
        snr_R = snr_R(:,tMin:tMax);
        snr_M = snr_M(:,tMin:tMax);
        pos_M = pos_M(:,tMin:tMax);
        Eph = Eph(:,:,tMin:tMax);
        date = date(tMin:tMax,:);
        loss_R = loss_R(tMin:tMax);
        loss_M = loss_M(tMin:tMax);
        
    end
    
else %real-time
    
    %for the Kalman filter execution in real-time
    iono_R = zeros(8,1);
    iono_M = zeros(8,1);
    pr2_M = zeros(32,1);
    pr2_R = zeros(32,1);
    ph2_M = zeros(32,1);
    ph2_R = zeros(32,1);
    
end

%check if MASTER position is available
if (mode < 12) & (~flag_ms_rtcm)
    pos_M(1,1:length(time_GPS)) = pos_M(1);
    pos_M(2,1:length(time_GPS)) = pos_M(2);
    pos_M(3,1:length(time_GPS)) = pos_M(3);
elseif (mode < 10) & (sum(abs(pos_M)) == 0 | isempty(pos_M))
    pos_M(1,1:length(time_GPS)) = pos_M(1);
    pos_M(2,1:length(time_GPS)) = pos_M(2);
    pos_M(3,1:length(time_GPS)) = pos_M(3);
    fprintf('Master position not found, thus it is set to user-defined values:\n');
    fprintf(' X=%.4f, Y=%.4f, Z=%.4f km\n', pos_M(1,1)/1000, pos_M(2,1)/1000, pos_M(3,1)/1000);
end

%----------------------------------------------------------------------------------------------
% POST-PROCESSING: KALMAN FILTER ON PHASE AND CODE DOUBLE DIFFERENCES WITHOUT A CONSTRAINT
%----------------------------------------------------------------------------------------------

if (mode == 1) & (mode_vinc == 0)
    
    fid_kal = fopen([filerootOUT '_kal_00.bin'],'w+');
    fid_sat = fopen([filerootOUT '_sat_00.bin'],'w+');
    fid_conf = fopen([filerootOUT '_conf_00.bin'],'w+');
    
    if (mode_data == 0)
        Eph_t = rt_find_eph (Eph, time_GPS(1));
    else
        Eph_t = Eph(:,:,1);
    end
    
    kalman_goGPS_init (pos_M(:,1), time_GPS(1), Eph_t, iono_R, pr1_R(:,1), pr1_M(:,1), ph1_R(:,1), ph1_M(:,1), pr2_R(:,1), pr2_M(:,1), ph2_R(:,1), ph2_M(:,1), snr_R(:,1), snr_M(:,1), 1);
    
    fwrite(fid_kal, [Xhat_t_t; Cee(:)], 'double');
    fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
    fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');

    if (flag_cov == 0)
        if (flag_ge == 1), rtplot_googleearth (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), date(1,:)), end;
        rtplot_matlab (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), 0, 0, 0, 0, flag_ms, ref_path, mat_path, flag_amb);
    else
        if (flag_ge == 1), rtplot_googleearth_cov (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date(1,:)), end;
        rtplot_matlab_cov (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), 0, 0, 0, 0, flag_ms, ref_path, mat_path, flag_amb);
    end
    if (flag_amb == 1)
        rtplot_amb (1, window, Xhat_t_t(o3+1:o3+32), sqrt(diag(Cee(o3+1:o3+32,o3+1:o3+32))), conf_cs)
    else
        if (flag_cpu == 0)
            rtplot_skyplot (1, azR, elR, conf_sat, pivot);
            rtplot_snr (snr_R(:,1));
        else
            rttext_sat (1, azR, elR, snr_R(:,1), conf_sat, pivot);
        end
    end
    
    for t = 2 : length(time_GPS)
        
        if (mode_data == 0)
            Eph_t = rt_find_eph (Eph, time_GPS(t));
        else
            Eph_t = Eph(:,:,t);
        end
        
        [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_loop (pos_M(:,t), time_GPS(t), Eph_t, iono_R, pr1_R(:,t), pr1_M(:,t), ph1_R(:,t), ph1_M(:,t), pr2_R(:,t), pr2_M(:,t), ph2_R(:,t), ph2_M(:,t), snr_R(:,t), snr_M(:,t), 1);
        
        fwrite(fid_kal, [Xhat_t_t; Cee(:)], 'double');
        fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
        fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');
        
        if (flag_cov == 0)
            if (flag_ge == 1), rtplot_googleearth (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), date(t,:)), end;
            rtplot_matlab (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), check_on, check_off, check_pivot, check_cs, flag_ms, ref_path, mat_path, flag_amb);
        else
            if (flag_ge == 1), rtplot_googleearth_cov (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date(t,:)), end;
            rtplot_matlab_cov (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), check_on, check_off, check_pivot, check_cs, flag_ms, ref_path, mat_path, flag_amb);
        end
        if (flag_amb == 1)
            rtplot_amb (t, window, Xhat_t_t(o3+1:o3+32), sqrt(diag(Cee(o3+1:o3+32,o3+1:o3+32))), conf_cs);
            pause(0.1);
        else
            if (flag_cpu == 0)
                rtplot_skyplot (t, azR, elR, conf_sat, pivot);
                rtplot_snr (snr_R(:,t));
            else
                rttext_sat (t, azR, elR, snr_R(:,t), conf_sat, pivot);
            end
            pause(0.01);
        end

    end
    
    fclose(fid_kal);
    fclose(fid_sat);
    fclose(fid_conf);
    
%----------------------------------------------------------------------------------------------
% POST-PROCESSING: KALMAN FILTER ON PHASE AND CODE DOUBLE DIFFERENCES WITH A CONSTRAINT
%----------------------------------------------------------------------------------------------
    
elseif (mode == 1) & (mode_vinc == 1)
    
    fid_kal = fopen([filerootOUT '_kal_00.bin'],'w+');
    fid_sat = fopen([filerootOUT '_sat_00.bin'],'w+');
    fid_conf = fopen([filerootOUT '_conf_00.bin'],'w+');
    
    if (mode_data == 0)
        Eph_t = rt_find_eph (Eph, time_GPS(1));
    else
        Eph_t = Eph(:,:,1);
    end
    
    %repeat more than once the reference loop
    %(this constrained mode works only for circuits)
    ref_loop = [ref_path; ref_path];
    
    kalman_goGPS_vinc_init (pos_M(:,1), time_GPS(1), Eph_t, iono_R, pr1_R(:,1), pr1_M(:,1), ph1_R(:,1), ph1_M(:,1), pr2_R(:,1), pr2_M(:,1), ph2_R(:,1), ph2_M(:,1), snr_R(:,1), snr_M(:,1), 1, ref_loop);
    
    fwrite(fid_kal, [Xhat_t_t; Yhat_t_t; Cee(:)], 'double');
    fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
    fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');
    
    if (flag_ge == 1), rtplot_googleearth (1, [Yhat_t_t(1); Yhat_t_t(2); Yhat_t_t(3)], date(1,:)), end;
    rtplot_matlab (1, [Yhat_t_t(1); Yhat_t_t(2); Yhat_t_t(3)], pos_M(:,1), 0, 0, 0, 0, flag_ms, ref_path, mat_path, flag_amb);
    if (flag_amb == 1)
        rtplot_amb (1, window, Xhat_t_t(o1+1:o1+32), sqrt(diag(Cee(o1+1:o1+32,o1+1:o1+32))), conf_cs);
    else
        if (flag_cpu == 0)
            rtplot_skyplot (1, azR, elR, conf_sat, pivot);
            rtplot_snr (snr_R(:,1));
        else
            rttext_sat (1, azR, elR, snr_R(:,1), conf_sat, pivot);
        end
    end
    
    for t = 2 : length(time_GPS)
        
        if (mode_data == 0)
            Eph_t = rt_find_eph (Eph, time_GPS(t));
        else
            Eph_t = Eph(:,:,t);
        end
        
        [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_vinc_loop (pos_M(:,t), time_GPS(t), Eph_t, iono_R, pr1_R(:,t), pr1_M(:,t), ph1_R(:,t), ph1_M(:,t), pr2_R(:,t), pr2_M(:,t), ph2_R(:,t), ph2_M(:,t), snr_R(:,t), snr_M(:,t), 1, ref_loop);
        
        fwrite(fid_kal, [Xhat_t_t; Yhat_t_t; Cee(:)], 'double');
        fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
        fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');
        
        if (flag_ge == 1), rtplot_googleearth (t, [Yhat_t_t(1); Yhat_t_t(2); Yhat_t_t(3)], date(t,:)), end;
        rtplot_matlab (t, [Yhat_t_t(1); Yhat_t_t(2); Yhat_t_t(3)], pos_M(:,t), check_on, check_off, check_pivot, check_cs, flag_ms, ref_path, mat_path, flag_amb);
        if (flag_amb == 1)
            rtplot_amb (t, window, Xhat_t_t(o1+1:o1+32), sqrt(diag(Cee(o1+1:o1+32,o1+1:o1+32))), conf_cs);
            pause(0.1);
        else
            if (flag_cpu == 0)
                rtplot_skyplot (t, azR, elR, conf_sat, pivot);
                rtplot_snr (snr_R(:,t));
            else
                rttext_sat (t, azR, elR, snr_R(:,t), conf_sat, pivot);
            end
            pause(0.01);
        end

    end
    
    fclose(fid_kal);
    fclose(fid_sat);
    fclose(fid_conf);
    
%----------------------------------------------------------------------------------------------
% POST-PROCESSING: KALMAN FILTER ON PHASE AND CODE, WITHOUT INTERNET CONNECTION AND WITHOUT A CONSTRAINT
%----------------------------------------------------------------------------------------------
    
elseif (mode == 2)
    
% TO BE IMPLEMENTED
    
%----------------------------------------------------------------------------------------------
% POST-PROCESSING: LEAST SQUARES ADJ. ON CODE DOUBLE DIFFERENCES, NO CONSTRAINT
%----------------------------------------------------------------------------------------------
    
elseif (mode == 3)
    
    fid_kal = fopen([filerootOUT '_kal_00.bin'],'w+');
    fid_sat = fopen([filerootOUT '_sat_00.bin'],'w+');
    fid_conf = fopen([filerootOUT '_conf_00.bin'],'w+');
    
    nN = 32;
    check_on = 0;
    check_off = 0;
    check_pivot = 0;
    check_cs = 0;
    
    if (mode_data == 0)
        Eph_t = rt_find_eph (Eph, time_GPS(1));
    else
        Eph_t = Eph(:,:,1);
    end
    
    MQ_goGPS_loop (time_GPS(1), Eph_t, pos_M(:,1), pr1_R(:,1), pr1_M(:,1), pr2_R(:,1), pr2_M(:,1), snr_R(:,1), snr_M(:,1), 1);
    
    Xhat_t_t_dummy = [Xhat_t_t; zeros(nN,1)];
    Cee_dummy = [Cee zeros(o3,nN); zeros(nN,o3) zeros(nN,nN)];
    fwrite(fid_kal, [Xhat_t_t_dummy; Cee_dummy(:)], 'double');
    fwrite(fid_sat, zeros(32,6), 'double');
    fwrite(fid_conf, zeros(65,1), 'int8');
    
    if (flag_cov == 0)
        if (flag_ge == 1), rtplot_googleearth (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), date(1,:)), end;
        rtplot_matlab (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), 0, 0, 0, 0, flag_ms, ref_path, mat_path);
    else
        if (flag_ge == 1), rtplot_googleearth_cov (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date(1,:)), end;
        rtplot_matlab_cov (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), 0, 0, 0, 0, flag_ms, ref_path, mat_path);
    end
    if (flag_cpu == 0)
        rtplot_skyplot (1, azR, elR, conf_sat, pivot);
        rtplot_snr (snr_R(:,1));
    else
        rttext_sat (1, azR, elR, snr_R(:,1), conf_sat, pivot);
    end
    
    for t = 2 : length(time_GPS)
        
        if (mode_data == 0)
            Eph_t = rt_find_eph (Eph, time_GPS(t));
        else
            Eph_t = Eph(:,:,t);
        end
        
        MQ_goGPS_loop (time_GPS(t), Eph_t, pos_M(:,t), pr1_R(:,t), pr1_M(:,t), pr2_R(:,t), pr2_M(:,t), snr_R(:,t), snr_M(:,t), 1);
        
        Xhat_t_t_dummy = [Xhat_t_t; zeros(nN,1)];
        Cee_dummy = [Cee zeros(o3,nN); zeros(nN,o3) zeros(nN,nN)];
        fwrite(fid_kal, [Xhat_t_t_dummy; Cee_dummy(:)], 'double');
        fwrite(fid_sat, zeros(32,6), 'double');
        fwrite(fid_conf, zeros(65,1), 'int8');
        
        if (flag_cov == 0)
            if (flag_ge == 1), rtplot_googleearth (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), date(t,:)), end;
            rtplot_matlab (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), check_on, check_off, check_pivot, check_cs, flag_ms, ref_path, mat_path);
        else
            if (flag_ge == 1), rtplot_googleearth_cov (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date(t,:)), end;
            rtplot_matlab_cov (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), check_on, check_off, check_pivot, check_cs, flag_ms, ref_path, mat_path);
        end
        if (flag_cpu == 0)
            rtplot_skyplot (t, azR, elR, conf_sat, pivot);
            rtplot_snr (snr_R(:,t));
        else
            rttext_sat (t, azR, elR, snr_R(:,t), conf_sat, pivot);
        end
        
        pause(0.01);
    end
    
    fclose(fid_kal);
    fclose(fid_sat);
    fclose(fid_conf);
    
%----------------------------------------------------------------------------------------------
% POST-PROCESSING: LEAST SQUARES ADJ. ON CODE, NO CONSTRAINT
%----------------------------------------------------------------------------------------------
    
elseif (mode == 4)
    
    fid_kal = fopen([filerootOUT '_kal_00.bin'],'w+');
    fid_sat = fopen([filerootOUT '_sat_00.bin'],'w+');
    fid_conf = fopen([filerootOUT '_conf_00.bin'],'w+');
    
    nN = 32;
    check_on = 0;
    check_off = 0;
    check_pivot = 0;
    check_cs = 0;
    
    if (mode_data == 0)
        Eph_t = rt_find_eph (Eph, time_GPS(1));
    else
        Eph_t = Eph(:,:,1);
    end
    
    MQ_goGPS_SA_loop (time_GPS(1), Eph_t, pr1_R(:,1), pr2_R(:,1), 1);
    
    Xhat_t_t_dummy = [Xhat_t_t; zeros(nN,1)];
    Cee_dummy = [Cee zeros(o3,nN); zeros(nN,o3) zeros(nN,nN)];
    fwrite(fid_kal, [Xhat_t_t_dummy; Cee_dummy(:)], 'double');
    fwrite(fid_sat, zeros(32,6), 'double');
    fwrite(fid_conf, zeros(65,1), 'int8');
    
    if (flag_cov == 0)
        if (flag_ge == 1), rtplot_googleearth (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), date(1,:)), end;
        rtplot_matlab (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), 0, 0, 0, 0, flag_ms, ref_path, mat_path);
    else
        if (flag_ge == 1), rtplot_googleearth_cov (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date(1,:)), end;
        rtplot_matlab_cov (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), 0, 0, 0, 0, flag_ms, ref_path, mat_path);
    end
    if (flag_cpu == 0)
        rtplot_skyplot (1, azR, elR, conf_sat, pivot);
        rtplot_snr (snr_R(:,1));
    else
        rttext_sat (1, azR, elR, snr_R(:,1), conf_sat, pivot);
    end
    
    for t = 2 : length(time_GPS)
        
        if (mode_data == 0)
            Eph_t = rt_find_eph (Eph, time_GPS(t));
        else
            Eph_t = Eph(:,:,t);
        end
        
        MQ_goGPS_SA_loop (time_GPS(t), Eph_t, pr1_R(:,t), pr2_R(:,t), 1);
        
        Xhat_t_t_dummy = [Xhat_t_t; zeros(nN,1)];
        Cee_dummy = [Cee zeros(o3,nN); zeros(nN,o3) zeros(nN,nN)];
        fwrite(fid_kal, [Xhat_t_t_dummy; Cee_dummy(:)], 'double');
        fwrite(fid_sat, zeros(32,6), 'double');
        fwrite(fid_conf, zeros(65,1), 'int8');
        
        if (flag_cov == 0)
            if (flag_ge == 1), rtplot_googleearth (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), date(t,:)), end;
            rtplot_matlab (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), check_on, check_off, check_pivot, check_cs, flag_ms, ref_path, mat_path);
        else
            if (flag_ge == 1), rtplot_googleearth_cov (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date(t,:)), end;
            rtplot_matlab_cov (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), check_on, check_off, check_pivot, check_cs, flag_ms, ref_path, mat_path);
        end
        if (flag_cpu == 0)
            rtplot_skyplot (t, azR, elR, conf_sat, pivot);
            rtplot_snr (snr_R(:,t));
        else
            rttext_sat (t, azR, elR, snr_R(:,t), conf_sat, pivot);
        end
        
        pause(0.01);
    end
    
    fclose(fid_kal);
    fclose(fid_sat);
    fclose(fid_conf);
    
%----------------------------------------------------------------------------------------------
% POST-PROCESSING: KALMAN FILTER ON CODE DOUBLE DIFFERENCES, NO CONSTRAINT
%----------------------------------------------------------------------------------------------
    
elseif (mode == 5)
    
    fid_kal = fopen([filerootOUT '_kal_00.bin'],'w+');
    fid_sat = fopen([filerootOUT '_sat_00.bin'],'w+');
    fid_conf = fopen([filerootOUT '_conf_00.bin'],'w+');
    
    nN = 32;
    check_on = 0;
    check_off = 0;
    check_pivot = 0;
    check_cs = 0;
    
    if (mode_data == 0)
        Eph_t = rt_find_eph (Eph, time_GPS(1));
    else
        Eph_t = Eph(:,:,1);
    end
    
    kalman_goGPS_cod_init (pos_M(:,1), time_GPS(1), Eph_t, iono_R, pr1_R(:,1), pr1_M(:,1), pr2_R(:,1), pr2_M(:,1), snr_R(:,1), snr_M(:,1), 1);
    
    Xhat_t_t_dummy = [Xhat_t_t; zeros(nN,1)];
    Cee_dummy = [Cee zeros(o3,nN); zeros(nN,o3) zeros(nN,nN)];
    fwrite(fid_kal, [Xhat_t_t_dummy; Cee_dummy(:)], 'double');
    fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
    fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');
    
    if (flag_cov == 0)
        if (flag_ge == 1), rtplot_googleearth (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), date(1,:)), end;
        rtplot_matlab (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), 0, 0, 0, 0, flag_ms, ref_path, mat_path);
    else
        if (flag_ge == 1), rtplot_googleearth_cov (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date(1,:)), end;
        rtplot_matlab_cov (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), 0, 0, 0, 0, flag_ms, ref_path, mat_path);
    end
    if (flag_cpu == 0)
        rtplot_skyplot (1, azR, elR, conf_sat, pivot);
        rtplot_snr (snr_R(:,1));
    else
        rttext_sat (1, azR, elR, snr_R(:,1), conf_sat, pivot);
    end
    
    for t = 2 : length(time_GPS)
        
        if (mode_data == 0)
            Eph_t = rt_find_eph (Eph, time_GPS(t));
        else
            Eph_t = Eph(:,:,t);
        end
        
        [check_on, check_off, check_pivot, check_cs] = kalman_goGPS_cod_loop (pos_M(:,t), time_GPS(t), Eph_t, iono_R, pr1_R(:,t), pr1_M(:,t), pr2_R(:,t), pr2_M(:,t), snr_R(:,t), snr_M(:,t), 1);
        
        Xhat_t_t_dummy = [Xhat_t_t; zeros(nN,1)];
        Cee_dummy = [Cee zeros(o3,nN); zeros(nN,o3) zeros(nN,nN)];
        fwrite(fid_kal, [Xhat_t_t_dummy; Cee_dummy(:)], 'double');
        fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
        fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');
        
        if (flag_cov == 0)
            if (flag_ge == 1), rtplot_googleearth (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), date(t,:)), end;
            rtplot_matlab (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), check_on, check_off, check_pivot, check_cs, flag_ms, ref_path, mat_path);
        else
            if (flag_ge == 1), rtplot_googleearth_cov (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date(t,:)), end;
            rtplot_matlab_cov (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), check_on, check_off, check_pivot, check_cs, flag_ms, ref_path, mat_path);
        end
        if (flag_cpu == 0)
            rtplot_skyplot (t, azR, elR, conf_sat, pivot);
            rtplot_snr (snr_R(:,t));
        else
            rttext_sat (t, azR, elR, snr_R(:,t), conf_sat, pivot);
        end
        
        pause(0.01);
    end
    
    fclose(fid_kal);
    fclose(fid_sat);
    fclose(fid_conf);
    
%----------------------------------------------------------------------------------------------
% POST-PROCESSING: KALMAN FILTER ON CODE, NO CONSTRAINT
%----------------------------------------------------------------------------------------------
    
elseif (mode == 6)
    
    fid_kal = fopen([filerootOUT '_kal_00.bin'],'w+');
    fid_sat = fopen([filerootOUT '_sat_00.bin'],'w+');
    fid_conf = fopen([filerootOUT '_conf_00.bin'],'w+');
    
    nN = 32;
    check_on = 0;
    check_off = 0;
    check_pivot = 0;
    check_cs = 0;
    
    if (mode_data == 0)
        Eph_t = rt_find_eph (Eph, time_GPS(1));
    else
        Eph_t = Eph(:,:,1);
    end
    
    kalman_goGPS_SA_cod_init (time_GPS(1), Eph_t, iono_R, pr1_R(:,1), pr2_R(:,1), 1);
    
    Xhat_t_t_dummy = [Xhat_t_t; zeros(nN,1)];
    Cee_dummy = [Cee zeros(o3,nN); zeros(nN,o3) zeros(nN,nN)];
    fwrite(fid_kal, [Xhat_t_t_dummy; Cee_dummy(:)], 'double');
    fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
    fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');
    
    if (flag_cov == 0)
        if (flag_ge == 1), rtplot_googleearth (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), date(1,:)), end;
        rtplot_matlab (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), 0, 0, 0, 0, flag_ms, ref_path, mat_path);
    else
        if (flag_ge == 1), rtplot_googleearth_cov (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date(1,:)), end;
        rtplot_matlab_cov (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), 0, 0, 0, 0, flag_ms, ref_path, mat_path);
    end
    if (flag_cpu == 0)
        rtplot_skyplot (1, azR, elR, conf_sat, pivot);
        rtplot_snr (snr_R(:,1));
    else
        rttext_sat (1, azR, elR, snr_R(:,1), conf_sat, pivot);
    end
    
    for t = 2 : length(time_GPS)
        
        if (mode_data == 0)
            Eph_t = rt_find_eph (Eph, time_GPS(t));
        else
            Eph_t = Eph(:,:,t);
        end
        
        kalman_goGPS_SA_cod_loop (time_GPS(t), Eph_t, iono_R, pr1_R(:,t), pr2_R(:,t), snr_R(:,t), 1);
        
        Xhat_t_t_dummy = [Xhat_t_t; zeros(nN,1)];
        Cee_dummy = [Cee zeros(o3,nN); zeros(nN,o3) zeros(nN,nN)];
        fwrite(fid_kal, [Xhat_t_t_dummy; Cee_dummy(:)], 'double');
        fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
        fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');
        
        if (flag_cov == 0)
            if (flag_ge == 1), rtplot_googleearth (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), date(t,:)), end;
            rtplot_matlab (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), check_on, check_off, check_pivot, check_cs, flag_ms, ref_path, mat_path);
        else
            if (flag_ge == 1), rtplot_googleearth_cov (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date(t,:)), end;
            rtplot_matlab_cov (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), check_on, check_off, check_pivot, check_cs, flag_ms, ref_path, mat_path);
        end
        if (flag_cpu == 0)
            rtplot_skyplot (t, azR, elR, conf_sat, pivot);
            rtplot_snr (snr_R(:,t));
        else
            rttext_sat (t, azR, elR, snr_R(:,t), conf_sat, pivot);
        end
        pause(0.01);
    end
    
    fclose(fid_kal);
    fclose(fid_sat);
    fclose(fid_conf);
    
%----------------------------------------------------------------------------------------------
% REAL-TIME: KALMAN FILTER ON PHASE AND CODE DOUBLE DIFFERENCES WITH/WITHOUT A CONSTRAINT
%----------------------------------------------------------------------------------------------
    
elseif (mode == 11)
    
    goGPS_realtime(filerootOUT, mode_vinc, flag_ms, flag_ge, flag_cov, flag_NTRIP, flag_ms_rtcm, ref_path, mat_path, pos_M, iono_R, pr2_M, pr2_R, ph2_M, ph2_R);
    
%----------------------------------------------------------------------------------------------
% REAL-TIME: UBLOX MONITORING
%----------------------------------------------------------------------------------------------
    
elseif (mode == 12)
    
    goGPS_ublox_monitor(filerootOUT);
    
    %interrupt execution
    break
    
%----------------------------------------------------------------------------------------------
% REAL-TIME: MASTER MONITORING
%----------------------------------------------------------------------------------------------
    
elseif (mode == 13)
    
    goGPS_master_monitor(filerootOUT, flag_NTRIP);
    
    %interrupt execution
    break
    
end

%----------------------------------------------------------------------------------------------
% INPUT/OUTPUT DATA FILE READING
%----------------------------------------------------------------------------------------------

%stream reading
% [time_GPS, time_R, time_M, pr1_R, pr1_M, ph1_R, ph1_M, snr_R, snr_M, ...
%  Eph, loss_R, loss_M, stream_R, stream_M] = load_stream (filerootIN);

%---------------------------------

%observation file (OBS) and ephemerides file (EPH) reading
%[time_GPS, time_R, time_M, pr1_R, pr1_M, ph1_R, ph1_M, snr_R, snr_M, ...
% pos_M, Eph, delay, loss_R, loss_M] = load_goGPSinput (filerootOUT);

%---------------------------------

%reading of the file with Kalman filter results
[Xhat_t_t, Yhat_t_t, Cee, azM, azR, elM, elR, distM, distR, ...
    conf_sat, conf_cs, pivot] = load_goGPSoutput (filerootOUT, mode_vinc);

%variable saving for final graphical representations
nObs = size(Xhat_t_t,2);
pos_KAL = zeros(3,nObs);
estim_amb = zeros(32,nObs);
sigma_amb = zeros(32,nObs);
for i = 1 : nObs
    if (mode_vinc == 0)
        pos_KAL(:,i) = [Xhat_t_t(1,i); Xhat_t_t(o1+1,i); Xhat_t_t(o2+1,i)];
        estim_amb(:,i) = Xhat_t_t(o3+1:o3+32,i);
        sigma_amb(:,i) = sqrt(diag(Cee(o3+1:o3+32,o3+1:o3+32,i)));
    else
        pos_KAL(:,i) = [Yhat_t_t(1,i); Yhat_t_t(2,i); Yhat_t_t(3,i)];
        estim_amb(:,i) = Xhat_t_t(o1+1:o1+32,i);
        sigma_amb(:,i) = sqrt(diag(Cee(o1+1:o1+32,o1+1:o1+32,i)));
    end
end

%----------------------------------------------------------------------------------------------

%computation and visualization time reading
[dt_acqR, dt_decR, dt_acqM, dt_decM, dt_saveI, dt_kal, dt_saveO, ...
    dt_plot, dt_ge, dt_sky, dt_snr] = load_goGPStime (filerootOUT);

%----------------------------------------------------------------------------------------------
% GEODETIC COORDINATES SAVING (TEXT FILE)
%----------------------------------------------------------------------------------------------

%cartesian coordinates (X,Y,Z)
X_KAL = pos_KAL(1,:)';
Y_KAL = pos_KAL(2,:)';
Z_KAL = pos_KAL(3,:)';

%coordinate transformation
[phi_KAL, lam_KAL, h_KAL] = cart2geod(X_KAL, Y_KAL, Z_KAL);
phi_KAL = phi_KAL * 180/pi;
lam_KAL = lam_KAL * 180/pi;

%file saving
fid_geod = fopen([filerootOUT '_geod.txt'], 'wt');
for i = 1 : length(phi_KAL)
    fprintf(fid_geod, '%.8f\t%.8f\t%.4f\n', phi_KAL(i), lam_KAL(i), h_KAL(i));
end
fclose(fid_geod);

%----------------------------------------------------------------------------------------------
% GOOGLE EARTH FILE SAVING (KML FILE)
%----------------------------------------------------------------------------------------------

%"clampedToGround" plots the points attached to the ground
%"absolute" uses the height defined in the tag <coordinates>;
%N.B. Google Earth uses orthometric heights
z_pos = 'clampedToGround';
%z_pos = 'absolute';
%URL to load the icon for the points
iconR = 'http://maps.google.com/mapfiles/kml/pal2/icon26.png';
point_colorR = 'FFF5005A';
%point size
scaleR = 0.4;

%file saving (Google Earth KML)
fid_kml = fopen([filerootOUT '.kml'], 'wt');
fprintf(fid_kml, '<?xml version="1.0" standalone="yes"?>\n');
fprintf(fid_kml, '<kml creator="goGPS" xmlns="http://earth.google.com/kml/2.2">\n');
fprintf(fid_kml, '  <Document>\n');
fprintf(fid_kml, '    <name><![CDATA[%s]]></name>\n', [filerootOUT '.kml']);
fprintf(fid_kml, '    <Snippet><![CDATA[created by goGPS]]></Snippet>\n');
for i = 1 : length(phi_KAL)
    fprintf(fid_kml, '      <Placemark>\n');
    fprintf(fid_kml, '        <Point>\n');
    fprintf(fid_kml, '          <altitudeMode>%s</altitudeMode>\n',z_pos);
    fprintf(fid_kml, '          <coordinates>%.6f,%.6f,%.6f</coordinates>\n',lam_KAL(i),phi_KAL(i),h_KAL(i));
    fprintf(fid_kml, '        </Point>\n');
    fprintf(fid_kml, '        <Snippet></Snippet>\n');
    fprintf(fid_kml, '        <Style>\n');
    fprintf(fid_kml, '          <IconStyle>\n');
    fprintf(fid_kml, '            <Icon>\n');
    fprintf(fid_kml, '              <href>%s</href>\n',iconR);
    fprintf(fid_kml, '            </Icon>\n');
    fprintf(fid_kml, '            <color>%s</color>\n',point_colorR);
    fprintf(fid_kml, '            <colorMode>normal</colorMode>\n');
    fprintf(fid_kml, '            <scale>%.2f</scale>\n',scaleR);
    fprintf(fid_kml, '          </IconStyle>\n');
    fprintf(fid_kml, '        </Style>\n');
    fprintf(fid_kml, '      </Placemark>\n');
end
fprintf(fid_kml, '  </Document>\n</kml>');
fclose(fid_kml);

%----------------------------------------------------------------------------------------------
% REPRESENTATION OF THE ESTIMATED TRAJECTORY (AND TEXT FILE SAVING)
%----------------------------------------------------------------------------------------------

%cartesian coordinates (X,Y,Z)
X_KAL = pos_KAL(1,:)';
Y_KAL = pos_KAL(2,:)';
Z_KAL = pos_KAL(3,:)';

%coordinate transformation
[EST_KAL, NORD_KAL, h_KAL] = cart2plan(X_KAL, Y_KAL, Z_KAL);

%trajectory plotting
figure
plot(EST_KAL, NORD_KAL, '.r');
xlabel('EST [m]'); ylabel('NORD [m]'); grid on;

%data saving
fid_plan = fopen([filerootOUT '_plan.txt'], 'wt');
for i = 1 : length(EST_KAL)
    fprintf(fid_plan, '%.8f\t%.8f\t%.4f\n', EST_KAL(i), NORD_KAL(i), h_KAL(i));
end
fclose(fid_plan);

%----------------------------------------------------------------------------------------------
% RINEX FILE SAVING
%----------------------------------------------------------------------------------------------

% ublox2RINEX(stream_R, [filerootOUT '_rover_RINEX.txt']);

%----------------------------------------------------------------------------------------------
% REPRESENTATION OF THE ESTIMATED ERROR COVARIANCE (AND TEXT FILE SAVING)
%----------------------------------------------------------------------------------------------

if (flag_cov == 1) & (mode_vinc == 0)
    
    %covariance propagation
    Cee_ENU = global2localCov(Cee([1 o1+1 o2+1],[1 o1+1 o2+1],:), Xhat_t_t([1 o1+1 o2+1],:));
    
    %trajectory plotting
    figure
    plot(EST_KAL, NORD_KAL, '.r'); axis equal
    xlabel('EST [m]'); ylabel('NORD [m]'); grid on;
    
    hold on
    for i = 1 : size(Cee_ENU,3)         % ellipse definition
        T = chol(Cee_ENU(1:2,1:2,i));   % Cholesky decomposition
        for j = 1 : size(x_circle,1)    % ellipse definition
            x_ellipse(j,:) = x_circle(j,:) * T + [EST_KAL(i), NORD_KAL(i)];
        end;
        plot(x_ellipse(:,1),x_ellipse(:,2));
    end;
    hold off
    
    %file saving
    fid_cov = fopen([filerootOUT '_cov.txt'], 'wt');
    for i = 1 : length(phi_KAL)
        fprintf(fid_cov, '%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n', Cee_ENU(1,1,i), Cee_ENU(1,2,i), ...
            Cee_ENU(1,3,i), Cee_ENU(2,2,i), Cee_ENU(2,3,i), Cee_ENU(3,3,i));
    end
    fclose(fid_cov);
    
end

%----------------------------------------------------------------------------------------------
% REPRESENTATION OF THE REFERENCE TRAJECTORY
%----------------------------------------------------------------------------------------------

%[EST_ref, NORD_ref, h_ref] = cart2plan(ref_path(:,1), ref_path(:,2),ref_path(:,3));

%%reference data plot
%plot(EST_ref, NORD_ref, 'm', 'LineWidth', 2)

%hold off

%----------------------------------------------------------------------------------------------
% REPRESENTATION OF THE 3D TRAJECTORY
%----------------------------------------------------------------------------------------------

%3D plot
figure
plot3(EST_KAL, NORD_KAL, h_KAL, '.r');
xlabel('EST [m]'); ylabel('NORD [m]'); zlabel('h [m]'); grid on

%----------------------------------------------------------------------------------------------
% REPRESENTATION OF THE VISIBLE SATELLITES CONFIGURATION
%----------------------------------------------------------------------------------------------

% if (mode == 1)
%
%    %figure
%    %imagesc(abs(conf_sat)), grid;
%    %colormap(1-gray);
%
%    figure
%    subplot('position',[0.1 0.35 0.8 0.55]);
%    hold on; grid on;
%    title('Satellite configuration')
%    for i = 1 : 32
%       index = find(abs(conf_sat(i,:)) == 1);
%       index_cs = intersect(index, find(conf_cs(i,:) == 1));
%       index_pivot = intersect(index, find(pivot == i));
%       if ~isempty(index)
%          plot(index,i*ones(size(index)),'b.-');
%          plot(index_pivot,i*ones(size(index_pivot)),'r.');
%          plot(index_cs,i*ones(size(index_cs)),'g.');
%       end
%    end
%    axis([1 size(conf_sat,2) 0.5 32.5])
%    hold off;
%    clear i index index_cs index_pivot
%
%    subplot('position',[0.1 0.1 0.8 0.2]);
%    hold on; grid on;
%    s1 = sum(abs(conf_sat));
%    plot(s1,'b.-');
%    s2 = [0; pivot(2:end) - pivot(1:end-1)];
%    plot(find(s2>0),s1(find(s2>0)),'r.')
%    s3 = sum(conf_cs);
%    plot(find(s3>0),s3(find(s3>0)),'g.');
%    axis([1 size(conf_sat,2) 0 max(s1)])
%    hold off;
%    clear s1 s2 s3
%
% end

%----------------------------------------------------------------------------------------------
% REPRESENTATION OF AZIMUTH, ELEVATION AND DISTANCE FOR VISIBILE SATELLITES
%----------------------------------------------------------------------------------------------

% if (mode == 1)
%
%    coltab = jet;
%
%    f1 = figure; hold on; grid on; title('Azimuth')
%    f2 = figure; hold on; grid on; title('Elevation')
%    f3 = figure; hold on; grid on; title('Distance')
%    k = 1;
%    for i = 1 : 32
%       index = find(abs(conf_sat(i,:)) == 1)';
%       if ~isempty(index)
%          %azimuth
%          figure(f1)
%          h = plot(index,azR(i,index),'b.-'); grid on;
%          set(h,'Color',coltab(2*i-1,:));
%          %elevation
%          figure(f2)
%          h = plot(index,elR(i,index),'r.-');
%          set(h,'Color',coltab(2*i-1,:));
%          %distance
%          figure(f3)
%          h = plot(index,distR(i,index)*1e-6,'g.-');
%          set(h,'Color',coltab(2*i-1,:));
%          %legend
%          list{k} = num2str(i);
%          k = k+1;
%       end
%    end
%    figure(f1); hold off; legend(list)
%    figure(f2); hold off; legend(list)
%    figure(f3); hold off; legend(list)
%    clear f1 f2 f3
%    clear i k h
%    clear coltab
%
% end

%----------------------------------------------------------------------------------------------
% REPRESENTATION OF THE S/N RATIO FOR MASTER AND ROVER
%----------------------------------------------------------------------------------------------

% if (mode == 1)
% 
%    coltab = jet;
%    coltab = [1 1 1; coltab([1 16 32 48 56],:)];
% 
%    figure
%    subplot(2,1,1);
%    imagesc(snr_R.*abs(conf_sat),[-0.5 60.5]);
%    title('Rover S/N ratio');
%    axis xy; colormap(coltab);
%    h = colorbar; set(h,'YTick',0:10:60);
%    subplot(2,1,2);
%    imagesc(snr_M.*abs(conf_sat),[-0.5 60.5]);
%    title('Master S/N ratio');
%    axis xy; colormap(coltab);
%    h = colorbar; set(h,'YTick',0:10:60);
% 
%    clear h coltab
% end

%----------------------------------------------------------------------------------------------
% REPRESENTATION OF THE COMBINATIONS OF ESTIMATED AMBIGUITIES
%----------------------------------------------------------------------------------------------

% if (mode == 1)
% 
%    for i = 1 : 32
%       index = find(conf_sat(i,:) == 1)';
%       index_cs = find(conf_cs(i,:) == 1)';
%       if ~isempty(index)
%          j = [1; find(index(2:end) - index(1:end-1) > 1)+1];
%          figure
%          %combination of estimated ambiguities
%          plot(index,estim_amb(i,index),'b.-'); grid on;
%          hold on
%          %cycle-slip
%          plot(index_cs,estim_amb(i,index_cs),'mo');
%          %combination of estimated ambiguities for new satellites
%          plot(index(j),estim_amb(i,index(j)),'g.');
%          %acceptance interval
%          plot(index, estim_amb(i,index) + sigma_amb(i,index),'r:');
%          plot(index, estim_amb(i,index) - sigma_amb(i,index),'r:');
%          hold off
%          title(['Combination of estimated ambiguities between PIVOT and SATELLITE ',num2str(i)]);
%       end
%    end
% 
% end

%----------------------------------------------------------------------------------------------
% STATISTICS COMPUTATION AND VISUALIZATION
%----------------------------------------------------------------------------------------------

if (mode < 10) & (mode_vinc == 0) & (~isempty(ref_path))
    %coordinate transformation
    [EST_REF, NORD_REF, h_REF] = cart2plan(ref_path(:,1), ref_path(:,2), ref_path(:,3));
    
    ref = [EST_REF NORD_REF h_REF];
    
    [dist2D, proj] = ref_2d_projection(ref,EST_KAL,NORD_KAL); %#ok<NASGU>
    
    fprintf('\n');
    fprintf('-------- STATISTICS ----------');
    fprintf('\n');
    fprintf('Mean2D: %7.4f m\n',mean(dist2D));
    fprintf('Std2D:  %7.4f m\n',std(dist2D));
    fprintf('RMS2D:  %7.4f m\n\n',sqrt(std(dist2D)^2+mean(dist2D)^2));
    
    [dist3D,proj] = ref_3d_projection(ref,EST_KAL,NORD_KAL,h_KAL);
    
    fprintf('Mean3D: %7.4f m\n',mean(dist3D));
    fprintf('Std3D:  %7.4f m\n',std(dist3D));
    fprintf('RMS3D:  %7.4f m\n',sqrt(std(dist3D)^2+mean(dist3D)^2));
    fprintf('--------------------------------\n\n');
end

%----------------------------------------------------------------------------------------------

%re-enable MATLAB warnings
warning on

%evaluate computation time
toc

%__________________________________/\\\\\\\\\\\\___/\\\\\\\\\\\\\______/\\\\\\\\\\\___
% ________________________________/\\\//////////___\/\\\/////////\\\__/\\\/////////\\\_
%  ___/\\\\\\\\___________________/\\\______________\/\\\_______\/\\\_\//\\\______\///__
%   __/\\\////\\\______/\\\\\_____\/\\\____/\\\\\\\__\/\\\\\\\\\\\\\/___\////\\\_________
%    _\//\\\\\\\\\____/\\\///\\\___\/\\\___\/////\\\__\/\\\/////////________\////\\\______
%     __\///////\\\___/\\\__\//\\\__\/\\\_______\/\\\__\/\\\____________________\////\\\___
%      __/\\_____\\\__\//\\\__/\\\___\/\\\_______\/\\\__\/\\\_____________/\\\______\//\\\__
%       _\//\\\\\\\\____\///\\\\\/____\//\\\\\\\\\\\\/___\/\\\____________\///\\\\\\\\\\\/___
%        __\////////_______\/////_______\////////////_____\///_______________\///////////_____
%
%(ASCII-art from http://patorjk.com/software/taag/)