% =========================================================================
%   OBJECT goKalmanFilter
% =========================================================================
%
% DESCRIPTION:
%   Object to manage a kalman filter solution
%
% EXAMPLE:
%   goKF = goKalmanFilter();
%
% REQUIRES:
%   goObservation:    goGPS object to store the input observations
%    - iniReader:     class to pass the parameters of the receivers
%      - cprintf:     http://www.mathworks.com/matlabcentral/fileexchange/24093-cprintf-display-formatted-colored-text-in-the-command-window
%
% LIST of METHODS
%
%  setDefaultVariances(obj)
%  setCurrentParameters(obj)
%  allocateMemory(obj, nRec, nFreq)
%  init(obj, goObs)
%  init_T(obj, mode)
%  init_Xhat_t_t(obj, goObs, mode)
%  init_X_t1_t(obj)
%  init_Cee(obj, nRec, mode)
%  init_doppler(obj, goObs)
%  init_KxDOP(obj, mode)
%
%----------------------------------------------------------------------------------------------
% Copyright (C) 2009-2013 Mirko Reguzzoni, Eugenio Realini
%----------------------------------------------------------------------------------------------
%
%    Code contributed by Andrea Gatti and Lisa Pertusini
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

classdef goKalmanFilter < handle
    
    properties (GetAccess = 'public', SetAccess = 'public')
        
        % Initialization maximum number of satellites
        nSat = 32;
        
        % Initialization variances
        sigmaq0;
        
        % User defined variances
        sigmaq;
        
        %minimum number of satellites to be used in the Kalman filter
        min_nsat = 2;
        
        %cut-off [degrees]
        cutoff = 0;
        
        %signal-to-noise ratio threshold [dB]
        snr_threshold = 0;
        
        %cycle slip threshold [cycles]
        cs_threshold = 10;
        
        %parameter used to select the weight mode for GPS observations
        %          - weights=0: same weight for all the observations
        %          - weights=1: weight based on satellite elevation
        %          - weights=2: weight based on signal-to-noise ratio
        %          - weights=3: weight based on combined elevation and signal-to-noise ratio
        weights = 2;
        
        %weight function parameters
        snr_a = 30;
        snr_0 = 10;
        snr_1 = 50;
        snr_A = 30;
        
        % KF mode
        % possible values:
        % 1) static
        % 2) constant velocity
        % 3) constant acceleration
        % 4) constant velocity + attitude estimation
        % 5) constant acceleration + attitude estimation
        mode = 1;
        
        %order of the dynamic model polynomial
        order = 2;
        
        % number of parameters to be estimated
        nPar = 3;
        
        %ambiguity restart method
        amb_restart_method = 1;
        
        %transition matrix
        T = [];
        
        %identity matrix
        I = [];
        
        %tmp state estimation at time t for each receiver
        Xhat_t_t_R = {};
        
        %state estimation at time t
        Xhat_t_t = [];
        
        %state estimation at time t+1 (using dynamics only)
        X_t1_t = [];
        
        %estimation error covariance matrix at time t
        Cee = [];
        
        %satellite configuration (1: visible, 0: not visible) at time t
        conf_sat = [];
        
        %cycle-slip configuration (1: cs, 0: no cs) at time t
        conf_cs = [];
        
        %index of the current pivot satellite for each receiver
        pivot = [];
        
        %index of the previous pivot satellite for each receiver
        pivot_old = [];
        
        %number of unknown phase ambiguities
        nN = [];
        
        %method used to estimate phase ambiguities
        %          - amb_estim_method=0: observed code - phase difference
        %          - amb_estim_method=1: Kalman-predicted code - phase difference
        %          - amb_estim_method=2: Least squares adjustment
        amb_estim_method = 0;
        
        %estimation sampling rate
        interval = 1; %default 1 Hz (to avoid problems with real-time modes)
        
        %threshold on the condition number on the eigenvalues of the N matrix (least squares)
        cond_num_threshold = 1e6;
        
        %azimuth, elevation and distance of satellites with respect to the
        %ROVER and MASTER
        satCoordR = struct('az',zeros(32,1),'el',zeros(32,1),'dist',zeros(32,1)); % for each receiver: azimuth (az), elevation (el), distance (dist)
        satCoordM = struct('az',zeros(32,1),'el',zeros(32,1),'dist',zeros(32,1)); % azimuth (az), elevation (el), distance (dist)
%         satCoordR; % for each receiver: azimuth (az), elevation (el), distance (dist)
%         satCoordM; % azimuth (az), elevation (el), distance (dist)
        
        % DILUTION OF PRECISION
        xDOP;  % P, H, V, KP, KH, KV
        
        %Doppler-predicted range (ROVER and MASTER)
        doppler_pred_range_R;
        doppler_pred_range_M;
        
        %flag containing initialization status
        initKF = false;
    end
    
    methods
        % Creator (Brahma)
        function obj = goKalmanFilter(goObs, mode, sampling_rate)
            if nargin < 3
                sampling_rate = 1;
            end
            if nargin < 2
                mode = 1;
            end
            if mode > 5
                mode = 1;
            end
            switch mode
                %if there is only one receiver, do not estimate the attitude
                case 1, obj.nPar = 3;   % static filter (3 positions)
                case 2, obj.nPar = 6;   % const.velocity filter (3 positions+3 velocities)
                case 3, obj.nPar = 9;   % const.acceleration filter (3 positions+3 velocities+3 accelerations)
                    %more then one receiver
                case 4, obj.nPar = 12;  % const.velocity filter + attitude angles and variations
                case 5, obj.nPar = 15;  % const.acceleration filter + attitude angles and variations
            end
            obj.mode = mode;
            obj.interval = 1/sampling_rate;	% Init estimation sampling rate
            obj.setDefaultVariances();      % Init variances
            obj.setCurrentParameters();     % Init current parameters
            obj.allocateMemory(goObs.getNumRec(), goObs.getGNSSnFreq(goObs.idGPS)); % only GPS observations
            
            obj.init(goObs)
        end
        
        % Destructor (Shiva)
        function delete(obj)
        end
        
        % return the status of the KF initialization
        function isI = isInitialized(obj)
            isI = obj.initKF;
        end
    end
    
    % Initialization functions
    
    methods (Access = 'private')
        % Function to initialize all the variances used in the KF
        function setDefaultVariances(obj)
            %variance of initial state %big variances also for velocities
            %because the receivers are assumed to be still at the first
            %epoch, but then they move
            obj.sigmaq0.pos = 9;        %[m^2]
            obj.sigmaq0.vel = 9;        %[m^2/s^2]
            obj.sigmaq0.acc = 9;        %[m^4/s^4]
            obj.sigmaq0.ang = 9;        %[rad^2]
            obj.sigmaq0.ang_vel = 9;    %[rad^2/s^2]
            
            %variance of ambiguity combinations [cycles]
            obj.sigmaq0.N = 1000;
            
            %variance of velocity coordinates [m^2/s^2]
            obj.sigmaq.vE = 1e-1;
            obj.sigmaq.vN = 1e-1;
            obj.sigmaq.vU = 1e-1;
            obj.sigmaq.vel = 1e-0;
            
            %variance of code observations [m^2]
            % sigmaq.cod1 = 0.36;
            obj.sigmaq.cod1 = 9;
            obj.sigmaq.cod2 = 0.16;
            
            %variance of phase observations [m^2]
            %(maximize to obtain a code-based solution)
            % sigmaq.ph = 0.000004;
            obj.sigmaq.ph = 0.001;
            % sigmaq.ph = 0.001e30;
            
            %variance of DEM height [m^2]
            %(maximize to disable DEM usage)
            % sigmaq_dtm = 0.09;
            obj.sigmaq.dtm = 1e30;
        end
        
        % Function to get the current global parameters
        function setCurrentParameters(obj)
            global sigmaq0 sigmaq0_N
            global cutoff snr_threshold cond_num_threshold
            obj.sigmaq0.pos = sigmaq0;
            obj.sigmaq0.N = sigmaq0_N;
            obj.cutoff = cutoff;
            obj.snr_threshold = snr_threshold;
            obj.cond_num_threshold = cond_num_threshold;
        end
        
        % Function to preallocate memory for the KF matrices
        function allocateMemory(obj, nRec, nFreq)
            nPar = obj.nPar;
            nSat = obj.nSat;
            obj.nN = (nSat*nFreq*nRec);                 %number of unknown phase ambiguities
            obj.T = eye(nPar+obj.nN);                   %transition matrix
            obj.I = eye(nPar+obj.nN);                   %identity matrix
            obj.Xhat_t_t_R = cell(nRec,1);              %tmp state estimation at time t for each receiver
            for r=1:(nRec),
                obj.Xhat_t_t_R{r} = zeros(nPar+nSat*nFreq,1);
            end
            obj.Xhat_t_t = zeros(nPar+obj.nN,1);        %state estimation at time t
            obj.X_t1_t = zeros(nPar+obj.nN,1);          %state estimation at time t+1 (using dynamics only)
            obj.Cee = zeros(nPar+obj.nN);               %estimation error covariance matrix at time t
            obj.conf_sat = zeros(nSat, nRec);           %satellite configuration (1: visible, 0: not visible) at time t
            obj.conf_cs = zeros(nSat, nRec);            %cycle-slip configuration (1: cs, 0: no cs) at time t
            obj.pivot = zeros(1,nRec);                  %index of the current pivot satellite
            obj.pivot_old = zeros(1,nRec);              %index of the previous pivot satellite
            
            % spherical coordinates initialization
            
            % pre-allocate the structure for the vectors obj.satCoordR(1|2|3).az,
            % obj.satCoordR(1|2|3).el obj.satCoordR(1|2|3).dist
            obj.satCoordR(nRec) = struct('az',zeros(nSat,1),'el',zeros(nSat,1),'dist',zeros(nSat,1));
            for r=1:(nRec-1),
                obj.satCoordR(r) = struct('az',zeros(nSat,1),'el',zeros(nSat,1),'dist',zeros(nSat,1));
            end
            obj.satCoordM = struct('az',zeros(nSat,1),'el',zeros(nSat,1),'dist',zeros(nSat,1));
            
            % DILUTION OF PRECISION
            obj.xDOP = struct('P', zeros(nRec,1), 'H', zeros(nRec,1), 'V', zeros(nRec,1), 'KP', zeros(nRec,1), 'KH', zeros(nRec,1), 'KV', zeros(nRec,1));
            
            % Doppler-predicted range
            obj.doppler_pred_range_R = zeros(nSat,nRec,nFreq);
            obj.doppler_pred_range_M = zeros(nSat,1,nFreq);
        end
        
        % Function to fill KF initial matrices
        function init(obj, goObs)
            obj.init_T(obj.mode);
            obj.init_Xhat_t_t(goObs);
            obj.init_X_t1_t();
            obj.init_Cee(goObs.getNumRec());
            obj.init_doppler(goObs);
            obj.init_KxDOP();
            % set the status of the KF initialization
            obj.initKF = true;
        end
        
        % Function to fill Transition matrix according to the number of parameters
        function init_T(obj, mode)
            nPar = obj.nPar;
            % Transition matrix filling
            switch(mode)
                case 1,
                    % obj.T(1:nPar,1:nPar) = eye(nPar);
                    % redundant, becuse T is already an identity matrix
                case 2,
                    % obj.T(1:nPar,1:nPar) = eye(nPar);
                    obj.T(1:2:nPar,2:2:nPar) = diag(ones(3,1)*obj.interval);
                case 3,
                    % obj.T(1:nPar,1:nPar) = eye(nPar);
                    obj.T(1:3:nPar,2:3:nPar) = diag(ones(3,1)*obj.interval);
                    obj.T(2:3:nPar,3:3:nPar) = diag(ones(3,1)*obj.interval);
                case 4,
                    % obj.T(1:nPar,1:nPar) = eye(nPar);
                    obj.T(1:2:(nPar-6),2:2:(nPar-6)) = diag(ones(3,1)*obj.interval);
                case 5,
                    % obj.T(1:nPar,1:nPar) = eye(nPar);
                    obj.T(1:3:(nPar-6),2:3:(nPar-6)) = diag(ones(3,1)*obj.interval);
                    obj.T(2:3:(nPar-6),3:3:(nPar-6)) = diag(ones(3,1)*obj.interval);
            end
            % note that: the remaining part of the T matrix has already
            % been created as an identity matrix sized nN.
        end
        
        % initialization of the parameter vector for all receivers
        function init_Xhat_t_t(obj, goObs, mode)   %% to initialize Xhat_t_t_R: cell of [nPar+nSat*nFreq,1] and Xhat_t_t;
            nSat = obj.nSat;
            nRec = goObs.getNumRec();
            nPar = obj.nPar;
            nN = obj.nN;
            nFreq = goObs.getGNSSnFreq(goObs.idGPS);
            
            % define logical matrices for the satellites in view
            commonSat_pr = false(nSat,nRec);
            commonSat_pr_ph = false(nSat,nRec);
            
            % select only the satellites in common between master and rover
            pr_R = reshape(goObs.getGNSSpr_R(goObs.idGPS, 0, 0, 1, 0),nSat,nRec,nFreq);
            ph_R = reshape(goObs.getGNSSph_R(goObs.idGPS, 0, 0, 1, 0),nSat,nRec,nFreq);
            pr_M = goObs.getGNSSpr_M(goObs.idGPS, 0, 1, 0);
            ph_M = goObs.getGNSSph_M(goObs.idGPS, 0, 1, 0);
            
            if (nFreq == 2) % double frequency
                % reshape the pr and ph observations to have nSat rows,
                % nRec columns, nFreq planes
                commonSat_pr = (pr_R(:,:,1) ~= 0) & repmat((pr_M(:,:,1) ~= 0),1,nRec) & (pr_R(:,:,2) ~= 0) & repmat((pr_M(:,:,2) ~= 0),1,nRec);
                commonSat_pr_ph = (pr_R(:,:,1) ~= 0) & repmat((pr_M(:,:,1) ~= 0),1,nRec) & (ph_R(:,:,1) ~= 0) & repmat((ph_M(:,:,1) ~= 0),1,nRec) & ...
                    (pr_R(:,:,2) ~= 0) & repmat((pr_M(:,:,2) ~= 0),1,nRec) & (ph_R(:,:,2) ~= 0) & repmat((ph_M(:,:,2) ~= 0),1,nRec);
            else
                if (nFreq == 1) % single frequency
                    commonSat_pr = (pr_R(:,:,1) ~= 0) & repmat((pr_M(:,:,1) ~= 0),1,nRec);
                    commonSat_pr_ph = (pr_R(:,:,1) ~= 0) & repmat((pr_M(:,:,1) ~= 0),1,nRec) & ...
                        (ph_R(:,:,1) ~= 0) & repmat((ph_M(:,:,1) ~= 0),1,nRec);
                else
                    % to be used for nFreq>2
                end
            end
            
            
            %------------------------------------------------------------------------------------
            % APPROXIMATE POSITION for the rovers and the master
            %-----------------------------------------------------------------------------------
            
            [XR0 flag_XR] = goObs.getX0_R(0);  %[sized 3xnRec]
            [XM0 flag_M] = goObs.getX0_M();    %[sized 3x1]
            
            %--------------------------------------------------------------------------------------------
            % KALMAN FILTER INITIAL STATE
            %--------------------------------------------------------------------------------------------
            
            %%zeroes vectors useful in the matrices definition
            %Z_om_1 = zeros(o1-1,1);
            
            %variances of the phase ambiguity estimate
            sigma2_N = zeros(nN,1);
            
            %define vectors for the SP3 ephemerides
            SP3_time = goObs.getGNSS_SP3time();
            SP3_coor = goObs.getGNSS_SP3coordinates();
            SP3_clck = goObs.getGNSS_SP3clock();
            
            %logical indexes of the remaining satellites after cutoff
            goodSat_pr_M = false(nSat,1);
            goodSat_pr_R = false(nSat,nRec);
            
            % initialize the covariance matrix for the rover positions
            cov_XR = zeros(3,3,nRec);
            
            % select only the satellites in view for the master
            % (NOT a logical vector!)
            % _init because it is before applying the cutoff
            if (nFreq == 2) % double frequency
                sat_pr_M_init = find((pr_M(:,:,1) ~= 0) & (pr_M(:,:,2) ~= 0));
            else
                if (nFreq == 1)
                    sat_pr_M_init = find(pr_M(:,:,1) ~= 0);
                else
                    % to be used for nFreq>2
                end
            end
            
            % Compute initial receiver and satellite position and clocks
            % using the cutoff value for the master
            [XM, dtM, XS, dtS, XS_tx, VS_tx, time_tx, ...
                err_tropo_M, err_iono_M, sat_pr_M, ...
                elM, azM, distM, ...
                cov_XM, var_dtM] ...
                = init_positioning(goObs.getTime_Ref(), pr_M(sat_pr_M_init,1,1), goObs.getGNSSsnr_M(goObs.idGPS, sat_pr_M_init,1,1), ...
                goObs.getGNSSeph(goObs.idGPS), SP3_time, SP3_coor, SP3_clck, ...
                goObs.getIono(), XM0, [], [], sat_pr_M_init, obj.cutoff, obj.snr_threshold, flag_M, 0);
            obj.satCoordM(sat_pr_M).el = elM;
            obj.satCoordM(sat_pr_M).az = azM;
            obj.satCoordM(sat_pr_M).dist = distM;
            clear elM azM distM
            goodSat_pr_M(sat_pr_M) = 1;
            
            %having at least 4 satellites in common in view
            %if (sum(commonSat_pr) >= 4)
            % initialization of variables outside the loop
            sat_pr_R_init = zeros(obj.nSat, obj.nRec);
            err_tropo_R = zeros(obj.nSat, obj.nRec);
            err_iono_R = zeros(obj.nSat, obj.nRec);
            var_dtR = zeros(nRec,1);
            dtR = zeros(nRec,1);
            cond_num = zeros(nRec,1);
            for r=1:nRec
                %having at least 4 satellites in view from the master
                %station after applying the cutoff
                if (sum(goodSat_pr_M) < 4); return; end
                % select only the satellites in view for the rover
                % after the Master cutoff
                % _init because it is before applying the rover cutoff
                % (it is a logical vector)
                sat_pr_R_init(:,r) = (goodSat_pr_M ~= 0) & (commonSat_pr(:,r) ~= 0);
                
                [obj.XR(:,r), dtR(r), ~, ~, ~, ~, ~, ...
                    err_tropo_R(:,r), err_iono_R(:,r), sat_pr_R, ...
                    obj.satCoordR(sat_pr_R,r).el, obj.satCoordR(sat_pr_R,r).az, obj.satCoordR(sat_pr_R,r).dist, ...
                    cov_XR(:,:,r), var_dtR(r), obj.xDOP(r).P, obj.xDOP(r).H, obj.xDOP(r).V, cond_num(r)] ...
                    = init_positioning(goObs.getTime_Ref(), pr_R(sat_pr_R_init(:,r),r,1), goObs.getGNSSsnr_R(goObs.idGPS, sat_pr_R_init(:,r),r,1), ...
                    goObs.getGNSSeph(goObs.idGPS), SP3_time, SP3_coor, SP3_clck, ...
                    goObs.getIono(), XR0(:,r), XS, dtS, sat_pr_R_init(:,r), obj.cutoff, obj.snr_threshold, flag_XR(r), 1);
                
                goodSat_pr_R(sat_pr_R,r) = 1;
            end
            %keep only satellites that rover and master have in common
            %over the cutoff threshold
            obj.goodSat_pr = goodSat_pr_R;
            
            % atmospheric errors of the receivers for the satellites in
            % common between the master and each rover (nRec columns)
            err_tropo_R = err_tropo_R*obj.goodSat_pr;
            err_iono_R  = err_iono_R*obj.goodSat_pr;
            err_tropo_M = repmat(err_tropo_M,1,nRec)*obj.goodSat_pr;
            err_iono_M  = repmat(err_iono_M,1,nRec)*obj.goodSat_pr;
            
            %apply cutoffs also to phase satellites
            obj.goodSat_pr_ph = commonSat_pr_ph & obj.goodSat_pr;
            
            % fill doppler variables
            %for i = 1:sum(obj.goodSat_pr,r)
            if (~isempty(goObs.getClockDrift_M()) && goObs.getGNSSdop_M(goObs.idGPS, goodSat_pr_M, 1, 1, 1) == 0 && any(goObs.getGNSSeph(goObs.idGPS)))
                %satObs = goObs.getSatObservation(goObs.idGPS, goodSat_pr_M);
                [goObs.getGNSSdop_M(goObs.idGPS, goodSat_pr_M, 1, 1, 2), goObs.getGNSSdop_M(goObs.idGPS, goodSat_pr_M, 1, 1, 1)] ...
                    = doppler_shift_approx(goObs.getPos_M(), zeros(3,1), ...
                    XS_tx', VS_tx', time_tx, ... %satObs.X(goodSat_pr_M,:)', satObs.V(goodSat_pr_M,:)', satObs.time(goodSat_pr_M,1), ...
                    goObs.getClockDrift_M(), goodSat_pr_M, goObs.getGNSSeph(goObs.idGPS));
            end
            
            %--------------------------------------------------------------------------------------------
            % SATELLITE CONFIGURATION SAVING AND PIVOT SELECTION
            %--------------------------------------------------------------------------------------------
            
            %satellites configuration: code only (-1), both code and phase (+1);
            conf_sat = zeros(nSat,nRec);
            conf_sat(obj.goodSat_pr) = -1;
            conf_sat(obj.goodSat_pr_ph) = +1;
            
            %cycle-slip configuration (no cycle-slip)
            conf_cs = zeros(nSat,nRec);
            
            % initialize pivot
            obj.pivot = zeros(1,nRec);
            pivot_index = zeros(1,nRec);
            %previous pivot
            obj.pivot_old = zeros(1,nRec);
            
            %current pivot
            for r=1:nRec
                % if there is at least one visible satellite for each
                % receiver (with both code and phase obervations)
                if ~isempty(obj.goodSat_pr_ph(:,r))
                    % find the index for the most elevated satellite (PIVOT)
                    [null_max_elR, pivot_index(r)] = max(obj.satCoordR(obj.goodSat_pr_ph(:,r),r).el);
                    s_id = 1:nSat;                          % all the satellites
                    s_id = s_id(obj.goodSat_pr_ph(:,r));        % extract available satellites
                    obj.pivot(r) = s_id(pivot_index(r));       % get the pivot satellite
                else %only code observations
                    [null_max_elR, pivot_index(r)] = max(obj.satCoordR(obj.goodSat_pr(:,r),r).el);
                    s_id = 1:obj.nSat;                      % all the satellites
                    s_id = s_id(obj.goodSat_pr(:,r));           % extract available satellites
                    obj.pivot(r) = s_id(pivot_index(r));       % get the pivot satellite
                end
                
                %if at least 4 satellites are available after the cutoffs, and if the
                % condition number in the least squares does not exceed the threshold
                if (sum(obj.goodSat_pr,r) >= 4 && cond_num(r) < obj.cond_num_threshold)
                    
                    if isempty(cov_XR(:,:,r)) %if it was not possible to compute the covariance matrix
                        cov_XR(:,:,r) = obj.sigmaq0 * eye(3);
                    end
                    obj.sigma2_XR_R(:,r) = diag(cov_XR(:,:,r));
                else
                    return
                end
            end
            
            %do not use least squares ambiguity estimation
            % NOTE: LS amb. estimation is automatically switched off if the number of
            % satellites with phase available is not sufficient
            
            %ambiguity initialization: initialized value if the satellite is visible,
            %0 if the satellite is not visible
            N1 = zeros(nSat,nRec);
            N2 = zeros(nSat,nRec);
            sigma2_N1 = zeros(nSat,nRec);
            sigma2_N2 = zeros(nSat,nRec);
            cov_N1 = zeros(3,3,nRec);
            cov_N2 = zeros(3,3,nRec);
            N = zeros(nSat, nRec);
            for r=1:nRec
                if (sum(obj.goodSat_pr(:,r)) + sum(obj.goodSat_pr_ph(:,r)) - 2 <= 3 + sum(obj.goodSat_pr_ph(:,r)) - 1 || sum(obj.goodSat_pr_ph(:,r)) <= 4)
                    
                    %computation of the phase double differences in order to estimate N
                    if ~isempty(obj.goodSat_pr_ph(:,r))
                        [N1(obj.goodSat_pr_ph(:,r),r), sigma2_N1(obj.goodSat_pr_ph(:,r),r)] = amb_estimate_observ(pr_R(obj.goodSat_pr_ph(:,r),r,1), pr_M(obj.goodSat_pr_ph(:,r),1,1), ph_R(obj.goodSat_pr_ph(:,r),r,1), ph_M(obj.goodSat_pr_ph(:,r),1,1), obj.pivot(r), obj.goodSat_pr_ph(:,r), 1);
                        [N2(obj.goodSat_pr_ph(:,r),r), sigma2_N2(obj.goodSat_pr_ph(:,r),r)] = amb_estimate_observ(pr_R(obj.goodSat_pr_ph(:,r),r,2), pr_M(obj.goodSat_pr_ph(:,r),1,2), ph_R(obj.goodSat_pr_ph(:,r),r,2), ph_M(obj.goodSat_pr_ph(:,r),1,2), obj.pivot(r), obj.goodSat_pr_ph(:,r), 2);
                    end
                    
                    if (nFreq == 2)
                        N(:,r) = [N1(:,r); N2(:,r)];
                        sigma2_N(:,r) = [sigma2_N1(:,r); sigma2_N2(:,r)];
                    else
                        if (nFreq == 1)
                            N(:,r) = N1(:,r);
                            sigma2_N(:,r) = sigma2_N1(:,r);
                        else
                            % to be used for nFreq > 2
                        end
                    end
                    
                    %use least squares ambiguity estimation
                else
                    
                    %ROVER positioning improvement with code and phase double differences
                    if ~isempty(obj.goodSat_pr_ph(:,r))
                        [obj.XR(:,r), N1(obj.goodSat_pr_ph(:,r),r), cov_XR(:,:,r), cov_N1(:,:,r), obj.xDOP(r).P, obj.xDOP(r).H, obj.xDOP(r).V] ...
                            = LS_DD_code_phase(obj.XR(:,r), XM, XS(obj.goodSat_pr_ph,:), ...
                            pr_R(obj.goodSat_pr_ph(:,r),r,1), ph_R(obj.goodSat_pr_ph(:,r),r,1), snr_R(obj.goodSat_pr_ph(:,r),r,1), ...
                            pr_M(obj.goodSat_pr_ph(:,r),1,1), ph_M(obj.goodSat_pr_ph(:,r),1,1), snr_M(obj.goodSat_pr_ph(:,r),1,1), ...
                            obj.satCoordR(obj.goodSat_pr_ph(:,r),r).el, obj.satCoordM(obj.goodSat_pr_ph(:,r),1).el, ...
                            err_tropo_R(obj.goodSat_pr_ph(:,r),r), err_iono_R(obj.goodSat_pr_ph(:,r),r), ...
                            err_tropo_M(obj.goodSat_pr_ph(:,r),r), err_iono_M(obj.goodSat_pr_ph(:,r),r), ...
                            pivot_index(r), nFreq);
                        [null_XR, N2(obj.goodSat_pr_ph(:,r),r), null_cov_XR, cov_N2(:,:,r)] ...
                            = LS_DD_code_phase(obj.XR(:,r), XM, XS(obj.goodSat_pr_ph,:), ...
                            pr_R(obj.goodSat_pr_ph(:,r),r,2), ph_R(obj.goodSat_pr_ph(:,r),r,2), snr_R(obj.goodSat_pr_ph(:,r),r,2), ...
                            pr_M(obj.goodSat_pr_ph(:,r),1,2), ph_M(obj.goodSat_pr_ph(:,r),1,2), snr_M(obj.goodSat_pr_ph(:,r),1,2), ...
                            obj.satCoordR(obj.goodSat_pr_ph(:,r),r).el, obj.satCoordM(obj.goodSat_pr_ph(:,r),1).el, ...
                            err_tropo_R(obj.goodSat_pr_ph(:,r),r), err_iono_R(obj.goodSat_pr_ph(:,r),r), ...
                            err_tropo_M(obj.goodSat_pr_ph(:,r),r), err_iono_M(obj.goodSat_pr_ph(:,r),r), ...
                            pivot_index(r), nFreq);
                    end
                    
                    if isempty(cov_XR(:,:,r)) %if it was not possible to compute the covariance matrix
                        cov_XR(:,:,r) = obj.sigmaq0 * eye(3);
                    end
                    obj.sigma2_XR_R(:,r) = diag(cov_XR(:,:,r));
                    
                    if isempty(cov_N1(:,:,r)) %if it was not possible to compute the covariance matrix
                        cov_N1(:,:,r) = obj.sigmaq0.N * eye(length(obj.goodSat_pr_ph(:,r),r));
                    end
                    
                    if isempty(cov_N2(:,:,r)) %if it was not possible to compute the covariance matrix
                        cov_N2(:,:,r) = obj.sigmaq0.N * eye(length(obj.goodSat_pr_ph(:,r),r));
                    end
                    
                    if (nFreq == 2)
                        N(:,r) = [N1(:,r); N2(:,r)];
                        obj.sigma2_N(obj.goodSat_pr_ph(:,r),r) = diag(cov_N1(:,:,r));
                        obj.sigma2_N(nSat+obj.goodSat_pr_ph(:,r),r) = diag(cov_N2(:,:,r));
                    else
                        if (nFreq == 1)
                            N(:,r) = N1(:,r);
                            obj.sigma2_N(obj.goodSat_pr_ph(:,r),r) = diag(cov_N1(:,:,r));
                        else
                            % to be used for nFreq > 2
                        end
                    end
                end
                %initialization of the initial point with 3/6/9(positions/velocities/accelerations) +
                % 3 attitude angles + 3 attitude variations +
                %32 or 64 (N combinations) variables
                switch(mode)
                    case 1,
                        obj.Xhat_t_t_R{r} = [obj.XR(1,r); obj.XR(2,r); obj.XR(3,r); N(:,r)];
                    case {2,4},
                        obj.Xhat_t_t_R{r} = [obj.XR(1,r); 0; obj.XR(2,r); 0; obj.XR(3,r); 0; N(:,r)];
                    case {3,5},
                        obj.Xhat_t_t_R{r} = [obj.XR(1,r); 0; 0; obj.XR(2,r); 0; 0; obj.XR(3,r); 0; 0; N(:,r)];
                end
            end
            %in the state vector the coordinates of the baricenter are
            %considered
            
            % obj.XR contains in the first column the baricenter coordinates,
            % then the coordinates of each receiver
            obj.XR = [mean(obj.XR,2) obj.XR];
            
            %if there is only one receiver, delete the baricenter column
            %(useless)
            if nRec == 1;
                obj.XR(:,1) = [];
            end
            
            %build the state vector
            switch(mode)
                case 1,
                    obj.Xhat_t_t(1:3) = [obj.XR(1,1); obj.XR(2,1); obj.XR(3,1)];
                case {2,4},
                    obj.Xhat_t_t(1:6) = [obj.XR(1,1); 0; obj.XR(2,1); 0; obj.XR(3,1); 0];
                case {3,5},
                    obj.Xhat_t_t(1:9) = [obj.XR(1,1); 0; 0; obj.XR(2,1); 0; 0; obj.XR(3,1); 0; 0];
            end
            
            switch(mode)
                case {1,2,3},   %when not estimating the attitude
                    obj.Xhat_t_t (nPar+1:end) = N(:);
                case {4,5},     % when estimating the attitude (roll, pitch, yaw angles)
                    obj.attitude = goObs.getInitialAttitude();
                    obj.Xhat_t_t (nPar-6+1:nPar) = [obj.attitude.roll; 0; obj.attitude.pitch; 0; obj.attitude.yaw; 0];
                    obj.Xhat_t_t (nPar+1:end) = N(:);
            end
        end
        
        % initialization of point estimation at step t+1 ==
        % estimation at step t, because the initial velocity is equal to 0
        function init_X_t1_t(obj)
            X_t1_t = obj.T*obj.Xhat_t_t;
        end
        
        % initialization of state covariance matrix
        function init_Cee(obj, nRec, mode)
            obj.sigma2_XR = (1/nRec)^2*obj.sigma2_XR_R;     % variance propagation
            %positions
            obj.Cee(1,1) = obj.sigma2_XR(1);
            switch(mode)
                case {1,2,3}
                    o1 = obj.nPar/3;
                    obj.Cee(1+o1,1+o1) = obj.sigma2_XR(2);
                    obj.Cee(1+o1*2,1+o1*2) = obj.sigma2_XR(3);
                case {4,5}
                    o1 = (obj.nPar-6)/3;
                    obj.Cee(1+o1,1+o1) = obj.sigma2_XR(2);
                    obj.Cee(1+o1*2,1+o1*2) = obj.sigma2_XR(3);
            end
            %velocities
            switch(mode)
                case {2,3}
                    o1 = obj.nPar/3;
                    obj.Cee(2:o1:obj.nPar,2:o1:obj.nPar) = obj.sigmaq0.vel;
                case {4,5}
                    o1 = (obj.nPar-6)/3;
                    obj.Cee(2:o1:(obj.nPar-6),2:o1:(obj.nPar-6)) = obj.sigmaq0.vel;
            end
            %acceleration
            switch(mode)
                case 3
                    o1 = obj.nPar/3;
                    obj.Cee(3:o1:obj.nPar,3:o1:obj.nPar) = obj.sigmaq0.acc;
                case 5
                    o1 = (obj.nPar-6)/3;
                    obj.Cee(3:o1:(obj.nPar-6),3:o1:(obj.nPar-6)) = obj.sigmaq0.acc;
            end
            
            switch(mode)
                case {4,5}
                    % angular attitude
                    obj.Cee((obj.nPar-5):(obj.nPar-3),(obj.nPar-5):(obj.nPar-3)) = obj.sigmaq0.ang;
                    %angulare velocities
                    obj.Cee((obj.nPar-23):obj.nPar,(obj.nPar-2):obj.nPar) = obj.sigmaq0.ang_vel;
            end
            %initial ambiguities
            obj.Cee(obj.nPar+1:end,obj.nPar+1:end) = diag(obj.sigma2_N(:));
        end
        
        % doppler-based prediction of phase ranges
        function init_doppler(obj, goObs)
            
            %extract the doppler observations at the first epoch
            %for the rover
            %reshape to have nSat rows, nRec columns, nFreq planes
            dop_R = reshape(goObs.getGNSSdop_R(goObs.idGPS, 0, 0, 1, 0),obj.nSat,nRec,nFreq);
            %for the master
            dop_M = goObs.getGNSSdop_M(goObs.idGPS, 0, 1, 0);
            
            %extract phase observations at the first epoch
            %for the rover
            %reshape to have nSat rows, nRec columns, nFreq planes
            ph_R = reshape(goObs.getGNSSph_R(goObs.idGPS, 0, 0, 1, 0),obj.nSat,nRec,nFreq);
            %for the master
            ph_M = goObs.getGNSSph_M(goObs.idGPS, 0, 1, 0);
            
            for r=1:nRec
                if (dop_R(obj.goodSat_pr,r,1))
                    obj.doppler_pred_range_R(obj.goodSat_pr,r) = ph_R(obj.goodSat_pr,r,1) - dop_R(obj.goodSat_pr,r,1);
                end
                if (dop_R(obj.goodSat_pr,r,2))
                    obj.doppler_pred_range_R(obj.goodSat_pr,r) = ph_R(obj.goodSat_pr,r,2) - dop_R(obj.goodSat_pr,r,2);
                end
            end
            if (dop_M(obj.goodSat_pr,1))
                obj.doppler_pred_range_M(obj.goodSat_pr,1) = ph_M(obj.goodSat_pr,1) - dop_M(obj.goodSat_pr,1);
            end
            if (dop_M(obj.goodSat_pr,2))
                obj.doppler_pred_range_M(obj.goodSat_pr,1) = ph_M(obj.goodSat_pr,2) - dop_M(obj.goodSat_pr,2);
            end
        end
        
        % initial kalman filter DOP (DILUTION OF PRECISION)
        function init_KxDOP(obj, mode)
            
            %covariance propagation
            switch(mode)
                case {1,2,3}
                    o1 = obj.nPar/3;
                    obj.Cee_XYZ = obj.Cee(1:o1:obj.nPar,1:o1:obj.nPar);
                    obj.Cee_ENU = global2localCov(obj.Cee_XYZ, obj.Xhat_t_t(1:o1:obj.nPar));
                case {4,5}
                    o1 = (obj.nPar-6)/3;
                    obj.Cee_XYZ = obj.Cee(1:o1:(obj.nPar-6),1:o1:(obj.nPar-6));
                    obj.Cee_ENU = global2localCov(obj.Cee_XYZ, obj.Xhat_t_t(1:o1:(obj.nPar-6)));
            end
            
            %KF DOP computation
            obj.xDOP.KP = sqrt(obj.Cee_XYZ(1,1) + obj.Cee_XYZ(2,2) + obj.Cee_XYZ(3,3));
            obj.xDOP.KH = sqrt(obj.Cee_ENU(1,1) + obj.Cee_ENU(2,2));
            obj.xDOP.KV = sqrt(obj.Cee_ENU(3,3));
        end
        
    end
    
    methods (Access = 'private')
        
    end
    
end

