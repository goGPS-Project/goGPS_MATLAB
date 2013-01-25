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
        satCoordR; % for each receiver: azimuth (az), elevation (el), distance (dist)
        satCoordM; % azimuth (az), elevation (el), distance (dist)
        
        % DILUTION OF PRECISION
        xDOP;  % P, H, V, KP, KH, KV
        
        %Doppler-predicted range (ROVER and MASTER)
        doppler_pred_range_R;
        doppler_pred_range_M;
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
                case 1, obj.nPar = 3;
                case 2, obj.nPar = 6;
                case {3,4}, obj.nPar = 9;
                case 5, obj.nPar = 12;
            end
            obj.mode = mode;
            obj.interval = 1/sampling_rate;	% Init estimation sampling rate
            obj.setDefaultVariances();  % Init variances
            obj.setCurrentParameters(); % Init current parameters
            obj.allocateMemory(goObs.getNumRec(), goObs.getGNSSNumFreq(goObs.idGPS)); % only GPS observations
        end
        
        % Distructor (Shiva)
        function delete(obj)
        end
    end
    
    methods (Access = 'private')
        % Function to initialize all the variances used in the KF
        function setDefaultVariances(obj)
            %variance of initial state
            obj.sigmaq0.pos = 9;
            obj.sigmaq0.vel = 0;
            obj.sigmaq0.ang = 0;
            
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
        
        % Function to fill KF matrices
        function init(obj, goObs)
            obj.init_T(obj.mode);
            obj.init_Xhat_t_t_R(goObs);
        end
        
        % Function to fill Transition matrix according to the number of parameters
        function init_T(obj, mode)
            nPar = obj.nPar;
            % Transition matrix filling
            switch(mode)
                case 1,
                    obj.T(1:nPar,1:nPar) = eye(nPar);
                case 2,
                    obj.T(1:nPar,1:nPar) = eye(nPar);
                    obj.T(1:2:nPar,2:2:nPar) = diag(ones(nPar-3,1)*obj.interval);
                case 3,
                    obj.T(1:nPar,1:nPar) = eye(nPar);
                    obj.T(1:3:nPar,2:3:nPar) = diag(ones(nPar-6,1)*obj.interval);
                    obj.T(2:3:nPar,3:3:nPar) = diag(ones(nPar-6,1)*obj.interval);
                case 4,
                    obj.T(1:nPar,1:nPar) = eye(nPar);
                    obj.T(1:2:nPar-3,2:2:nPar-3) = diag(ones(nPar-6,1)*obj.interval);
                    obj.T(nPar-2:nPar,nPar-2:nPar) = diag(ones(nPar-6,1)*obj.interval);
                case 5,
                    obj.T(1:nPar,1:nPar) = eye(nPar);
                    obj.T(1:2:nPar-3,2:2:nPar-3) = diag(ones(nPar-9,1)*obj.interval);
                    obj.T(nPar-2:nPar,nPar-2:nPar) = diag(ones(nPar-9,1)*obj.interval);
            end
            % note that: the remaining part of the T matrix has already
            % been created as an identity matrix sized nN.
        end
        
        % initialization of the parameter vector for all receivers
        function init_Xhat_t_t_R(obj, goObs)   %% to initialize Xhat_t_t_R: cell of [nPar+nSat*nFreq,1];
            nSat = obj.nSat;
            nRec = goObs.getNumRec();
            nPar = obj.nPar;
            nN = obj.nN;
            nFreq = goObs.getGNSSNumFreq(goObs.idGPS);
            
            % define logical matrices for the satellites in view
            commonSat_pr = logical(zeros(nSat,nRec));
            commonSat_pr_ph = logical(zeros(nSat,nRec));
            
            % select only the satellites in common between master and rover
            if (nFreq == 2) % double frequency
                % reshape the pr and ph observations to have nSat rows,
                % nRec columns, nFreq planes
                pr_R = reshape(goObs.getGNSSpr_R(goObs.idGPS, 0, 0, 1, 0),nSat,nRec,nFreq);
                ph_R = reshape(goObs.getGNSSph_R(goObs.idGPS, 0, 0, 1, 0),nSat,nRec,nFreq);
                pr_M = goObs.getGNSSpr_M(goObs.idGPS, 0, 1, 0);
                ph_M = goObs.getGNSSph_M(goObs.idGPS, 0, 1, 0);
                
                commonSat_pr = (pr_R(:,:,1) ~= 0) & repmat((pr_M(:,:,1) ~= 0),1,nRec) & (pr_R(:,:,2) ~= 0) & repmat((pr_M(:,:,2) ~= 0),1,nRec);
                commonSat_pr_ph = (pr_R(:,:,1) ~= 0) & repmat((pr_M(:,:,1) ~= 0),1,nRec) & (ph_R(:,:,1) ~= 0) & repmat((ph_M(:,:,1) ~= 0),1,nRec) & ...
                    (pr_R(:,:,2) ~= 0) & repmat((pr_M(:,:,2) ~= 0),1,nRec) & (ph_R(:,:,2) ~= 0) & repmat((ph_M(:,:,2) ~= 0),1,nRec);
            else
                if (nFreq == 1) % single frequency
                    sat_pr = (pr_R(:,:,1) ~= 0) & repmat((pr_M(:,:,1) ~= 0),1,nRec);
                    sat = (pr_R(:,:,1) ~= 0) & repmat((pr_M(:,:,1) ~= 0),1,nRec) & ...
                        (ph_R(:,:,1) ~= 0) & repmat((ph_M(:,:,1) ~= 0),1,nRec);
                else
                    % to be used for nFreq>2
                end
            end
            %            end
            
            %------------------------------------------------------------------------------------
            % APPROXIMATE POSITION for the rovers and the master
            %-----------------------------------------------------------------------------------
            
            [XR0 flag_XR] = goObs.getPos_R(0);  %[sized 3xnRec]
            [XM0 flag_M] = goObs.getPos_M();    %[sized 3x1]
            
            %--------------------------------------------------------------------------------------------
            % KALMAN FILTER INITIAL STATE
            %--------------------------------------------------------------------------------------------
            
            %%zeroes vectors useful in the matrices definition
            %Z_om_1 = zeros(o1-1,1);
            
            %variances of the phase ambiguity estimate
            sigma2_N = zeros(nN,1);
            
            %define a vector for the SP3 ephemerides (loaded as a structure)
            ephSP3 = goObs.getSP3();
            
            %logical indexes of the remaining satellites after cutoff
            goodSat_pr_M = logical(zeros(nSat,1));
            goodSat_pr_R = logical(zeros(nSat,nRec));
            
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
             obj.satCoordM(sat_pr_M).el, obj.satCoordM(sat_pr_M).az, obj.satCoordM(sat_pr_M).dist, ...
             cov_XM, var_dtM] ...
             = init_positioning(goObs.getTime_Ref(), pr_M(sat_pr_M_init,1,1), goObs.getGNSSsnr_M(goObs.idGPS, sat_pr_M_init,1,1), ...
             goObs.getGNSSeph(goObs.idGPS), ephSP3.time, ephSP3.coor, ephSP3.clck, ...
             goObs.getIono(), XM0, [], [], sat_pr_M_init, obj.cutoff, obj.snr_threshold, flag_M, 0);
            
         goodSat_pr_M(sat_pr_M) = 1;
            
            %having at least 4 satellites in common in view
            %if (sum(commonSat_pr) >= 4)
            for r=1:nRec
                %having at least 4 satellites in view from the master
                %station after applying the cutoff
                if (sum(goodSat_pr_M) < 4); return; end
                % select only the satellites in view for the rover
                % after the Master cutoff
                % _init because it is before applying the rover cutoff
                % (it is a logical vector)
                sat_pr_R_init(:,r) = (goodSat_pr_M ~= 0) & (commonSat_pr(:,r) ~= 0);
                
                [XR(:,r), dtR(r), ~, ~, ~, ~, ~, err_tropo_R(:,r), err_iono_R(:,r), sat_pr_R, obj.satCoordR(sat_pr_R,r).el, obj.satCoordR(sat_pr_R,r).az, obj.satCoordR(sat_pr_R,r).dist, cov_XR(:,:,r), var_dtR(r), PDOP(r), HDOP(r), VDOP(r), cond_num(r)] = init_positioning(goObs.getTime_Ref(), pr_R(sat_pr_R_init(:,r),r,1), goObs.getGNSSsnr_R(goObs.idGPS, sat_pr_R_init(:,r),r,1), goObs.getGNSSeph(goObs.idGPS), ephSP3.time, ephSP3.coor, ephSP3.clck, goObs.getIono(), XR0(:,r), [], [], sat_pr_R_init(:,r), obj.cutoff, obj.snr_threshold, flag_XR(r), 1);
                goodSat_pr_R(sat_pr_R,r) = 1;
            end
            %keep only satellites that rover and master have in common
            %over the cutoff threshold
            goodSat_pr = goodSat_pr_R;
            
            % atmospheric errors of the receivers for the satellites in
            % common between the master and each rover (nRec columns)
            err_tropo_R = err_tropo_R*goodSat_pr;
            err_iono_R  = err_iono_R*goodSat_pr;
            err_tropo_M = repmat(err_tropo_M,1,nRec)*goodSat_pr;
            err_iono_M  = repmat(err_tropo_M,1,nRec)*goodSat_pr;
            
            %apply cutoffs also to phase satellites
            goodSat_pr_ph = commonSat_pr_ph & goodSat_pr;
            
            % fill doppler variables
            %for i = 1:sum(goodSat_pr,r)
            if (~isempty(goObs.getClockDrift_M()) & goObs.getGNSSdop_M(goObs.idGPS, goodSat_pr_M, 1, 1, 1) == 0 & any(goObs.getGNSSeph(goObs.idGPS)))
                satObs = goObs.getSatObservation(goObs.idGPS, goodSat_pr_M);
                [goObs.getGNSSdop_M(goObs.idGPS, goodSat_pr_M, 1, 1, 2), goObs.getGNSSdop_M(goObs.idGPS, goodSat_pr_M, 1, 1, 1)] ...
                    = doppler_shift_approx(goObs.getPos_M(), zeros(3,1), ...
                                           satObs.X(goodSat_pr_M,:)', satObs.V(goodSat_pr_M,:)', satObs.time(goodSat_pr_M,1), ...
                                           goObs.getClockDrift_M(), goodSat_pr_M, goObs.getGNSSeph(goObs.idGPS));
            end
            
            
            %--------------------------------------------------------------------------------------------
            % SATELLITE CONFIGURATION SAVING AND PIVOT SELECTION
            %--------------------------------------------------------------------------------------------
            
            %satellites configuration: code only (-1), both code and phase (+1);
            conf_sat = zeros(nSat,nRec);
            conf_sat(goodSat_pr) = -1;
            conf_sat(goodSat_pr_ph) = +1;
            
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
                if ~isempty(goodSat_pr_ph(:,r))
                    % find the index for the most elevated satellite (PIVOT)
                    [null_max_elR, pivot_index(r)] = max(obj.satCoordR(goodSat_pr_ph(:,r),r).el);
                    s_id = 1:nSat;                          % all the satellites
                    s_id = s_id(goodSat_pr_ph(:,r));        % extract available satellites
                    obj.pivot(r) = s_id(pivot_index(r));       % get the pivot satellite
                else %only code observations
                    [null_max_elR, pivot_index(r)] = max(obj.satCoordR(goodSat_pr(:,r),r).el);
                    s_id = 1:obj.nSat;                      % all the satellites
                    s_id = s_id(goodSat_pr(:,r));           % extract available satellites
                    obj.pivot(r) = s_id(pivot_index(r));       % get the pivot satellite
                end
                
                %if at least 4 satellites are available after the cutoffs, and if the
                % condition number in the least squares does not exceed the threshold
                if (sum(goodSat_pr,r) >= 4 & cond_num(r) < obj.cond_num_threshold)
                    
                    if isempty(cov_XR(:,:,r)) %if it was not possible to compute the covariance matrix
                        cov_XR(:,:,r) = obj.sigmaq0 * eye(3);
                    end
                    sigma2_XR(r) = diag(cov_XR(:,:,r));
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
            
            for r=1:nRec
                if (sum(goodSat_pr(:,r)) + sum(goodSat_pr_ph(:,r)) - 2 <= 3 + sum(goodSat_pr_ph(:,r)) - 1 | sum(goodSat_pr_ph(:,r)) <= 4)
                    
                    %computation of the phase double differences in order to estimate N
                    if ~isempty(goodSat_pr_ph(:,r))
                        [N1(goodSat_pr_ph(:,r),r), sigma2_N1(goodSat_pr_ph(:,r),r)] = amb_estimate_observ(pr_R(goodSat_pr_ph(:,r),r,1), pr_M(goodSat_pr_ph(:,r),1,1), ph_R(goodSat_pr_ph(:,r),r,1), ph_M(goodSat_pr_ph(:,r),1,1), obj.pivot(r), goodSat_pr_ph(:,r), 1);
                        [N2(goodSat_pr_ph(:,r),r), sigma2_N2(goodSat_pr_ph(:,r),r)] = amb_estimate_observ(pr_R(goodSat_pr_ph(:,r),r,2), pr_M(goodSat_pr_ph(:,r),1,2), ph_R(goodSat_pr_ph(:,r),r,2), ph_M(goodSat_pr_ph(:,r),1,2), obj.pivot(r), goodSat_pr_ph(:,r), 2);
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
                    if ~isempty(goodSat_pr_ph(:,r))
                        [     XR(:,r), N1(goodSat_pr_ph(:,r),r),      cov_XR(:,:,r), cov_N1(:,:,r), PDOP(r), HDOP(r), VDOP(r)] ...
                        = LS_DD_code_phase(XR(:,r), XM, XS(goodSat_pr_ph,:), ...
                                           pr_R(goodSat_pr_ph(:,r),r,1), ph_R(goodSat_pr_ph(:,r),r,1), snr_R(goodSat_pr_ph(:,r),r,1), ...
                                           pr_M(goodSat_pr_ph(:,r),1,1), ph_M(goodSat_pr_ph(:,r),1,1), snr_M(goodSat_pr_ph(:,r),1,1), ...
                                           obj.satCoordR(goodSat_pr_ph(:,r),r).el, obj.satCoordM(goodSat_pr_ph(:,r),1).el, ...
                                           err_tropo_R(goodSat_pr_ph(:,r),r), err_iono_R(goodSat_pr_ph(:,r),r), ...
                                           err_tropo_M(goodSat_pr_ph(:,r),r), err_iono_M(goodSat_pr_ph(:,r),r), ...
                                           pivot_index(r), nFreq);
                        [null_XR, N2(goodSat_pr_ph(:,r),r), null_cov_XR, cov_N2(:,:,r)] ...
                        = LS_DD_code_phase(XR(:,r), XM, XS(goodSat_pr_ph,:), ...
                                           pr_R(goodSat_pr_ph(:,r),r,2), ph_R(goodSat_pr_ph(:,r),r,2), snr_R(goodSat_pr_ph(:,r),r,2), ...
                                           pr_M(goodSat_pr_ph(:,r),1,2), ph_M(goodSat_pr_ph(:,r),1,2), snr_M(goodSat_pr_ph(:,r),1,2), ...
                                           obj.satCoordR(goodSat_pr_ph(:,r),r).el, obj.satCoordM(goodSat_pr_ph(:,r),1).el, ...
                                           err_tropo_R(goodSat_pr_ph(:,r),r), err_iono_R(goodSat_pr_ph(:,r),r), ...
                                           err_tropo_M(goodSat_pr_ph(:,r),r), err_iono_M(goodSat_pr_ph(:,r),r), ...
                                           pivot_index(r), nFreq);
                    end
                    
                    if isempty(cov_XR(:,:,r)) %if it was not possible to compute the covariance matrix
                        cov_XR(:,:,r) = obj.sigmaq0 * eye(3);
                    end
                    sigma2_XR(:,r) = diag(cov_XR(:,:,r));
                    
                    if isempty(cov_N1(:,:,r)) %if it was not possible to compute the covariance matrix
                        cov_N1(:,:,r) = obj.sigmaq0.N * eye(length(goodSat_pr_ph(:,r),r));
                    end
                    
                    if isempty(cov_N2(:,:,r)) %if it was not possible to compute the covariance matrix
                        cov_N2(:,:,r) = obj.sigmaq0.N * eye(length(goodSat_pr_ph(:,r),r));
                    end
                    
                    if (nFreq == 2)
                        N(:,r) = [N1(:,r); N2(:,r)];
                        sigma2_N(goodSat_pr_ph(:,r),r) = diag(cov_N1(:,:,r));
                        sigma2_N(nSat+goodSat_pr_ph(:,r),r) = diag(cov_N2(:,:,r));
                    else
                        if (nFreq == 1)
                            N(:,r) = N1(:,r);
                            sigma2_N(goodSat_pr_ph(:,r),r) = diag(:,:,r);
                        else
                            % to be used for nFreq > 2
                        end
                    end
                end
                %initialization of the initial point with 6(positions and velocities) +
                %32 or 64 (N combinations) variables
                Xhat_t_t = [XR(1); Z_om_1; XR(2); Z_om_1; XR(3); Z_om_1; N];
            end
        end
    end
end