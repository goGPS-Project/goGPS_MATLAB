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
        % Initialization variances
        sigmaq0;
        
        % User defined variances
        sigmaq;
        
        % Minimum number of satellites to be used in the Kalman filter
        min_nsat = 2;
        
        % Cut-off [degrees]
        cutoff = 0;
        
        % Signal-to-noise ratio threshold [dB]
        snr_threshold = 0;
        
        % Cycle slip threshold [cycles]
        cs_threshold = 10;
        
        % Parameter used to select the weight mode for GPS observations
        %          - weights=0: same weight for all the observations
        %          - weights=1: weight based on satellite elevation
        %          - weights=2: weight based on signal-to-noise ratio
        %          - weights=3: weight based on combined elevation and signal-to-noise ratio
        weights = 2;
        
        % Weight function parameters
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
        
        % Order of the dynamic model polynomial
        order = 2;
        
        % Number of parameters to be estimated
        nPar = 3;
        
        % Ambiguity restart method
        amb_restart_method = 1;
        
        % Transition matrix
        T = [];
        
        % Identity matrix
        I = [];
        
        % Receiver(s) coordinates
        XR = [];
        
        % Tmp state estimation at time t for each receiver
        Xhat_t_t_R = {};
        
        % State estimation at time t
        Xhat_t_t = [];
        
        % State estimation at time t+1 (using dynamics only)
        X_t1_t = [];
        
        % Estimation error covariance matrix at time t
        Cee = [];
        
        % Tmp estimated coordinates variance for each receiver
        sigma2_XR_R = [];
        
        % Estimated coordinates variance for the cluster
        sigma2_XR = [];
        
        % Estimated variance for the initial ambiguities
        sigma2_N
        
        % Satellite configuration (1: visible, 0: not visible) at time t:
        conf_sat = [];
        
        % Previous configuration (1: visible, 0: not visible) at time t:
        conf_sat_old = [];
        
        % Cycle-slip configuration (1: cs, 0: no cs) at time t
        conf_cs = [];
        
        % Index of the current pivot satellite for each receiver
        pivot = [];
        
        % Index of the previous pivot satellite for each receiver
        pivot_old = [];
        
        % Number of unknown phase ambiguities
        nN = [];
        
        % Method used to estimate phase ambiguities
        %          - amb_estim_method=0: observed code - phase difference
        %          - amb_estim_method=1: Kalman-predicted code - phase difference
        %          - amb_estim_method=2: Least squares adjustment
        amb_estim_method = 0;
        
        % Estimation sampling rate
        interval = 1; %default 1 Hz (to avoid problems with real-time modes)
        
        % Threshold on the condition number on the eigenvalues of the N matrix (least squares)
        cond_num_threshold = 1e6;
        
        % Azimuth, elevation and distance of satellites with respect to the
        % ROVER and MASTER
        satCoordR = struct('az',zeros(32,1),'el',zeros(32,1),'dist',zeros(32,1)); % for each receiver: azimuth (az), elevation (el), distance (dist)
        satCoordM = struct('az',zeros(32,1),'el',zeros(32,1),'dist',zeros(32,1)); % azimuth (az), elevation (el), distance (dist)
        
        % Satellites in common between master and rover, over the cutoff
        goodSat_pr;     %only code
        goodSat_pr_ph;  %phase and code
        
        % DILUTION OF PRECISION
        xDOP;  % P, H, V, KP, KH, KV
        
        % Doppler-predicted range (ROVER and MASTER)
        doppler_pred_range_R;
        doppler_pred_range_M;
        
        % Flag containing initialization status
        initKF = false;
    end
    
    methods
        % Creator (Brahma)
        function obj = goKalmanFilter(goObs, goIni, mode, sampling_rate)
            if nargin < 4
                sampling_rate = 1;
            end
            if nargin < 3
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
                case 5, obj.nPar = 12;  % const.acceleration filter + attitude angles without variations
            end
            obj.mode = mode;
            obj.interval = 1/sampling_rate;	% Init estimation sampling rate
            obj.setDefaultVariances();      % Init variances
            obj.setCurrentParameters();     % Init current parameters
            obj.allocateMemory(goObs.getNumRec(), goObs.getGNSSnFreq(goGNSS.ID_GPS)); % only GPS observations
            
            obj.init(goObs, goIni)
        end
        
        % Destructor (Shiva)
        function delete(obj)
        end
        
        % return the status of the KF initialization
        function isI = isInitialized(obj)
            isI = obj.initKF;
        end
    end
    
    % Function to fill KF initial matrices
    methods (Access = 'public')
        function init(obj, goObs, goIni)
            obj.init_T(obj.mode);
            obj.init_Xhat_t_t(goObs, goIni, obj.mode);
            obj.init_X_t1_t(); 
            obj.init_Cee(goObs.getNumRec(), obj.mode);
       %     %obj.init_doppler(goObs, goObs.getNumRec(), goObs.getGNSSnFreq(goGNSS.ID_GPS));
       %     obj.init_KxDOP(obj.mode);
            obj.MR_loop(goObs, goIni, obj.mode);
       %     % set the status of the KF initialization
       %     obj.initKF = true;
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
            nP = obj.nPar;
            nSat = goGNSS.MAX_SAT;
            obj.nN = (nSat*nFreq*nRec);                 %number of unknown phase ambiguities
            obj.T = eye(nP+obj.nN);                   %transition matrix
            obj.I = eye(nP+obj.nN);                   %identity matrix
            obj.XR = zeros(3,nRec);                   %receiver(s) coordinates
            obj.Xhat_t_t_R = cell(nRec,1);              %tmp state estimation at time t for each receiver
            for r=1:(nRec),
                obj.Xhat_t_t_R{r} = zeros(nP+nSat*nFreq,1);
            end
            obj.Xhat_t_t = zeros(nP+obj.nN,1);        %state estimation at time t
            obj.X_t1_t = zeros(nP+obj.nN,1);          %state estimation at time t+1 (using dynamics only)
            obj.Cee = zeros(nP+obj.nN);               %estimation error covariance matrix at time t
            obj.sigma2_XR_R = zeros(3,nRec);            %estimation positioning variance for each receiver
            obj.sigma2_XR_R = zeros(3,1);               %estimation positioning variance for the cluster
            obj.sigma2_N = zeros(nSat*nFreq,nRec);      %variances of the phase ambiguity estimate
            obj.conf_sat = zeros(nSat);                 %DD satellite configuration (1: visible, 0: not visible) at time t, 
            obj.conf_sat_old = zeros(nSat);                 %DD satellite configuration (1: visible, 0: not visible) at time t-1, 
            obj.conf_cs = zeros(nSat, nRec);            %cycle-slip configuration (1: cs, 0: no cs) at time t
            obj.pivot = zeros(1,nRec);                  %index of the current pivot satellite
            obj.pivot_old = zeros(1,nRec);              %index of the previous pivot satellite
            
            % spherical coordinates initialization
            
            % pre-allocate the structure for the vectors obj.satCoordR(1|2|3).az,
            % obj.satCoordR(1|2|3).el obj.satCoordR(1|2|3).dist
            obj.satCoordR = struct('az',zeros(nSat,nRec),'el',zeros(nSat,nRec),'dist',zeros(nSat,nRec));
            obj.satCoordM = struct('az',zeros(nSat,1),'el',zeros(nSat,1),'dist',zeros(nSat,1));
            
            % satellites in common between master and rover, over the cutoff
            obj.goodSat_pr = zeros(nSat,nRec);
            obj.goodSat_pr_ph = zeros(nSat,nRec);
            
            % DILUTION OF PRECISION
            obj.xDOP = struct('P', zeros(nRec,1), 'H', zeros(nRec,1), 'V', zeros(nRec,1), 'KP', zeros(nRec,1), 'KH', zeros(nRec,1), 'KV', zeros(nRec,1));
            
            % Doppler-predicted range
            obj.doppler_pred_range_R = zeros(nSat,nRec,nFreq);
            obj.doppler_pred_range_M = zeros(nSat,1,nFreq);
        end
        
        % Function to fill Transition matrix according to the number of parameters
        function init_T(obj, mode)
            nP = obj.nPar;
            % Transition matrix filling
            switch(mode)
                case 1,
                    % obj.T(1:nPar,1:nPar) = eye(nPar);
                    % redundant, becuse T is already an identity matrix
                case 2,
                    % obj.T(1:nPar,1:nPar) = eye(nPar);
                    obj.T(1:2:nP,2:2:nP) = diag(ones(3,1)*obj.interval);
                case 3,
                    % obj.T(1:nPar,1:nPar) = eye(nPar);
                    obj.T(1:3:nP,2:3:nP) = diag(ones(3,1)*obj.interval);
                    obj.T(2:3:nP,3:3:nP) = diag(ones(3,1)*obj.interval);
                case 4,
                    % obj.T(1:nPar,1:nPar) = eye(nPar);
                    obj.T(1:2:(nP-6),2:2:(nP-6)) = diag(ones(3,1)*obj.interval);
                case 5,
                    % obj.T(1:nPar,1:nPar) = eye(nPar);
                    obj.T(1:3:(nP-3),2:3:(nP-3)) = diag(ones(3,1)*obj.interval);
                    obj.T(2:3:(nP-3),3:3:(nP-3)) = diag(ones(3,1)*obj.interval);
            end
            % note that: the remaining part of the T matrix has already
            % been created as an identity matrix sized nN.
        end
        
        % initialization of the parameter vector for all receivers
        function init_Xhat_t_t(obj, goObs, goIni, mode)   %% to initialize Xhat_t_t_R: cell of [nPar+nSat*nFreq,1] and Xhat_t_t;

            goObs.doPreProcessing()
            
            %load of the observation data after pre-processing
            
            nP = obj.nPar;
            
            ID_GNSS=1; % <- must be taken from the object!
            
            nRec=goObs.getNumRec();
            
            time_GPS=goObs.getTime_Ref();
            
            nFreq=goObs.getGNSSnFreq(ID_GNSS);
            Eph=goObs.getGNSSeph(ID_GNSS);
            iono=goObs.getIono();
            SP3_time=goObs.getGNSS_SP3time();
            SP3_coor=goObs.getGNSS_SP3coordinates();
            SP3_clck=goObs.getGNSS_SP3clock();
            
    
            
            global Xhat_t_t  % forse conviene aggiungere output alla funzione goGPS_LS_DD_code_phase
            
            
            global azR azM elR elM distR distM;
            satCoord = struct('az',zeros(goGNSS.MAX_SAT,nRec+1),'el',zeros(goGNSS.MAX_SAT,nRec+1),'dist',zeros(goGNSS.MAX_SAT,nRec+1)); %first column: MASTER, further columns: ROVERS
            err_iono = NaN(goGNSS.MAX_SAT, nRec+1);
            err_tropo = NaN(goGNSS.MAX_SAT, nRec+1);
            
           
            for t = 1 : 1                
                pr1_M=goObs.getGNSSpr_M(ID_GNSS,0,t,1);   %pr = getGNSSpr_M(obj, idGNSS, idSat, idObs, nFreq)
                ph1_M=goObs.getGNSSph_M(ID_GNSS,0,t,1);   %ph = getGNSSph_M(obj, idGNSS, idSat, idObs, nFreq)
                snr_M=goObs.getGNSSsnr_M(ID_GNSS,0,t,1);  %snr = getGNSSsnr_M(obj, idGNSS, idSat, idObs, nFreq)
                dop1_M=goObs.getGNSSdop_M(ID_GNSS,0,t,1); %dop = getGNSSdop_M(obj, idGNSS, idSat, idObs, nFreq)
                
                
                pr2_M=zeros(size(pr1_M));
                ph2_M=zeros(size(pr1_M));
                dop2_M=zeros(size(pr1_M));
                if nFreq==2
                    pr2_M=goObs.getGNSSpr_M(ID_GNSS,0,t,2);   %pr = getGNSSpr_M(obj, idGNSS, idSat, idObs, nFreq)
                    ph2_M=goObs.getGNSSph_M(ID_GNSS,0,t,2);   %ph = getGNSSph_M(obj, idGNSS, idSat, idObs, nFreq)
                    dop2_M=goObs.getGNSSdop_M(ID_GNSS,0,t,2); %dop = getGNSSdop_M(obj, idGNSS, idSat, idObs, nFreq)
                end
                

                pr1_R=zeros(goGNSS.MAX_SAT,nRec);
                ph1_R=zeros(goGNSS.MAX_SAT,nRec);
                snr_R=zeros(goGNSS.MAX_SAT,nRec);
                dop1_R=zeros(goGNSS.MAX_SAT,nRec);
                pr2_R=zeros(goGNSS.MAX_SAT,nRec);
                ph2_R=zeros(goGNSS.MAX_SAT,nRec);
                dop2_R=zeros(goGNSS.MAX_SAT,nRec);
                
                pr1_R=reshape(goObs.getGNSSpr_R(ID_GNSS,0,0,t,1),goGNSS.MAX_SAT,nRec);       %pr = getGNSSpr_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
                ph1_R=reshape(goObs.getGNSSph_R(ID_GNSS,0,0,t,1),goGNSS.MAX_SAT,nRec);       %ph = getGNSSph_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
                snr_R=reshape(goObs.getGNSSsnr_R(ID_GNSS,0,0,t,1),goGNSS.MAX_SAT,nRec);      %snr = getGNSSsnr_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
                dop1_R=reshape(goObs.getGNSSdop_R(ID_GNSS,0,0,t,1),goGNSS.MAX_SAT,nRec);     %dop = getGNSSdop_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
                
                if nFreq==2
                    pr2_R=reshape(goObs.getGNSSpr_R(ID_GNSS,0,0,t,2),goGNSS.MAX_SAT,nRec);       %pr = getGNSSpr_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
                    ph2_R=reshape(goObs.getGNSSph_R(ID_GNSS,0,0,t,2),goGNSS.MAX_SAT,nRec);       %ph = getGNSSph_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
                    dop2_R=reshape(goObs.getGNSSdop_R(ID_GNSS,0,0,t,2),goGNSS.MAX_SAT,nRec);     %dop = getGNSSdop_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
                end
                
                
                
                %  QUESTO VA FATTO ADESSO? PENSO DI SI'
                %if (~flag_SP3)  <-- sistemare il flag_SP3 prendenolo dall'obj
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
                dop1_R(delsat,:) = 0;
                dop1_M(delsat,:) = 0;
                dop2_R(delsat,:) = 0;
                dop2_M(delsat,:) = 0;
                snr_R(delsat,:) = 0;
                snr_M(delsat,:) = 0;
                %end
                
                
                
                % Processing
                % ----------
                %  - for each rover receiver:
                %       - enhance coordinates with code and phase DD in single
                %         epoch (lambda)
                %  - estimation of apriori attitude
                %  - enhance solution with constrained least squares (DD in
                %         single epoch (lambda))
                
                
                % for each rover receiver: enhance coordinates with code and phase DD in single epoch (lambda)
                % --------------------------------------------------------------------------------------------
                
                % cosa sono???
                check_on = 0;
                check_off = 0;
                check_pivot = 0;
                check_cs = 0;
                plot_t = 1;
                %
                
               
                
                

                %cartesian to geodetic conversion of MASTER coordinates
                [pos_M flag_M]= goObs.getPos_M(t);
                [phiM, lamM, hM] = cart2geod(pos_M(1,t), pos_M(2,t), pos_M(3,t));
                
                %radians to degrees
                phiM = phiM * 180 / pi;
                lamM = lamM * 180 / pi;
                
                
                % Eph_t = Eph(:,:,t);
                Eph_t=rt_find_eph(goObs.getGNSSeph(goGNSS.ID_GPS), goObs.getTime_Ref(t));
                               
                X_sat=zeros(goGNSS.MAX_SAT,3,nRec);
                %conf_sat=zeros(goGNSS.MAX_SAT,nRec);
                
                for i=1:nRec
                    statistic = zeros(2,length(time_GPS)); % <-- VA DIMENSIONATO IN 3D!!!!?
                    ambiguity = 0;                         % <-- VA DIMENSIONATO IN 3D!!!!?
                    
                    Xhat_t_t=zeros(size(Xhat_t_t));
                    
                    [X_sat_i obj.goodSat_pr_ph(:,i)]=goGPS_LS_DD_code_phase(time_GPS(t), pos_M, pr1_R(:,i), pr1_M, pr2_R(:,i), pr2_M, ph1_R(:,i), ph1_M, ph2_R(:,i), ph2_M, snr_R(:,i), snr_M, Eph_t, SP3_time, SP3_coor, SP3_clck, iono, 1);
                    
                    X_sat(find(obj.goodSat_pr_ph(:,i)==1),1:3,i)=X_sat_i;
                    
                    obj.XR(1,i)=Xhat_t_t(1);
                    obj.XR(1,i)=Xhat_t_t(3);
                    obj.XR(1,i)=Xhat_t_t(5);
                    
                    % compute elevation and atmospheric corrections from XR_DD coordinates
                    [phiR, lamR, hR] = cart2geod(obj.XR(1,i), obj.XR(2,i), obj.XR(3,i));
                    phiR = phiR * 180 / pi;
                    lamR = lamR * 180 / pi;
                    
                    distM(~distM)=NaN;
                    azM(isnan(distM))=NaN;
                    elM(isnan(distM))=NaN;
                    azR(isnan(distM))=NaN;
                    elR(isnan(distM))=NaN;
                    distR(isnan(distM))=NaN;
                    
                    obj.satCoordM.az(:,1)= azM;
                    obj.satCoordM.el(:,1)= elM;
                    obj.satCoordM.dist(:,1)= distM;
                    
                    obj.satCoordR.az(:,i)= azR;
                    obj.satCoordR.el(:,i)= elR;
                    obj.satCoordR.dist(:,i)= distR;
                    
                    %computation of atmospheric errors of Master receiver
                    err_tropo(:,1) = tropo_error_correction(obj.satCoordM.el(:,1), hM);
                    err_iono(:,1) = iono_error_correction(phiM, lamM, obj.satCoordM.az(:,1), obj.satCoordM.el(:,1), goObs.getTime_Ref(t), goObs.getIono(), goObs.getSBAS());
                    err_iono(isnan(distM),1)=NaN;
                    
                    %computation of atmospheric errors of Rover receiver
                    err_tropo(:,i+1) = tropo_error_correction(obj.satCoordR.el(:,i), hR);
                    
                    %computation of ionospheric errors
                    err_iono(:,i+1) = iono_error_correction(phiR, lamR, obj.satCoordR.az(:,i), obj.satCoordR.el(:,i), goObs.getTime_Ref(t), goObs.getIono(), goObs.getSBAS());
                    err_iono(isnan(distM),i+1)=NaN;
                    
                end
                
            end
            
            
            % leftin zeros the apriori attitude
            attitude_approx=[0 0 0]';
            
            
            % instrumental RS coordinates
            % ---------------------------
            % get geometry
            [geometry ev_point]=goIni.getGeometry();
            
            % barycenter definition
            xb=mean(geometry,2);
            % barycentric instrumental RS coordinates
            xR=geometry-repmat(xb,1,nRec);
            
            % non barycentric! the origin is the first receiver
            %xR=geometry;
            
            
            % estimated local coordinates of baycenter
            %-----------------------------------------
            Xb_apriori=[mean(obj.XR(1,:)), mean(obj.XR(2,:)),mean(obj.XR(3,:))]';
            
            
            
            
            % code + phase double differenecs with Xb and attitude
            % ----------------------------------------------------
            obj.conf_sat=ones(goGNSS.MAX_SAT,1);
            
            for i=1:nRec
                obj.conf_sat=obj.conf_sat(:,1)&obj.goodSat_pr_ph(:,i);
            end
            sat=find(obj.conf_sat==1); %only common sat! it would be better to use all receiver independently
            
            
            XS=NaN(goGNSS.MAX_SAT,3);
            XS(sat,1:3)=X_sat(sat,1:3,1);
            
            %actual pivot
            [null_max_elR pivot_index]=max(obj.satCoordM.el(sat,1));
            pivot_r = sat(pivot_index);
            
            
            
            %% compute diff --------- must be put outside to avoid the recomputation every epoch
            syms s_Xb s_Yb s_Zb s_phi_b s_lam_b s_roll s_pitch s_yaw s_x s_y s_z s_XS s_YS s_ZS s_XS_Piv s_YS_Piv s_ZS_Piv s_XM s_YM s_ZM
            syms r11 r12 r13 r21 r22 r23 r31 r32 r33
            
            % rotation from local to global
            s_Rgl=[-sin(s_lam_b) cos(s_lam_b) 0; ...
                -sin(s_phi_b)*cos(s_lam_b) -sin(s_phi_b)*sin(s_lam_b) cos(s_phi_b); ...
                cos(s_phi_b)*cos(s_lam_b) cos(s_phi_b)*sin(s_lam_b) sin(s_phi_b)];
            
            
            % rotation from instrumental to local
            s_Rli=[cos(s_roll)*cos(s_pitch) -sin(s_pitch)*cos(s_yaw)+cos(s_roll)*sin(s_pitch)+sin(s_yaw) sin(s_roll)*sin(s_yaw)+cos(s_roll)*sin(s_pitch)*cos(s_yaw); ...
                sin(s_roll)+cos(s_pitch) sin(s_roll)*sin(s_pitch)*sin(s_yaw) sin(s_roll)*sin(s_pitch)*cos(s_yaw)-cos(s_roll)*sin(s_yaw); ...
                -sin(s_pitch) cos(s_pitch)*sin(s_yaw) cos(s_pitch)*cos(s_yaw)];
            
            % rotation from instrumental to local
            %s_Rli=[r11 r12 r13; r21 r22 r23; r31 r32 r33];
            
            
            s_X=[s_Xb; s_Yb; s_Zb]+s_Rgl*s_Rli*[s_x;s_y;s_z];
            
            s_PR_DD=sqrt((s_X(1)-s_XS)^2+(s_X(2)-s_YS)^2 + (s_X(3)-s_ZS)^2) - sqrt((s_XM-s_XS)^2+(s_YM-s_YS)^2 + (s_ZM-s_ZS)^2) - ...
                (sqrt((s_X(1)-s_XS_Piv)^2+(s_X(2)-s_YS_Piv)^2 + (s_X(3)-s_ZS_Piv)^2) - sqrt((s_XM-s_XS_Piv)^2+(s_YM-s_YS_Piv)^2 + (s_ZM-s_ZS_Piv)^2));
            
            
            s_Ai=[diff(s_PR_DD,s_Xb) diff(s_PR_DD,s_Yb) diff(s_PR_DD,s_Zb) diff(s_PR_DD,s_roll) diff(s_PR_DD,s_pitch) diff(s_PR_DD,s_yaw)];
            %s_Ai=[diff(s_PR_DD,s_Xb) diff(s_PR_DD,s_Yb) diff(s_PR_DD,s_Zb) diff(s_PR_DD,r11) diff(s_PR_DD,r12) diff(s_PR_DD,r13) diff(s_PR_DD,r21) diff(s_PR_DD,r22) diff(s_PR_DD,r23) diff(s_PR_DD,r31) diff(s_PR_DD,r32) diff(s_PR_DD,r33)];
            
            
            
            F_Ai=inline(s_Ai); %% <-- colonne della matrice disegno corrispondenti alle incognite geometriche
            F_PR_DD=inline(s_PR_DD); %% <-- parte geometrica dei termini noti
            F_s_X=inline(s_X); %% <-- per calcolare le coordinate aggiornate di ogni ricevitore a partire da quelle nuove del punto 'baricentro'
            %%
            
            
            index_sat_without_pivot=sat;
            index_sat_without_pivot(pivot_index)=[];
            
            if (size(sat,1) >= 4) % & cond_num < cond_num_threshold)
                phase_1=1;
                
                %loop is needed to improve the atmospheric error correction
                for i = 1 : 3
                    %if (phase == 1)
                    [Xb_apriori, N1_hat, cov_Xb, cov_N1, cov_ATT, attitude_approx, obj.XR, PDOP, HDOP, VDOP] = LS_DD_code_phase_MR(Xb_apriori, obj.XR, pos_M(:,t), XS(sat,:), pr1_R(sat,:), ph1_R(sat,:), snr_R(sat,:), pr1_M(sat), ph1_M(sat,t), snr_M(sat,t), obj.satCoordR.el(sat,:), obj.satCoordM.el(sat,1), err_tropo(sat,2:nRec+1), err_iono(sat,2:nRec+1), err_tropo(sat,1), err_iono(sat,1), pivot_index, phase_1, attitude_approx, geometry, 0, F_Ai, F_PR_DD, F_s_X);
                    %else
                    %    [XR, N1_hat, cov_XR, cov_N1, PDOP, HDOP, VDOP, up_bound, lo_bound, posType] = LS_DD_code_phase(XR, XM, XS, pr2_R(sat), ph2_R(sat), snr_R(sat), pr2_M(sat), ph2_M(st), snr_M(sat), elR(sat), elM(sat), err_tropo_R, err_iono_R, err_tropo_M, err_iono_M, pivot_index, phase);
                    %end
                    
                    % compute elevation and atmospheric corrections from XR_DD coordinates
                    for r=1:nRec
                        [phiR, lamR, hR] = cart2geod(obj.XR(1,r), obj.XR(2,r), obj.XR(3,r));
                        phiR = phiR * 180 / pi;
                        lamR = lamR * 180 / pi;
                        [obj.satCoordR.az(sat,r),obj.satCoordR.el(sat,r), obj.satCoordR.dist(sat,r)] = topocent(obj.XR(:,r), XS(sat,:));
                        
                        %computation of atmospheric errors of Rover receiver
                        err_tropo(:,r+1) = tropo_error_correction(obj.satCoordR.el(:,r), hR);
                        
                        %computation of ionospheric errors
                        err_iono(:,r+1) = iono_error_correction(phiR, lamR, obj.satCoordR.az(:,r), obj.satCoordR.el(:,r), goObs.getTime_Ref(t), goObs.getIono(), goObs.getSBAS());
                        err_iono(isnan(distM),r+1)=NaN;        
                    end
                end
                
                
            else
            end

            
            %interfacing with old version
            %----------------------------
            % initialize pivot
            obj.pivot = zeros(1,nRec);
            %previous pivot
            obj.pivot_old = zeros(1,nRec);
            N1=zeros(goGNSS.MAX_SAT,nRec);  
            obj.sigma2_XR=diag(cov_Xb);
            
            
            cov_N1=diag(cov_N1);
            for r=1:nRec
                obj.pivot(r)=pivot_r;            
                N1(index_sat_without_pivot,r)=N1_hat((r-1)*length(index_sat_without_pivot)+1:r*length(index_sat_without_pivot));
                obj.sigma2_N(index_sat_without_pivot,r)=cov_N1((r-1)*length(index_sat_without_pivot)+1:r*length(index_sat_without_pivot));
                obj.xDOP.P(r)=PDOP(1,r);
                obj.xDOP.H(r)=HDOP(1,r);
                obj.xDOP.V(r)=VDOP(1,r);
            end
            
             switch(mode)
                 %case 1,
                 %    obj.Xhat_t_t(1:3) = [obj.XR(1,1); obj.XR(2,1); obj.XR(3,1)];
                 %case {2,4},
                 %    obj.Xhat_t_t(1:6) = [obj.XR(1,1); 0; obj.XR(2,1); 0; obj.XR(3,1); 0];
                 case {3,5},
                     obj.Xhat_t_t(1:9) = [Xb_apriori(1,1); 0; 0; Xb_apriori(2,1); 0; 0; Xb_apriori(3,1); 0; 0];
             end
                          
             switch(mode)
                 %case {1,2,3},   %when not estimating the attitude
                 %    obj.Xhat_t_t (nP+1:end) = N(:);
                 case {4,5},     % when estimating the attitude (roll, pitch, yaw angles)
                     %attitude = goObs.getInitialAttitude();
                     obj.Xhat_t_t (nP-3+1:nP) = [attitude_approx(1); attitude_approx(2); attitude_approx(3)];
                     obj.Xhat_t_t (nP+1:end) = N1(:); %% <-- sistemare per la doppia frequenza!
             end
             
%              return
%              
%             
%             
%             
%             
% 
%             %--------------------------------------------------------------------------------------------
%             % KALMAN FILTER INITIAL STATE
%             %--------------------------------------------------------------------------------------------
% 
%             %[X, dt, usableSat, satCoord, cov_X, var_dt, PDOP, HDOP, VDOP, cond_num] = goGNSS.LS_MR_C_SA(goObs, goIni);
%             
%             %--------------------------------------------------------------------------------------------
%             % SATELLITE CONFIGURATION SAVING AND PIVOT SELECTION
%             %--------------------------------------------------------------------------------------------
%             
%             %satellites configuration: code only (-1), both code and phase (+1);
%             conf_sat = zeros(nSat,nR);
%             conf_sat = obj.goodSat_pr*-1;
%             conf_sat(obj.goodSat_pr_ph) = ones(sum(obj.goodSat_pr_ph(:)),1);
%             
%             %cycle-slip configuration (no cycle-slip)
%             conf_cs = zeros(nSat,nR);
%             
%             % initialize pivot
%             obj.pivot = zeros(1,nR);
%             pivot_index = zeros(1,nR);
%             %previous pivot
%             obj.pivot_old = zeros(1,nR);
%             
%             %current pivot
%             for r=1:nR
%                 % if there is at least one visible satellite for each
%                 % receiver (with both code and phase obervations)
%                 if ~isempty(obj.goodSat_pr_ph(:,r))
%                     % find the index for the most elevated satellite (PIVOT)
%                     [null_max_elR, pivot_index(r)] = max(obj.satCoordR.el(obj.goodSat_pr_ph(:,r),r));
%                     s_id = 1:nSat;                          % all the satellites
%                     s_id = s_id(obj.goodSat_pr_ph(:,r));        % extract available satellites
%                     obj.pivot(r) = s_id(pivot_index(r));       % get the pivot satellite
%                 else %only code observations
%                     [null_max_elR, pivot_index(r)] = max(obj.satCoordR.el(obj.goodSat_pr(:,r),r));
%                     s_id = 1:goGNSS.MAX_SAT;                      % all the satellites
%                     s_id = s_id(obj.goodSat_pr(:,r));           % extract available satellites
%                     obj.pivot(r) = s_id(pivot_index(r));       % get the pivot satellite
%                 end
%                 
%                 %if at least 4 satellites are available after the cutoffs, and if the
%                 % condition number in the least squares does not exceed the threshold
%                 if (sum(obj.goodSat_pr(:,r)) >= 4 && cond_num(r) < obj.cond_num_threshold)
%                     
%                     if isempty(cov_XR(:,:,r)) %if it was not possible to compute the covariance matrix
%                         cov_XR(:,:,r) = obj.sigmaq0 * eye(3);
%                     end
%                     obj.sigma2_XR_R(:,r) = diag(cov_XR(:,:,r));
%                 else
%                     return
%                 end
%             end
%             
%             %do not use least squares ambiguity estimation
%             % NOTE: LS amb. estimation is automatically switched off if the number of
%             % satellites with phase available is not sufficient
%             
%             %ambiguity initialization: initialized value if the satellite is visible,
%             %0 if the satellite is not visible
%             N1 = zeros(nSat,nR);
%             N2 = zeros(nSat,nR);
%             sigma2_N1 = zeros(nSat,nR);
%             sigma2_N2 = zeros(nSat,nR);
%             N = zeros(nSat, nR);
%             for r=1:nR
%                 if (sum(obj.goodSat_pr(:,r)) + sum(obj.goodSat_pr_ph(:,r)) - 2 <= 3 + sum(obj.goodSat_pr_ph(:,r)) - 1 || sum(obj.goodSat_pr_ph(:,r)) <= 4)
%                     
%                     %computation of the phase double differences in order to estimate N
%                     if ~isempty(obj.goodSat_pr_ph(:,r))
%                         [N1(obj.goodSat_pr_ph(:,r),r), sigma2_N1(obj.goodSat_pr_ph(:,r),r)] = amb_estimate_observ(pr_R(obj.goodSat_pr_ph(:,r),r,1), pr_M(obj.goodSat_pr_ph(:,r),1,1), ph_R(obj.goodSat_pr_ph(:,r),r,1), ph_M(obj.goodSat_pr_ph(:,r),1,1), obj.pivot(r), obj.goodSat_pr_ph(:,r), 1);
%                         [N2(obj.goodSat_pr_ph(:,r),r), sigma2_N2(obj.goodSat_pr_ph(:,r),r)] = amb_estimate_observ(pr_R(obj.goodSat_pr_ph(:,r),r,2), pr_M(obj.goodSat_pr_ph(:,r),1,2), ph_R(obj.goodSat_pr_ph(:,r),r,2), ph_M(obj.goodSat_pr_ph(:,r),1,2), obj.pivot(r), obj.goodSat_pr_ph(:,r), 2);
%                     end
%                     
%                     if (nFreq == 2)
%                         N(:,r) = [N1(:,r); N2(:,r)];
%                         sigma2_N(:,r) = [sigma2_N1(:,r); sigma2_N2(:,r)];
%                     else
%                         if (nFreq == 1)
%                             N(:,r) = N1(:,r);
%                             sigma2_N(:,r) = sigma2_N1(:,r);
%                         else
%                             % to be used for nFreq > 2
%                         end
%                     end
%                     
%                     %use least squares ambiguity estimation
%                 else
%                     
%                     %ROVER positioning improvement with code and phase double differences
%                     if ~isempty(obj.goodSat_pr_ph(:,r))
%                         if nFreq == 1
%                             [obj.XR(:,r), N1(obj.goodSat_pr_ph(:,r),r), cov_XR(:,:,r), cov_N1(:,:,r), obj.xDOP.P(r), obj.xDOP.H(r), obj.xDOP.V(r)] ...
%                                 = LS_DD_code_phase(obj.XR(:,r), XM, XS, ...
%                                 pr_R(obj.goodSat_pr_ph(:,r),r,1), ph_R(obj.goodSat_pr_ph(:,r),r,1), goObs.getGNSSsnr_R(goGNSS.ID_GPS, obj.goodSat_pr_ph(:,r), r, 1, 1), ...
%                                 pr_M(obj.goodSat_pr_ph(:,r),1,1), ph_M(obj.goodSat_pr_ph(:,r),1,1), goObs.getGNSSsnr_M(goGNSS.ID_GPS, obj.goodSat_pr_ph(:,r), 1, 1), ...
%                                 obj.satCoordR.el(obj.goodSat_pr_ph(:,r),r), obj.satCoordM.el(obj.goodSat_pr_ph(:,r),1), ...
%                                 err_tropo_R(obj.goodSat_pr_ph(:,r),r), err_iono_R(obj.goodSat_pr_ph(:,r),r), ...
%                                 err_tropo_M(obj.goodSat_pr_ph(:,r),r), err_iono_M(obj.goodSat_pr_ph(:,r),r), ...
%                                 pivot_index(r), nFreq);
%                         elseif nFreq == 2
%                             [obj.XR(:,r), N1(obj.goodSat_pr_ph(:,r),r), cov_XR(:,:,r), cov_N1(:,:,r), obj.xDOP.P(r), obj.xDOP.H(r), obj.xDOP.V(r)] ...
%                                 = LS_DD_code_phase(obj.XR(:,r), XM, XS, ...
%                                 pr_R(obj.goodSat_pr_ph(:,r),r,1), ph_R(obj.goodSat_pr_ph(:,r),r,1), goObs.getGNSSsnr_R(goGNSS.ID_GPS, obj.goodSat_pr_ph(:,r), r, 1, 1), ...
%                                 pr_M(obj.goodSat_pr_ph(:,r),1,1), ph_M(obj.goodSat_pr_ph(:,r),1,1), goObs.getGNSSsnr_M(goGNSS.ID_GPS, obj.goodSat_pr_ph(:,r), 1, 1), ...
%                                 obj.satCoordR.el(obj.goodSat_pr_ph(:,r),r), obj.satCoordM.el(obj.goodSat_pr_ph(:,r),1), ...
%                                 err_tropo_R(obj.goodSat_pr_ph(:,r),r), err_iono_R(obj.goodSat_pr_ph(:,r),r), ...
%                                 err_tropo_M(obj.goodSat_pr_ph(:,r),r), err_iono_M(obj.goodSat_pr_ph(:,r),r), ...
%                                 pivot_index(r), nFreq);
%                             [null_XR, N2(obj.goodSat_pr_ph(:,r),r), null_cov_XR, cov_N2(:,:,r)] ...
%                                 = LS_DD_code_phase(obj.XR(:,r), XM, XS, ...
%                                 pr_R(obj.goodSat_pr_ph(:,r),r,2), ph_R(obj.goodSat_pr_ph(:,r),r,2), goObs.getGNSSsnr_R(goGNSS.ID_GPS, obj.goodSat_pr_ph(:,r), r, 1, 2), ...
%                                 pr_M(obj.goodSat_pr_ph(:,r),1,2), ph_M(obj.goodSat_pr_ph(:,r),1,2), goObs.getGNSSsnr_M(goGNSS.ID_GPS, obj.goodSat_pr_ph(:,r), 1, 2), ...
%                                 obj.satCoordR(obj.goodSat_pr_ph(:,r),r).el, obj.satCoordM(obj.goodSat_pr_ph(:,r),1).el, ...
%                                 err_tropo_R(obj.goodSat_pr_ph(:,r),r), err_iono_R(obj.goodSat_pr_ph(:,r),r), ...
%                                 err_tropo_M(obj.goodSat_pr_ph(:,r),r), err_iono_M(obj.goodSat_pr_ph(:,r),r), ...
%                                 pivot_index(r), nFreq);
%                         end
%                     end
%                     
%                     if isempty(cov_XR(:,:,r)) %if it was not possible to compute the covariance matrix
%                         cov_XR(:,:,r) = obj.sigmaq0 * eye(3);
%                     end
%                     obj.sigma2_XR_R(:,r) = diag(cov_XR(:,:,r));
%                     
%                     if isempty(cov_N1(:,:,r)) %if it was not possible to compute the covariance matrix
%                         cov_N1(:,:,r) = obj.sigmaq0.N * eye(length(obj.goodSat_pr_ph(:,r),r));
%                     end
%                     
%                     if (nFreq == 2)
%                         if isempty(cov_N2(:,:,r)) %if it was not possible to compute the covariance matrix
%                             cov_N2(:,:,r) = obj.sigmaq0.N * eye(length(obj.goodSat_pr_ph(:,r),r));
%                         end
%                     end
%                     
%                     if (nFreq == 2)
%                         N(:,r) = [N1(:,r); N2(:,r)];
%                         obj.sigma2_N(obj.goodSat_pr_ph(:,r),r) = diag(cov_N1(:,:,r));
%                         obj.sigma2_N(nSat+obj.goodSat_pr_ph(:,r),r) = diag(cov_N2(:,:,r));
%                     else
%                         if (nFreq == 1)
%                             N(:,r) = N1(:,r);
%                             obj.sigma2_N(obj.goodSat_pr_ph(:,r),r) = diag(cov_N1(:,:,r));
%                         else
%                             % to be used for nFreq > 2
%                         end
%                     end
%                 end
%                 %initialization of the initial point with 3/6/9(positions/velocities/accelerations) +
%                 % 3 attitude angles + 3 attitude variations +
%                 %32 or 64 (N combinations) variables
%                 switch(mode)
%                     case 1,
%                         obj.Xhat_t_t_R{r} = [obj.XR(1,r); obj.XR(2,r); obj.XR(3,r); N(:,r)];
%                     case {2,4},
%                         obj.Xhat_t_t_R{r} = [obj.XR(1,r); 0; obj.XR(2,r); 0; obj.XR(3,r); 0; N(:,r)];
%                     case {3,5},
%                         obj.Xhat_t_t_R{r} = [obj.XR(1,r); 0; 0; obj.XR(2,r); 0; 0; obj.XR(3,r); 0; 0; N(:,r)];
%                 end
%             end
%             %in the state vector the coordinates of the baricenter are
%             %considered
%             
%             % obj.XR contains in the first column the baricenter coordinates,
%             % then the coordinates of each receiver
%             obj.XR = [mean(obj.XR,2) obj.XR];
%             
%             %if there is only one receiver, delete the baricenter column
%             %(useless)
%             if nR == 1;
%                 obj.XR(:,1) = [];
%             end
%             
%             %build the state vector
%             switch(mode)
%                 case 1,
%                     obj.Xhat_t_t(1:3) = [obj.XR(1,1); obj.XR(2,1); obj.XR(3,1)];
%                 case {2,4},
%                     obj.Xhat_t_t(1:6) = [obj.XR(1,1); 0; obj.XR(2,1); 0; obj.XR(3,1); 0];
%                 case {3,5},
%                     obj.Xhat_t_t(1:9) = [obj.XR(1,1); 0; 0; obj.XR(2,1); 0; 0; obj.XR(3,1); 0; 0];
%             end
%             
%             switch(mode)
%                 case {1,2,3},   %when not estimating the attitude
%                     obj.Xhat_t_t (nP+1:end) = N(:);
%                 case {4,5},     % when estimating the attitude (roll, pitch, yaw angles)
%                     attitude = goObs.getInitialAttitude();
%                     obj.Xhat_t_t (nP-6+1:nP) = [attitude.roll; 0; attitude.pitch; 0; attitude.yaw; 0];
%                     obj.Xhat_t_t (nP+1:end) = N(:);
%             end
        end
        
        % initialization of point estimation at step t+1 ==
        % estimation at step t, because the initial velocity is equal to 0
        function init_X_t1_t(obj)
            X_t1_t = obj.T*obj.Xhat_t_t;
        end
        
        % initialization of state covariance matrix
        function init_Cee(obj, nRec, mode)
            %obj.sigma2_XR = sum((1/nRec)^2*obj.sigma2_XR_R);     % variance propagation
            %positions
            %obj.Cee(1,1) = obj.sigma2_XR(1);  
            Cee_diag=zeros(length(obj.Cee),1);
            Cee_diag(1) = obj.sigma2_XR(1);
            switch(mode)
                case {1,2,3}
                    o1 = obj.nPar/3;
                    obj.Cee(1+o1,1+o1) = obj.sigma2_XR(2);
                    obj.Cee(1+o1*2,1+o1*2) = obj.sigma2_XR(3);
                case 4
                    o1 = (obj.nPar-6)/3;
                    obj.Cee(1+o1,1+o1) = obj.sigma2_XR(2);
                    obj.Cee(1+o1*2,1+o1*2) = obj.sigma2_XR(3);
                case 5
                    o1 = (obj.nPar-3)/3;
                    %obj.Cee(1+o1,1+o1) = obj.sigma2_XR(2);
                    %obj.Cee(1+o1*2,1+o1*2) = obj.sigma2_XR(3);  
                    Cee_diag(1+o1) = obj.sigma2_XR(2);
                    Cee_diag(1+o1*2) = obj.sigma2_XR(3); 
            end
            %velocities
            switch(mode)
                case {2,3}
                    o1 = obj.nPar/3;
                    obj.Cee(2:o1:obj.nPar,2:o1:obj.nPar) = obj.sigmaq0.vel;
                case 4
                    o1 = (obj.nPar-6)/3;
                    obj.Cee(2:o1:(obj.nPar-6),2:o1:(obj.nPar-6)) = obj.sigmaq0.vel;
                case 5
                    o1 = (obj.nPar-3)/3;
                    %obj.Cee(2:o1:(obj.nPar-3),2:o1:(obj.nPar-3)) = obj.sigmaq0.vel;
                    Cee_diag(2:o1:(obj.nPar-3)) = obj.sigmaq0.vel;
            end
            %acceleration
            switch(mode)
                case 3
                    o1 = obj.nPar/3;
                    obj.Cee(3:o1:obj.nPar,3:o1:obj.nPar) = obj.sigmaq0.acc;
                case 5
                    o1 = (obj.nPar-3)/3;
                    %obj.Cee(3:o1:(obj.nPar-3),3:o1:(obj.nPar-3)) = obj.sigmaq0.acc;
                    Cee_diag(3:o1:(obj.nPar-3)) = obj.sigmaq0.acc;
            end
            
            switch(mode)
                case 4
                    % angular attitude
                    obj.Cee((obj.nPar-5):(obj.nPar-3),(obj.nPar-5):(obj.nPar-3)) = obj.sigmaq0.ang;
                    %angulare velocities
                    obj.Cee((obj.nPar-3):obj.nPar,(obj.nPar-2):obj.nPar) = obj.sigmaq0.ang_vel;
                case 5
                    % angular attitude
                    %obj.Cee((obj.nPar-2):(obj.nPar),(obj.nPar-2):(obj.nPar)) = obj.sigmaq0.ang;
                    Cee_diag((obj.nPar-2):(obj.nPar)) = obj.sigmaq0.ang;
            end
            %initial ambiguities
            %obj.Cee(obj.nPar+1:end,obj.nPar+1:end) = diag(obj.sigma2_N(:));
            Cee_diag(obj.nPar+1:end) = obj.sigma2_N(:);
            obj.Cee=diag(Cee_diag);
        end
        
        %         % doppler-based prediction of phase ranges
        %         function init_doppler(obj, goObs, nRec, nFreq)
        %
        %             %             %extract the doppler observations at the first epoch
        %             %             %for the rover
        %             %             %reshape to have nSat rows, nRec columns, nFreq planes
        %             %             dop_R = reshape(goObs.getGNSSdop_R(goGNSS.ID_GPS, 0, 0, 1, 0),goGNSS.MAX_SAT,nRec,nFreq);
        %             %             %for the master
        %             %             dop_M = goObs.getGNSSdop_M(goGNSS.ID_GPS, 0, 1, 0);
        %             %
        %             %             %extract phase observations at the first epoch
        %             %             %for the rover
        %             %             %reshape to have nSat rows, nRec columns, nFreq planes
        %             %             ph_R = reshape(goObs.getGNSSph_R(goGNSS.ID_GPS, 0, 0, 1, 0),goGNSS.MAX_SAT,nRec,nFreq);
        %             %             %for the master
        %             %             ph_M = goObs.getGNSSph_M(goGNSS.ID_GPS, 0, 1, 0);
        %             %
        %             %             for r=1:nRec
        %             %                 if (dop_R(obj.goodSat_pr(:,r),r,1))
        %             %                     obj.doppler_pred_range_R(obj.goodSat_pr(:,r),r) = ph_R(obj.goodSat_pr(:,r),r,1) - dop_R(obj.goodSat_pr(:,r),r,1);
        %             %                 end
        %             %                 if (dop_R(obj.goodSat_pr(:,r),r,2))
        %             %                     obj.doppler_pred_range_R(obj.goodSat_pr(:,r),r) = ph_R(obj.goodSat_pr(:,r),r,2) - dop_R(obj.goodSat_pr(:,r),r,2);
        %             %                 end
        %             %             end
        %             %             if (dop_M(obj.goodSat_pr(:,r),1))
        %             %                 obj.doppler_pred_range_M(obj.goodSat_pr(:,r),1) = ph_M(obj.goodSat_pr(:,r),1) - dop_M(obj.goodSat_pr(:,r),1);
        %             %             end
        %             %             if (dop_M(obj.goodSat_pr(:,r),2))
        %             %                 obj.doppler_pred_range_M(obj.goodSat_pr(:,r),1) = ph_M(obj.goodSat_pr(:,r),2) - dop_M(obj.goodSat_pr(:,r),2);
        %             %             end
        %         end
        
        % initial kalman filter DOP (DILUTION OF PRECISION)
        function init_KxDOP(obj, mode)
            
            %covariance propagation
            switch(mode)
                case {1,2,3}
                    o1 = obj.nPar/3;
                    Cee_XYZ = obj.Cee(1:o1:obj.nPar,1:o1:obj.nPar);
                    Cee_ENU = global2localCov(Cee_XYZ, obj.Xhat_t_t(1:o1:obj.nPar));
                case {4,5}
                    o1 = (obj.nPar-6)/3;
                    Cee_XYZ = obj.Cee(1:o1:(obj.nPar-6),1:o1:(obj.nPar-6));
                    Cee_ENU = global2localCov(Cee_XYZ, obj.Xhat_t_t(1:o1:(obj.nPar-6)));
            end
            
            %KF DOP computation
            obj.xDOP.KP = sqrt(Cee_XYZ(1,1) + Cee_XYZ(2,2) + Cee_XYZ(3,3));
            obj.xDOP.KH = sqrt(Cee_ENU(1,1) + Cee_ENU(2,2));
            obj.xDOP.KV = sqrt(Cee_ENU(3,3));
        end
        
        %Loop functions
        function MR_loop(obj, goObs, goIni, mode)
            
            nP = obj.nPar;
            
            ID_GNSS=1; % <- must be taken from the object!
            
            nRec=goObs.getNumRec();
            
            time_GPS=goObs.getTime_Ref();
            
            nFreq=goObs.getGNSSnFreq(ID_GNSS);
            Eph=goObs.getGNSSeph(ID_GNSS);
            iono=goObs.getIono();
            SP3_time=goObs.getGNSS_SP3time();
            SP3_coor=goObs.getGNSS_SP3coordinates();
            SP3_clck=goObs.getGNSS_SP3clock();
            
            global Xhat_t_t  % forse conviene aggiungere output alla funzione goGPS_LS_DD_code_phase
            
            global azR azM elR elM distR distM sigmaq0_N;
           
            
%             satCoord = struct('az',zeros(goGNSS.MAX_SAT,nRec+1),'el',zeros(goGNSS.MAX_SAT,nRec+1),'dist',zeros(goGNSS.MAX_SAT,nRec+1)); %first column: MASTER, further columns: ROVERS
            err_iono = NaN(goGNSS.MAX_SAT, nRec+1);
            err_tropo = NaN(goGNSS.MAX_SAT, nRec+1);
                     
                       
            % compute diff --------- must be put outside to avoid the recomputation every epoch
            % ---------------------------------------------------------------------------------
            syms s_Xb s_Yb s_Zb s_phi_b s_lam_b s_roll s_pitch s_yaw s_x s_y s_z s_XS s_YS s_ZS s_XS_Piv s_YS_Piv s_ZS_Piv s_XM s_YM s_ZM
            syms r11 r12 r13 r21 r22 r23 r31 r32 r33
            
            % rotation from local to global
            s_Rgl=[-sin(s_lam_b) cos(s_lam_b) 0; ...
                -sin(s_phi_b)*cos(s_lam_b) -sin(s_phi_b)*sin(s_lam_b) cos(s_phi_b); ...
                cos(s_phi_b)*cos(s_lam_b) cos(s_phi_b)*sin(s_lam_b) sin(s_phi_b)];            
            
            % rotation from instrumental to local
            s_Rli=[cos(s_roll)*cos(s_pitch) -sin(s_pitch)*cos(s_yaw)+cos(s_roll)*sin(s_pitch)+sin(s_yaw) sin(s_roll)*sin(s_yaw)+cos(s_roll)*sin(s_pitch)*cos(s_yaw); ...
                sin(s_roll)+cos(s_pitch) sin(s_roll)*sin(s_pitch)*sin(s_yaw) sin(s_roll)*sin(s_pitch)*cos(s_yaw)-cos(s_roll)*sin(s_yaw); ...
                -sin(s_pitch) cos(s_pitch)*sin(s_yaw) cos(s_pitch)*cos(s_yaw)];
            
            % rotation from instrumental to local
            %s_Rli=[r11 r12 r13; r21 r22 r23; r31 r32 r33];            
            
            s_X=[s_Xb; s_Yb; s_Zb]+s_Rgl*s_Rli*[s_x;s_y;s_z];
            
            s_PR_DD=sqrt((s_X(1)-s_XS)^2+(s_X(2)-s_YS)^2 + (s_X(3)-s_ZS)^2) - sqrt((s_XM-s_XS)^2+(s_YM-s_YS)^2 + (s_ZM-s_ZS)^2) - ...
                (sqrt((s_X(1)-s_XS_Piv)^2+(s_X(2)-s_YS_Piv)^2 + (s_X(3)-s_ZS_Piv)^2) - sqrt((s_XM-s_XS_Piv)^2+(s_YM-s_YS_Piv)^2 + (s_ZM-s_ZS_Piv)^2));
            
            s_Ai=[diff(s_PR_DD,s_Xb) diff(s_PR_DD,s_Yb) diff(s_PR_DD,s_Zb) diff(s_PR_DD,s_roll) diff(s_PR_DD,s_pitch) diff(s_PR_DD,s_yaw)];
            %s_Ai=[diff(s_PR_DD,s_Xb) diff(s_PR_DD,s_Yb) diff(s_PR_DD,s_Zb) diff(s_PR_DD,r11) diff(s_PR_DD,r12) diff(s_PR_DD,r13) diff(s_PR_DD,r21) diff(s_PR_DD,r22) diff(s_PR_DD,r23) diff(s_PR_DD,r31) diff(s_PR_DD,r32) diff(s_PR_DD,r33)];
            
            F_Ai=inline(s_Ai); %% <-- colonne della matrice disegno corrispondenti alle incognite geometriche
            F_PR_DD=inline(s_PR_DD); %% <-- parte geometrica dei termini noti
            F_s_X=inline(s_X); %% <-- per calcolare le coordinate aggiornate di ogni ricevitore a partire da quelle nuove del punto 'baricentro'
            % ---------------------------------------------------------------------------------            
            
            
            
            
            
            
            
            %for t = 2 : length(time_GPS)
            for t=2:2
                pr1_M=goObs.getGNSSpr_M(ID_GNSS,0,t,1);   %pr = getGNSSpr_M(obj, idGNSS, idSat, idObs, nFreq)
                ph1_M=goObs.getGNSSph_M(ID_GNSS,0,t,1);   %ph = getGNSSph_M(obj, idGNSS, idSat, idObs, nFreq)
                snr_M=goObs.getGNSSsnr_M(ID_GNSS,0,t,1);  %snr = getGNSSsnr_M(obj, idGNSS, idSat, idObs, nFreq)
                dop1_M=goObs.getGNSSdop_M(ID_GNSS,0,t,1); %dop = getGNSSdop_M(obj, idGNSS, idSat, idObs, nFreq)                
                
                pr2_M=zeros(size(pr1_M));
                ph2_M=zeros(size(pr1_M));
                dop2_M=zeros(size(pr1_M));
                if nFreq==2
                    pr2_M=goObs.getGNSSpr_M(ID_GNSS,0,t,2);   %pr = getGNSSpr_M(obj, idGNSS, idSat, idObs, nFreq)
                    ph2_M=goObs.getGNSSph_M(ID_GNSS,0,t,2);   %ph = getGNSSph_M(obj, idGNSS, idSat, idObs, nFreq)
                    dop2_M=goObs.getGNSSdop_M(ID_GNSS,0,t,2); %dop = getGNSSdop_M(obj, idGNSS, idSat, idObs, nFreq)
                end
                                
                pr1_R=zeros(goGNSS.MAX_SAT,nRec);
                ph1_R=zeros(goGNSS.MAX_SAT,nRec);
                snr_R=zeros(goGNSS.MAX_SAT,nRec);
                dop1_R=zeros(goGNSS.MAX_SAT,nRec);
                pr2_R=zeros(goGNSS.MAX_SAT,nRec);
                ph2_R=zeros(goGNSS.MAX_SAT,nRec);
                dop2_R=zeros(goGNSS.MAX_SAT,nRec);
                
                pr1_R=reshape(goObs.getGNSSpr_R(ID_GNSS,0,0,t,1),goGNSS.MAX_SAT,nRec);       %pr = getGNSSpr_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
                ph1_R=reshape(goObs.getGNSSph_R(ID_GNSS,0,0,t,1),goGNSS.MAX_SAT,nRec);       %ph = getGNSSph_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
                snr_R=reshape(goObs.getGNSSsnr_R(ID_GNSS,0,0,t,1),goGNSS.MAX_SAT,nRec);      %snr = getGNSSsnr_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
                dop1_R=reshape(goObs.getGNSSdop_R(ID_GNSS,0,0,t,1),goGNSS.MAX_SAT,nRec);     %dop = getGNSSdop_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
                
                if nFreq==2
                    pr2_R=reshape(goObs.getGNSSpr_R(ID_GNSS,0,0,t,2),goGNSS.MAX_SAT,nRec);       %pr = getGNSSpr_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
                    ph2_R=reshape(goObs.getGNSSph_R(ID_GNSS,0,0,t,2),goGNSS.MAX_SAT,nRec);       %ph = getGNSSph_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
                    dop2_R=reshape(goObs.getGNSSdop_R(ID_GNSS,0,0,t,2),goGNSS.MAX_SAT,nRec);     %dop = getGNSSdop_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
                end
                
                %  QUESTO VA FATTO ADESSO? PENSO DI SI'
                %if (~flag_SP3)  <-- sistemare il flag_SP3 prendenolo dall'obj
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
                dop1_R(delsat,:) = 0;
                dop1_M(delsat,:) = 0;
                dop2_R(delsat,:) = 0;
                dop2_M(delsat,:) = 0;
                snr_R(delsat,:) = 0;
                snr_M(delsat,:) = 0;
                %end
                

                % Processing
                % ----------
                %  - for each rover receiver:
                %       - enhance coordinates with code and phase DD in single
                %         epoch (lambda)
              
                
                % for each rover receiver: enhance coordinates with code and phase DD in single epoch (lambda)
                % --------------------------------------------------------------------------------------------
                
                % cosa sono???
                check_on = 0;
                check_off = 0;
                check_pivot = 0;
                check_cs = 0;
                plot_t = 1;
                %
                
                %cartesian to geodetic conversion of MASTER coordinates
                [pos_M flag_M]= goObs.getPos_M(t);
                [phiM, lamM, hM] = cart2geod(pos_M(1), pos_M(2), pos_M(3));
                
                %radians to degrees
                phiM = phiM * 180 / pi;
                lamM = lamM * 180 / pi;
                
                
                % Eph_t = Eph(:,:,t);
                Eph_t=rt_find_eph(goObs.getGNSSeph(goGNSS.ID_GPS), goObs.getTime_Ref(t));
                
                X_sat=zeros(goGNSS.MAX_SAT,3,nRec);
                %conf_sat=zeros(goGNSS.MAX_SAT,nRec);
                
                for i=1:nRec
                    statistic = zeros(2,length(time_GPS)); % <-- VA DIMENSIONATO IN 3D!!!!?
                    ambiguity = 0;                         % <-- VA DIMENSIONATO IN 3D!!!!?
                    
                    Xhat_t_t=zeros(size(Xhat_t_t));
                    azR=NaN(goGNSS.MAX_SAT,1);
                    azM=NaN(goGNSS.MAX_SAT,1);
                    elR=NaN(goGNSS.MAX_SAT,1);
                    elM=NaN(goGNSS.MAX_SAT,1);
                    distR=NaN(goGNSS.MAX_SAT,1);
                    distM=NaN(goGNSS.MAX_SAT,1);                 
            
                    
                    [X_sat_i obj.goodSat_pr_ph(:,i)]=goGPS_LS_DD_code_phase(time_GPS(t), pos_M, pr1_R(:,i), pr1_M, pr2_R(:,i), pr2_M, ph1_R(:,i), ph1_M, ph2_R(:,i), ph2_M, snr_R(:,i), snr_M, Eph_t, SP3_time, SP3_coor, SP3_clck, iono, 1);
                    
                    X_sat(find(obj.goodSat_pr_ph(:,i)==1),1:3,i)=X_sat_i;
                    
                    obj.XR(1,i)=Xhat_t_t(1);
                    obj.XR(2,i)=Xhat_t_t(3);
                    obj.XR(3,i)=Xhat_t_t(5);
                    
                    % compute elevation and atmospheric corrections from XR_DD coordinates
                    [phiR, lamR, hR] = cart2geod(obj.XR(1,i), obj.XR(2,i), obj.XR(3,i));
                    phiR = phiR * 180 / pi;
                    lamR = lamR * 180 / pi;
                    
                    distM(~distM)=NaN;
                    azM(isnan(distM))=NaN;
                    elM(isnan(distM))=NaN;
                    azR(isnan(distM))=NaN;
                    elR(isnan(distM))=NaN;
                    distR(isnan(distM))=NaN;
                    
                    obj.satCoordM.az(:,1)= azM;
                    obj.satCoordM.el(:,1)= elM;
                    obj.satCoordM.dist(:,1)= distM;
                    
                    obj.satCoordR.az(:,i)= azR;
                    obj.satCoordR.el(:,i)= elR;
                    obj.satCoordR.dist(:,i)= distR;
                    
                    %computation of atmospheric errors of Master receiver
                    err_tropo(:,1) = tropo_error_correction(obj.satCoordM.el(:,1), hM);
                    err_iono(:,1) = iono_error_correction(phiM, lamM, obj.satCoordM.az(:,1), obj.satCoordM.el(:,1), goObs.getTime_Ref(t), goObs.getIono(), goObs.getSBAS());
                    err_iono(isnan(distM),1)=NaN;
                    
                    %computation of atmospheric errors of Rover receiver
                    err_tropo(:,i+1) = tropo_error_correction(obj.satCoordR.el(:,i), hR);
                    
                    %computation of ionospheric errors
                    err_iono(:,i+1) = iono_error_correction(phiR, lamR, obj.satCoordR.az(:,i), obj.satCoordR.el(:,i), goObs.getTime_Ref(t), goObs.getIono(), goObs.getSBAS());
                    err_iono(isnan(distM),i+1)=NaN;
                    
                end
                
                % get apriori attitude (from previous step)
                switch(mode)
                    %case {1,2,3},   %when not estimating the attitude                    
                    case {5},                        
                        attitude_approx=[obj.Xhat_t_t(nP-3+1:nP)];                        
                end
                
                % instrumental RS coordinates
                % ---------------------------
                % get geometry
                [geometry ev_point]=goIni.getGeometry();
                
                % barycenter definition
                xb=mean(geometry,2);
                % barycentric instrumental RS coordinates
                xR=geometry-repmat(xb,1,nRec);
                
                % non barycentric! the origin is the first receiver
                %xR=geometry;
                
               
                
                
                % code + phase double differenecs with Xb and attitude
                % ----------------------------------------------------
                % old configuration
                obj.conf_sat_old=obj.conf_sat;
                obj.pivot_old=obj.pivot;
                
                % current configuration
                obj.conf_sat=ones(goGNSS.MAX_SAT,1);                
                for i=1:nRec
                    obj.conf_sat=obj.conf_sat(:,1)&obj.goodSat_pr_ph(:,i);
                end
                sat=find(obj.conf_sat==1); %only common sat! it would be better to use all receiver independently                
                
                XS=NaN(goGNSS.MAX_SAT,3);
                XS(sat,1:3)=X_sat(sat,1:3,1);
                
                %actual pivot
                [null_max_elR pivot_index]=max(obj.satCoordM.el(sat,1));
                pivot_r = sat(pivot_index);
                
                index_sat_without_pivot=sat;
                index_sat_without_pivot(pivot_index)=[];
                
                
                nsat=length(sat);
                
                
                %------------------------------------------------------------------------------------
                % LINEARIZATION POINT (APPROXIMATE COORDINATES)
                %------------------------------------------------------------------------------------                
                %approximate position
                switch(mode)
                    case {5},
                        XR0 = obj.Xhat_t_t([1,4,7]);
                        flag_XR = 2;
                        VR0 = obj.Xhat_t_t([2,5,8]);
                end
                
                %approximated coordinates X Y Z
                X_app = XR0(1);
                Y_app = XR0(2);
                Z_app = XR0(3);
                
                %----------------------------------------------------------------------------------------
                % CONVERSION FROM CARTESIAN TO GEODETIC COORDINATES
                %----------------------------------------------------------------------------------------
                [phiR_app, lamR_app, hR_app] = cart2geod(X_app, Y_app, Z_app);
             
                
%                 %----------------------------------------------------------------------------------------
%                 % EXTRACTION OF THE HEIGHT PSEUDO-OBSERVATION FROM THE DTM
%                 %----------------------------------------------------------------------------------------
%                 
%                 %projection to UTM coordinates
%                 [E_app, N_app] = geod2plan(phiR_app, lamR_app);
%                 
%                 %dtm tile detection (in which the approximated position lies)
%                 [tile_row,tile_col] = find ( (E_app > tile_georef(:,:,1)) & (E_app <= tile_georef(:,:,4)) & (N_app >= tile_georef(:,:,3)) & (N_app < tile_georef(:,:,2)));
%                 
%                 %tile buffer dimension
%                 tile_buffer_size = 3;
%                 %check if the approximated position lies within one of the available tiles, otherwise set nodata value
%                 if ( ~isempty(tile_row) & ~isempty(tile_col) )
%                     tile_buffer = cell(tile_buffer_size,tile_buffer_size);
%                     for i = -1 : 1
%                         for j = -1 : 1
%                             %definition of the path and the filename of the selected tile
%                             tile_path = strcat(dtm_dir,'/tiles/tile_',num2str(tile_row+i),'_',num2str(tile_col+j),'.mat');
%                             
%                             %check the existence of the file associated to the selected tile
%                             fid = fopen(tile_path,'r');
%                             if (fid ~= -1)
%                                 fclose(fid);
%                                 %load the selected tile
%                                 load(tile_path, 'tile');
%                             else
%                                 %load of a null tile
%                                 tile(1:tile_header.nrows, 1:tile_header.ncols) = tile_header.nodata;
%                             end
%                             %buffer creation around the selected tile
%                             tile_buffer{i+2,j+2} =  tile;
%                         end
%                     end
%                     
%                     %buffer conversion from cell to matrix
%                     tile_buffer = cell2mat(tile_buffer);
%                     
%                     %computation of the tile buffer dimension (cell number)
%                     [tile_height tile_width] = size(tile_buffer);
%                     
%                     %tile buffer lower left center coordinates extraction
%                     Ell = tile_georef(tile_row,tile_col,1) - tile_width/tile_buffer_size*tile_header.cellsize + tile_header.cellsize/2;
%                     Nll = tile_georef(tile_row,tile_col,3) - tile_height/tile_buffer_size*tile_header.cellsize + tile_header.cellsize/2;
%                     
%                     %extraction from the dtm of the height correspondent to the approximated position
%                     [h_dtm] = grid_bilin_interp(E_app, N_app, tile_buffer, tile_header.ncols*3, tile_header.nrows*3, tile_header.cellsize, Ell, Nll, tile_header.nodata);
%                     
%                     %antenna height addition
%                     h_dtm = h_dtm + h_antenna;
%                 else
                     h_dtm = tile_header.nodata;
%                 end


                %% questa cosa ??
%                 %----------------------------------------------------------------------------------------
%                 % MODEL ERROR COVARIANCE MATRIX
%                 %----------------------------------------------------------------------------------------
%                 
%                 %re-initialization of Cvv matrix of the model error
%                 % (if a static model is used, no noise is added)
%                 Cvv = zeros(o3+nN);
%                 if (o1 > 1)
%                     Cvv(o1,o1) = sigmaq_vE;
%                     Cvv(o2,o2) = sigmaq_vN;
%                     Cvv(o3,o3) = sigmaq_vU;
%                     
%                     %propagate diagonal local cov matrix to global cov matrix
%                     Cvv([o1 o2 o3],[o1 o2 o3]) = local2globalCov(Cvv([o1 o2 o3],[o1 o2 o3]), X_t1_t([1 o1+1 o2+1]));
%                 end
                %%
                
                
                
                %cycle-slip configuration
                obj.conf_cs = zeros(32,nRec);
                
                
                
                min_nsat=4;
                %if the number of available satellites after the cutoffs is equal or greater than min_nsat
                if (nsat >= min_nsat)                    
                    %------------------------------------------------------------------------------------
                    % SATELLITE ADDITION/LOSS
                    %------------------------------------------------------------------------------------                    
                    sat_dead = []; %#ok<NASGU>
                    sat_born = [];
                    sat_old=find(obj.conf_sat_old>0);
                    %search for a lost satellite
                        

                        phase=1; % da settare fuori                      
                        %find lost satellites
                        % if (length(sat) < length(sat_old))   %%% <-------------- perchè??? Se un sat muore e uno nasce, la lunghezza è cmq uguale!
                        sat_dead = setdiff(sat_old,sat);
                        if ~isempty(sat_dead)                             
                            %for lost satellites it is fundamental to set their N-PIVOT
                            % combinations to 0. Furthermore it could be convenient to raise
                            %their uncertainty (not necessary - done when a new satellite is
                            %added)
                            N1 = 0;
                            N2 = 0; 
                            check_off = 1;
                            for i=1:length(sat_dead)                                
%                             if (length(phase) == 2)
%                                 X_t1_t(o3+sat_dead,1) = N1;
%                                 X_t1_t(o3+32+sat_dead,1) = N2;
%                             else
%                                 if (phase == 1)
                                    obj.Xhat_t_t(nP+sat_dead(i):goGNSS.MAX_SAT:end,1) = N1;
%                                 else
%                                     X_t1_t(o3+sat_dead,1) = N2;
%                                 end                                
                            end
                        end
                        
                        
                        %search for a new satellite
                        % if (length(sat) > length(sat_old)) %%% <-------------- perchè??? Se un sat muore e uno nasce, la lunghezza è cmq uguale!
                        
                        check_on = 1;
                        
                        %new satellites
                        sat_born = setdiff(sat,sat_old);
                        
                        if ~isempty(sat_born)
                            %if a new satellite is going to be the pivot, its ambiguity needs
                            %to be estimated before applying the pivot change
                            if ~isempty(find(sat_born == obj.pivot(1), 1))
                                %%if it is not the only satellite with phase
                                %if (length(sat) > 1)
                                %if the former pivot is still among satellites with phase
                                if ~isempty(find(sat == obj.pivot_old(1), 1))
                                    %set the old pivot as temporary pivot
                                    pivot_tmp = obj.pivot_old(1);
                                else
                                    %find the best candidate as temporary pivot
                                    sat_tmp = setdiff(sat,obj.pivot);
                                    [max_elR, i] = max(obj.satCoordR.el(sat_tmp,1)); %#ok<ASGLU>
                                    pivot_tmp = sat_tmp(i);
                                    %reset the ambiguities of other satellites according to the temporary pivot
                                    sat_born = sat;
                                    obj.pivot_old(1,1:nRec) = pivot_tmp;
                                end
                                sat_born = setdiff(sat_born,pivot_tmp);
                                sat_slip1 = [];
                                sat_slip2 = [];
                                sat_slip = [];
                                for i=1:nRec
                                    %                                     if (length(phase) == 2)
                                    %                                         [N1_slip, N1_born] = ambiguity_init(XR0, XS, pr1_R(sat_pr), pr1_M(sat_pr), ph1_R(sat_pr), ph1_M(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip1, sat_born, distR(sat_pr), distM(sat_pr), err_tropo_R, err_tropo_M, err_iono_R, err_iono_M, pivot_tmp, phase, X_t1_t(o3+sat_pr), Cee(o3+sat_pr, o3+sat_pr)); %#ok<ASGLU>
                                    %                                         [N2_slip, N2_born] = ambiguity_init(XR0, XS, pr2_R(sat_pr), pr2_M(sat_pr), ph2_R(sat_pr), ph2_M(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip2, sat_born, distR(sat_pr), distM(sat_pr), err_tropo_R, err_tropo_M, (lambda2/lambda1)^2 * err_iono_R, (lambda2/lambda1)^2 * err_iono_M, pivot_tmp, phase, X_t1_t(o3+sat_pr), Cee(o3+sat_pr, o3+sat_pr)); %#ok<ASGLU>
                                    %                                         %[N1_slip, N1_born] = ambiguity_init(XR0, posS(:,sat_pr), pr1_Rsat(sat_pr), pr1_Msat(sat_pr), ph1_Rsat(sat_pr), ph1_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip1, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_R(sat_pr), err_tropo_M(sat_pr), err_iono_R(sat_pr), err_iono_M(sat_pr), pivot_tmp, phase, X_t1_t(o3+sat_pr)); %#ok<ASGLU>
                                    %                                         %[N2_slip, N2_born] = ambiguity_init(XR0, posS(:,sat_pr), pr2_Rsat(sat_pr), pr2_Msat(sat_pr), ph2_Rsat(sat_pr), ph2_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip2, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_R(sat_pr), err_tropo_M(sat_pr), (lambda2/lambda1)^2 * err_iono_R(sat_pr), (lambda2/lambda1)^2 * err_iono_M(sat_pr), pivot_tmp, phase, X_t1_t(o3+sat_pr)); %#ok<ASGLU>
                                    %
                                    %                                         X_t1_t(o3+sat_born,1) = N1_born;
                                    %                                         X_t1_t(o3+32+sat_born,1) = N2_born;
                                    %                                         %Cvv(o3+sat_born,o3+sat_born) = sigmaq_N1_born * eye(size(sat_born,1));
                                    %                                         %Cvv(o3+32+sat_born,o3+32+sat_born) = sigmaq_N2_born * eye(size(sat_born,1));
                                    %                                         Cvv(o3+sat_born,o3+sat_born) = sigmaq0_N * eye(size(sat_born,1));
                                    %                                         Cvv(o3+32+sat_born,o3+32+sat_born) = sigmaq0_N * eye(size(sat_born,1));
                                    %                                     else
                                    %                                         if (phase == 1)
                                    
                                    [N_slip, N_born] = ambiguity_init(XR0, X_sat(sat,:,i), pr1_R(sat,i), pr1_M(sat), ph1_R(sat,i), ph1_M(sat), snr_R(sat,i), snr_M(sat,i), obj.satCoordR.el(sat,i), obj.satCoordM.el(sat), sat, sat, sat_slip, sat_born,  obj.satCoordR.dist(sat,i), obj.satCoordM.dist(sat), err_tropo(sat,i+1), err_tropo(sat,1), err_iono(sat,i+1), err_iono(sat,1), pivot_tmp, phase, obj.Xhat_t_t(nP+goGNSS.MAX_SAT*(i-1)+sat), obj.Cee(nP+goGNSS.MAX_SAT*(i-1)+sat, nP+goGNSS.MAX_SAT*(i-1)+sat)); %#ok<ASGLU>
                                    %[N_slip, N_born] = ambiguity_init(XR0, posS(:,sat_pr), pr1_Rsat(sat_pr), pr1_Msat(sat_pr), ph1_Rsat(sat_pr), ph1_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_R(sat_pr), err_tropo_M(sat_pr), err_iono_R(sat_pr), err_iono_M(sat_pr), pivot_tmp, phase, X_t1_t(o3+sat_pr)); %#ok<ASGLU>
                                    %                                         else
                                    %                                             [N_slip, N_born] = ambiguity_init(XR0, XS, pr2_R(sat_pr), pr2_M(sat_pr), ph2_R(sat_pr), ph2_M(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip, sat_born, distR(sat_pr), distM(sat_pr), err_tropo_R, err_tropo_M, (lambda2/lambda1)^2 * err_iono_R, (lambda2/lambda1)^2 * err_iono_M, pivot_tmp, phase, X_t1_t(o3+sat_pr), Cee(o3+sat_pr, o3+sat_pr)); %#ok<ASGLU>
                                    %                                             %[N_slip, N_born] = ambiguity_init(XR0, posS(:,sat_pr), pr2_Rsat(sat_pr), pr2_Msat(sat_pr), ph2_Rsat(sat_pr), ph2_Msat(sat_pr), snr_R(sat_pr), snr_M(sat_pr), elR(sat_pr), elM(sat_pr), sat_pr, sat, sat_slip, sat_born, prRS_app(sat_pr), prMS_app(sat_pr), err_tropo_R(sat_pr), err_tropo_M(sat_pr), (lambda2/lambda1)^2 * err_iono_R(sat_pr), (lambda2/lambda1)^2 * err_iono_M(sat_pr), pivot_tmp, phase, X_t1_t(o3+sat_pr)); %#ok<ASGLU>
                                    %                                         end
                                    %
                                    obj.Xhat_t_t(nP+goGNSS.MAX_SAT*(i-1)+sat_born,1) = N_born;
                                    %Cvv(o3+sat_born,o3+sat_born) = sigmaq_N_born * eye(size(sat_born,1));
                                    Cvv(nP+goGNSS.MAX_SAT*(i-1)+sat_born,nP+goGNSS.MAX_SAT*(i-1)+sat_born) = sigmaq0_N * eye(size(sat_born,1));
                                end
                                sat_born = setdiff(sat_born,obj.pivot(1));
                                check_on = 0;
                            end
                            %end
                            
                        end
                        
                        keyboard
                        %------------------------------------------------------------------------------------
                        % PIVOT CHANGE
                        %------------------------------------------------------------------------------------
                        
                        %search for a possible PIVOT change
                        if (obj.pivot(1) ~= obj.pivot_old(1) & obj.pivot_old(1) ~= 0)
                            
                            check_pivot = 1;
                            
                            %matrix construction to update the PIVOT change
                            %sat: vector with the current visible satellites
                            %nsat: current satellites vector dimension
                            R = zeros(goGNSS.MAX_SAT);
                            R(sat,sat) = eye(length(sat));
                            R(sat,obj.pivot(1)) = -1;
                            R(obj.pivot_old,obj.pivot_old) = 0;
                            R(obj.pivot,obj.pivot) = 0;
                            
                            I0 = eye(nP);
                            Z_32_o3 = zeros(goGNSS.MAX_SAT,nP);
                            Z_o3_32 = zeros(nP,goGNSS.MAX_SAT);
                            Z_32_32 = zeros(goGNSS.MAX_SAT,goGNSS.MAX_SAT);
                            
                            %total matrix construction
                            %sat_old, sat
                            %if (length(phase) == 2)
                            %    A = [I0 Z_o3_32 Z_o3_32; Z_32_o3 R Z_32_32; Z_32_o3 Z_32_32 R];
                            %else
                                A = [I0 Z_o3_32; Z_32_o3 R];
                            %end
                            
                            %new state estimate
                            X_t1_t = A*X_t1_t;
                            
                            %re-computation of the Cee covariance matrix at the previous epoch
                            Cee = A*Cee*A';
                            
                        end
               
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
            end
            
            
            
            
            % obj........
        end
        
        
    end
    end

        
    
    % Function to fill KF loop matrices
%     methods (Access = 'public')
%         function MR_loop(obj, goObs, goIni, mode)
%             keyboard
%             
%             % obj........
%         end
%     end
    
    % Loop functions
    methods (Access = 'private')
        
    end
    
    
    
end
