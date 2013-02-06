% =========================================================================
%   OBJECT goObservation
% =========================================================================
%
% DESCRIPTION:
%   Object to read and collect observation of multiple recevers
%   (including master station)
%
% EXAMPLE:
%   goObs = goObservation();
%
% REQUIRES:
%   load_RINEX:    goGPS function to read RINEX files
%   iniReader:     class to pass the parameters of the receivers
%    - cprintf:    http://www.mathworks.com/matlabcentral/fileexchange/24093-cprintf-display-formatted-colored-text-in-the-command-window
%
% LIST of CONSTANT
%
%   idGPS          identifier of the GPS constellation inside this object
%   idGLONASS      identifier of the GLONASS constellation inside this object
%
% LIST of METHODS
%
%  GENERIC ------------------------------------------------------------
%
%   init(obj, nFreqG, nFreqR)
%   loadData(obj)
%
%   nRec = getNumRec(obj)
%   iono = getIono(obj)
%
%  REFERENCE FRAME ----------------------------------------------------
%   
%   att = getInitialAttitude(obj)
%
%  CONSTELLATION SPECIFIC ---------------------------------------------
%
%   nFreq = getGNSSnFreq(obj, idGNSS)
%
%   isAvailable = getGNSSstatus(obj, idGNSS)
%   ids = getGNSSlist(obj)
%
%  RECEIVER SPECIFIC --------------------------------------------------
%
%   setReceiverStatus(idRec, status)
%   status = getReceiverStatus(obj, idRec)
%
%   time = getTime_Ref(obj)
%   time = getTime_R(obj, idRec)
%   time = getTime_M(obj)
%
%   [X0R flagPos] = getX0_R(obj, idRec)
%   [X0M flagPos] = getX0_M(obj)
%   [XM flagPos] = getPos_M(obj, <idObs>)
%
%   setPos_M(obj, XM);
%
%   sr = getSamplingRate_R(obj, idRec)
%   sr = getSamplingRate_M(obj)
%
%  CLOCK SPECIFIC -----------------------------------------------------
%
%   initClockError_R(obj, idRec)
%   initClockError_M(obj, idRec)
%
%   dtMdot = getClockDrift_M(obj)
%   dtM    = getClockError_M(obj)
%
%  RECEIVERS GNSS SPECIFIC --------------------------------------------
%
%   eph = getGNSSeph(obj, idGNSS)
%
%   pr = getGNSSpr_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
%   pr = getGNSSpr_M(obj, idGNSS, idSat, idObs, nFreq)
%
%   ph = getGNSSph_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
%   ph = getGNSSph_M(obj, idGNSS, idSat, idObs, nFreq)
%
%   dop = getGNSSdop_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
%   dop = getGNSSdop_M(obj, idGNSS, idSat, idObs, nFreq)
%
%   snr = getGNSSsnr_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
%   snr = getGNSSsnr_M(obj, idGNSS, idSat, idObs, nFreq)
%
%   flag = getGNSSflag_R(obj, idGNSS, idRec, <idObs>)
%   flag = getGNSSflag_M(obj, idGNSS, <idObs>)
%
%  SATELLITE SPECIFIC -------------------------------------------------
%
%   isSP3 = isSP3(obj);
%   ephSP3 = X0_SP3(obj); % <deprecate> full struct containing (.time .coord .clock)
%   time  = getGNSS_SP3time(obj, <idObs>);
%   coord = getGNSS_SP3coordinates(obj <idObs>);
%   clock = getGNSS_SP3clock(obj, <idObs>);
%
%   satObs = getSatObservation(obj, idGNSS, idSat, <idObs>) % return a structure containing (.time .X .V)
%
%
%----------------------------------------------------------------------------------------------
% Copyright (C) 2009-2013 Mirko Reguzzoni, Eugenio Realini
%----------------------------------------------------------------------------------------------
%
%    Code contributed by Andrea Gatti
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

classdef goObservation < handle
    
    properties (GetAccess = 'public', SetAccess = 'private')
        % Constellations id => naming inside the goObservation object
        idGPS = 1;
        idGLONASS = 2;
        % idGALILEO = 3;
        % ...
        
        nSat = 32;      % max number of active satellites in a constellation

        idM = 1;            % id of the Master inside the various structures
    end
    
    properties (GetAccess = 'private', SetAccess = 'private')
        
    % =========================================================================
    %   General parameters
    % =========================================================================

        nObs = 0;       % number of observations avaiable
        
        iono = [];      % Ionosphere model observations

    % =========================================================================
    %   Configuration of the antennas
    % =========================================================================
        
        nRec = 1;       % number of receivers

        receiverOk = [] % for each receiver it contains the actual (enable/disable) status
        
        
        antennasRF;     % antennas Reference Frame it is a structure
                        %  .refRec => number of the remote that define the center of the RF
                        %  .pos    => [3, nRec] position of a remote antenna in the RF
                        
        attitude;       % structure containing roll, pitch, yawn                        
        
    % =========================================================================
    %   Used files
    % =========================================================================
    
        obsFile;        % array of nRec elements with field .name
        navFile;        % structure elements with field 
                        %  .name and 
                        %  .isSP3
        
        
    % =========================================================================
    %   Observations types
    % =========================================================================

        % Select constellation is logical an array that specify the available
        % GNSS constellation
        %   selectConstellation(1)  => GPS
        %   selectConstellation(2)  => GLONASS
        %   ecc...
        selectConstellation = logical([1 0]');
        
        % Sampling rate of the time table 
        % (e.g. the Master station is at 1Hz, Remotes at 10Hz => objFreq = 10)
        objSamplingRate = 1;
        samplingRate = [];      % sampling rate of each GNSS receiver  [nRec+1]
        
        navSP3;                 % structure containing precise navigation ephemerides
        
    % =========================================================================
    %   Receivers Observation (Master + Remotes)
    % =========================================================================

        % Time Chart 
        %   1st column    => time ref
        %   2nd column    => time GPS of the Master
        %   other columns => time GPS of the Remotes        
        timeChart = [];  % time table                        [nObs, (nRec+2)]        
        dateChart = [];  % date table                        [nObs, 6, (nRec+2)]        
        week = [];        % gps week of the reference time /date
        
        X0  = [];        % approximate initial position      [3, nRec+1]

        dt
        dtDot = [];      % clock drift of the receivers      [nRec+1]
        
        %   GPS SPECIFIC
        % =================================================================

        % Pseudo range and phase observations are stored in a new structure
        % they are 3D matrixes wheer each index is:
        % ? 1. index of the satellite (e.g. 1..32)
        % ? 2. index of the observation
        % ? 3. index of the frequencies (e.g. 1, 2)
        %
        % G =>  G is the RINEX id for GPS
        % R =>  R is the RINEX id for GLONASS   
        nFreqG = 2;        
        prxG  = []; % pseudo-range                          [nSat*(nRec+1), nObs, nFreqG]
        phxG  = []; % phase                                 [nSat*(nRec+1), nObs, nFreqG]
        dopxG = []; % doppler                               [nSat*(nRec+1), nObs, nFreqG]
        snrxG = []; % signal to noise ratio                 [nSat*(nRec+1), nObs, nFreqG]
        flagG = []; % it contains 0 => everything is ok     [nRec+1, nObs]
                    %             1 => missing observation
                    %             2 => interpolated observation
                
        %   GLONASS SPECIFIC
        % =================================================================
        nFreqR = 2;
        prxR  = []; % pseudo-range                          [nSat*(nRec+1), nObs, nFreqR]
        phxR  = []; % phase                                 [nSat*(nRec+1), nObs, nFreqR]
        dopxR = []; % doppler                               [nSat*(nRec+1), nObs, nFreqR]
        snrxR = []; % signal to noise ratio                 [nSat*(nRec+1), nObs, nFreqR]
        flagR = []; % it contains 0 => everything is ok     [nRec+1, nObs]
                    %             1 => missing observation  
                    %             2 => interpolated observation
        
    % =========================================================================
    %   Satellite Observations
    % =========================================================================
                            
        %   GPS SPECIFIC
        % =================================================================
        
        ephG = [];      % Ephemerides for GPS satellites
        
        timeG_tx = [];  % Time of first transmission              [nSat, 1]
        XSG_tx = [];    % Position at first trasmission           [nSat, 1]
        VSG_tx = [];    % Velocity at first trasmission           [nSat, 1]
                
        %   GLONASS SPECIFIC
        % =================================================================
        ephR = [];      % Ephemerides for GLONASS satellites

        timeR_tx = [];  % Time of first transmission              [nSat, 1]
        XSR_tx = [];    % Position at first trasmission           [nSat, 1]
        VSR_tx = [];    % Velocity at first trasmission           [nSat, 1]        
    end
    
    methods(Access = 'public')       
        % Creator 
        % the iniReader objuct MUST contains at least these sections/keys 
        % (the values are an example:)
        % - "Receivers"
        %     |- "data_path" = "../data/ublox/"
        %     |- "nRec" = 3
        %     |- "file_names" = [ "goBulgaro_01_UBX1_rover.obs" "goBulgaro_01_UBX2_rover.obs" "goBulgaro_01_UBX3_rover.obs" ]
        % - "Receiver configuration"
        %     |- "ref_rec" = 2
        %     |- "XYZ_rec1" = -3.95   0.62  0.12
        %     |- "XYZ_rec2" =  0      0     0
        %     |- "XYZ_rec3" = -3.95  -0.62  0.12
        %     |- "XYZ_ev_point" = -1.97     0     0.745
        % - "Master"
        %     |- "data_path" = "../data/permanent station/"
        %     |- "file_name" = "COMO1370.11O"
        % - "Navigational"
        %     |- "isSP3" = 0
        %     |- "data_path" = "../data/permanent station/"
        %     |- "file_name" = "COMO1370.11N"
        function obj = goObservation()            
            
        end
        
    % =========================================================================
    %  GENERIC
    % =========================================================================
        
        % Initialize the object
        function err = init(obj, ini, nFreqG, nFreqR)
            if (nargin < 3)      % nFreq hasn't been set up
                nFreqG = 2;      % Let's suppose we work in double frequencies (L1, L2, ...)
            end
            if (nargin < 4)      % nFreq hasn't been set up
                nFreqR = 2;      % Let's suppose we work in double frequencies (L1, L2, ...)
            end
            
            % Somewhere I need to read which GNSS constellation is
            % available, up to know GLONASS is not managed, so I disable it
            % in the observation object;
            obj.setGNSS(obj.idGPS, 1);
            obj.setGNSS(obj.idGLONASS, 0);

            % Set foundamental parameters
            if (obj.getGNSSstatus(obj.idGPS) == 1)
                obj.setGNSSnFreq(obj.idGPS,nFreqG);
            end
            if (obj.getGNSSstatus(obj.idGLONASS) == 1)
                obj.setGNSSnFreq(obj.idGLONASS, nFreqR);
            end    
            
            % Extract useful info from the ini file
            % ini file - easily modifiable
            % Notice that all the check of the ini parameters should be
            % done in the interface to prepare it. 
            % Up to now, no intarface is available, every error will show a
            % msgbox window and will close goGPS
            err = obj.importIniPar(ini);
            if (err.val > 0)
                msgbox(err.msg);
            end
            err = err.val;
        end
        
        % Reading files and importing observations
        function loadData(obj)
            % Loading RINEX data (Master + Receiver)            
            obj.readRINEXs(obj.obsFile, obj.navFile);
            
            % Eventually load precise SP3 ephemerides (when available)
            if obj.haveSP3
                obj.loadSP3(obj.navFile);
            end
            
            % Remove the observations of the satellites that are without ephemerides
            obj.cleanNoEphSat();
        end
        
        % Check file existence and parameters
        function inputOk = testInput(obj)
            
        end
                
        % Getter of the num of Receivers available
        function nRec = getNumRec(obj)
            nRec = obj.nRec;
        end
        
        % Get the iono observations
        function iono = getIono(obj)
            iono = obj.iono;
        end

    % =========================================================================
    %  REFERENCE FRAME 
    % =========================================================================
    
        %   Return the attitude of the structure at time t0
        function att = getInitialAttitude(obj)
            att = obj.attitude;
        end        
    
    % =========================================================================
    %  CONSTELLATION SPECIFIC
    % =========================================================================

        % Getter of the num of Frequencies available
        function nFreqR = getGNSSnFreq(obj, idGNSS)
            switch idGNSS
                case obj.idGPS
                    nFreqR = obj.nFreqR;
                case obj.idGLONASS
                    nFreqR = obj.nFreqG;
            end
        end
        
        % Get available constellations
        function isAvailable = getGNSSstatus(obj,idGNSS)
            isAvailable = (obj.selectConstellation(idGNSS) == 1);
        end
        
        % Get list of available constellations (ids)
        function ids = getGNSSlist(obj)
            ids = find(obj.selectConstellation(:) == 1);
        end
        
    % =========================================================================
    %  RECEIVER SPECIFIC
    % =========================================================================
    
        % Set the status of the receiver
        function setReceiverStatus(obj, idRec, status)
            obj.receiverOk(idRec+1) = logical(status);
        end
        
        % Get the status of the receiver
        function status = getReceiverStatus(obj, idRec)
            if (nargin == 1)              
                idRec = 0:obj.getNumRec();
            end
            idRec = idRec + 1;
            status = obj.receiverOk(idRec);
        end
        
        % Get reference time
        function time = getTime_Ref(obj)
            time = obj.timeChart(:,1);
        end
        
        % Get remote time
        function time = getTime_R(obj, idRec)
            if (nargin == 1)
                idRec = 0;
            end
            if (idRec == 0)
                idRec = 1:obj.nRec;
            end
            time = obj.timeChart(:,idRec+2); % the receiver having index 1 is the Referece and 2 is the Master
        end
        
        % Get master time
        function time = getTime_M(obj)
            time = obj.timeChart(:,2); % the receiver having index 1 is the Reference 2 is the master
        end
        
        % Get receiver approximate initial position
        % idRec => id of the receiver (if 0 get all receiver)
        function [X0R flagPos] = getX0_R(obj, idRec)
            if (nargin == 1)
                idRec = 0;
            end
            if (idRec == 0)
                idRec = 1:obj.nRec;
            end
            idRec = idRec + 1; % the receiver having index 1 is the Master
            X0R = obj.X0(:, idRec);
            
            flagPos = (sum(X0R,1) ~= 0);
        end
        
        % Get master initial position
        function [X0M flagPos] = getX0_M(obj)
            [X0M flagPos] = getPos_M(obj, 1);
        end
        
        % Get master position
        function [XM flagPos] = getPos_M(obj, idObs)
            if (nargin < 1) % if not specified set the entire position array to the value of XM 
                idObs = 1; % it should be 1 or nObs
            end
            if (idObs == 0)
                idObs = 1:obj.nObs;
            end
            % if the Master position is fixed
            if (size(obj.X0,2) == 1)
                XM = repmat(obj.X0(:),1,length(idObs));
            else
                XM = obj.X0(:, idObs);
            end
            
            flagPos = (sum(XM,1) ~= 0);
        end
        
        % Set receiver approximate position
        function setPos_M(obj, XM, idObs)
            if (nargin < 2) % if not specified set the entire position array to the value of XM 
                idObs = 1:size(obj.X0,2); % it should be 1 or nObs
            end
            if (idObs == 0)
                idObs = 1:obj.nObs;
            end
            obj.X0(1,idObs) = XM(1,:); % Set X
            obj.X0(2,idObs) = XM(2,:); % Set Y
            obj.X0(3,idObs) = XM(3,:); % Set Z            
        end
            
        % Get receiver sampling rate
        % idRec => id of the receiver (if 0 get all receiver)
        function sr = getSamplingRate_R(obj, idRec)
            if (nargin == 1)
                idRec = 0;
            end
            if (idRec == 0)
                idRec = 1:obj.nRec; % the receiver having index 1 is the Master
            end
            idRec = idRec + 1;
            sr = obj.samplingRate(idRec);
        end
        function sr = getSamplingRate_M(obj)
            sr = obj.samplingRate(1);
        end
        
    % =========================================================================
    %  CLOCK SPECIFIC 
    % =========================================================================
       
        % Compute receiver clock error and drift
        function initClockError_R(obj, idRec)
            fprintf('Computing receiver %d/%d clock error and drift...\n', idRec, obj.getNumRec());
            initClockError(obj, idRec);
            fprintf('Clock error and drift estimated\n');
        end
        
        % Compute master station clock error and drift
        function initClockError_M(obj)
            fprintf('Computing master station clock error and drift (needed to compute Doppler shift)...\n');
            initClockError(obj, 0);
            fprintf('Clock error and drift estimated\n');
        end
        
        % Compute clock error and drift (idRec == 0 is the Master)        
        function initClockError(obj, idRec)
            if (idRec == 0) % Master                
                [dt_tmp, dtDot_tmp] ...
                 = clock_error(obj.getX0_M(), obj.getTime_Ref(), ...
                               obj.getGNSSpr_M(obj.idGPS, 0, 0, 1), ...
                               obj.getGNSSsnr_M(obj.idGPS, 0, 0, 1), ...      
                               obj.getGNSSeph(obj.idGPS), ...
                               obj.getGNSS_SP3time(), ...
                               obj.getGNSS_SP3coordinates(), ...
                               obj.getGNSS_SP3clock(), ...
                               obj.getIono());
            else            % Rover
                [dt_tmp, dtDot_tmp] ...
                 = clock_error(obj.getX0_R(idRec), obj.getTime_Ref(), ...
                               obj.getGNSSpr_R(obj.idGPS, 0, idRec, 0, 1), ...
                               obj.getGNSSsnr_M(obj.idGPS, 0, idRec, 0, 1), ...
                               obj.getGNSSeph(obj.idGPS), ...
                               obj.getGNSS_SP3time(), ...
                               obj.getGNSS_SP3coordinates(), ...
                               obj.getGNSS_SP3clock(), ...
                               obj.getIono());
            end
            obj.dt(idRec+1,:) = dt_tmp';
            obj.dtDot(idRec+1,:) = dtDot_tmp';
        end
    
        % Get the drift of the clock of the Master receiver        
        function dtMdot = getClockDrift_M(obj, idObs)
            if (nargin == 1)
                idObs = 1;
            end
            if (idObs == 0)
                idObs = 1:obj.nObs;
            end
            dtMdot = obj.dtDot(1,idObs)';
        end
        
        % Get the error of the clock of the Master receiver
        function dtM = getClockError_M(obj, idObs)
            if (nargin == 1)
                idObs = 1;
            end
            if (idObs == 0)
                idObs = 1:obj.nObs;
            end
            dtM = obj.dt(1,idObs)';
        end        
        
    % =========================================================================
    %  RECEIVERS GNSS SPECIFIC 
    % =========================================================================
    
        % Get epemerides for GPS satellites
        function eph = getGNSSeph(obj, idGNSS)
            switch idGNSS
                case obj.idGPS
                    eph = obj.ephG;
                case obj.idGLONASS
                    eph = obj.ephR;
            end
        end
        
        % Get GPS pseudo-range observationa
        % idSat => id of the satellite      (if 0 get all the satellites)
        % idRec => id of the receiver       (if 0 get all receiver) 
        % idObs => id of the observation    (if 0 get all the available epocs)
        % nFreq => id of the frequency used (e.g. 1 = L1, 2 = L2, ...future frequencies...)
        function pr = getGNSSpr_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
            switch idGNSS
                case obj.idGPS
                    pr = goObservation.extractData(obj.prxG, idSat, idRec, idObs, nFreq, obj.nRec);
                case obj.idGLONASS
                    pr = goObservation.extractData(obj.prxR, idSat, idRec, idObs, nFreq, obj.nRec);
            end
        end
        function pr = getGNSSpr_M(obj, idGNSS, idSat, idObs, nFreq)
            switch idGNSS
                case obj.idGPS
                    pr = goObservation.extractData(obj.prxG, idSat, -1, idObs, nFreq, obj.nRec);
                case obj.idGLONASS
                    pr = goObservation.extractData(obj.prxR, idSat, -1, idObs, nFreq, obj.nRec);
            end
        end
        
        % Get GPS phase observationa
        % idSat => id of the satellite      (if 0 get all the satellites)
        % idRec => id of the receiver       (if 0 get all receiver) 
        % idObs => id of the observation    (if 0 get all the available epocs)
        % nFreq => id of the frequency used (e.g. 1 = L1, 2 = L2, ...future frequencies...)
        %                                   (if 0 get all the frequencies)
        function ph = getGNSSph_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
            switch idGNSS
                case obj.idGPS
                    ph = goObservation.extractData(obj.phxG, idSat, idRec, idObs, nFreq, obj.nRec);
                case obj.idGLONASS
                    ph = goObservation.extractData(obj.phxR, idSat, idRec, idObs, nFreq, obj.nRec);
            end
        end
        function ph = getGNSSph_M(obj, idGNSS, idSat, idObs, nFreq)
            switch idGNSS
                case obj.idGPS
                    ph = goObservation.extractData(obj.phxG, idSat, -1, idObs, nFreq, obj.nRec);
                case obj.idGLONASS
                    ph = goObservation.extractData(obj.phxR, idSat, -1, idObs, nFreq, obj.nRec);
            end
        end
        
        % Get GPS doppler observationa
        % idSat => id of the satellite      (if 0 get all the satellites)
        % idRec => id of the receiver       (if 0 get all receiver)
        % idObs => id of the observation    (if 0 get all the available epocs)
        % nFreq => id of the frequency used (e.g. 1 = L1, 2 = L2, ...future frequencies...)
        %                                   (if 0 get all the frequencies)
        function dop = getGNSSdop_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
            switch idGNSS
                case obj.idGPS
                    dop = goObservation.extractData(obj.dopxG, idSat, idRec, idObs, nFreq, obj.nRec);
                case obj.idGLONASS
                    dop = goObservation.extractData(obj.dopxR, idSat, idRec, idObs, nFreq, obj.nRec);
            end
        end
        function dop = getGNSSdop_M(obj, idGNSS, idSat, idObs, nFreq)
            switch idGNSS
                case obj.idGPS
                    dop = goObservation.extractData(obj.dopxG, idSat, -1, idObs, nFreq, obj.nRec);
                case obj.idGLONASS
                    dop = goObservation.extractData(obj.dopxR, idSat, -1, idObs, nFreq, obj.nRec);
            end
        end
        
        % Get GPS signal to noise ratio observationa
        % idSat => id of the satellite      (if 0 get all the satellites)
        % idRec => id of the receiver       (if 0 get all receiver)
        % idObs => id of the observation    (if 0 get all the available epocs)
        % nFreq => id of the frequency used (e.g. 1 = L1, 2 = L2, ...future frequencies...)
        %                                   (if 0 get all the frequencies)
        function snr = getGNSSsnr_R(obj, idGNSS, idSat, idRec, idObs, nFreq)
            switch idGNSS
                case obj.idGPS
                    snr = goObservation.extractData(obj.dopxG, idSat, idRec, idObs, nFreq, obj.nRec);
                case obj.idGLONASS
                    snr = goObservation.extractData(obj.dopxR, idSat, idRec, idObs, nFreq, obj.nRec);
            end
        end
        function snr = getGNSSsnr_M(obj, idGNSS, idSat, idObs, nFreq)
            switch idGNSS
                case obj.idGPS
                    snr = goObservation.extractData(obj.dopxG, idSat, -1, idObs, nFreq, obj.nRec);
                case obj.idGLONASS
                    snr = goObservation.extractData(obj.dopxR, idSat, -1, idObs, nFreq, obj.nRec);
            end
        end
        
        % Get receiver flag array
        % idRec => id of the receiver (if 0 get all receiver)
        function flag = getGNSSflag_R(obj, idGNSS, idRec, idObs)
            if (nargin == 2)
                idObs = 0;
            end
            
            if (idRec == 0)
                idRec = 1:obj.nRec; % the receiver having index 1 is the Master
            end
            idRec = idRec + 1;
            
            if (idObs == 0)
                idObs = 1:size(obj.flagG,2);
            end
            
            switch idGNSS
                case obj.idGPS
                    flag = obj.flagG(idRec,idObs);
                case obj.idGLONASS
                    flag = obj.flagR(idRec,idObs);
            end
        end
        function flag = getGNSSflag_M(obj, idGNSS, idObs)
            if (nargin == 2)
                idObs = 0;
            end
            
            if (idObs == 0)
                idObs = 1:size(obj.flagG,2);
            end
            
            switch idGNSS
                case obj.idGPS
                    flag = obj.flagG(1,idObs);
                case obj.idGLONASS
                    flag = obj.flagR(1,idObs);
            end
        end
        
     % ======================================================================
     %  SATELLITE SPECIFIC
     % ======================================================================
     
        % Return if the navigation file is SP3 (precise orbits) or not
        function haveSP3 = haveSP3(obj)
            haveSP3 = obj.navFile.isSP3;
        end
        
        % Init and read SP3 file
        function loadSP3(obj, navFile)
            if (navFile.isSP3())
                fprintf('Reading SP3 file...\n');
                
                [obj.navSP3.time, obj.navSP3.pos, obj.navSP3.clock] = load_SP3(navFile.name, obj.timeChart(:,1), obj.week);
                obj.navSP3.nObs = size(obj.navSP3.time,1);
            end
        end
        
        % get the full structure of the SP3 ephemerides => .time, .pos, .clock, .nObs
        function ephSP3 = X0_SP3(obj)
           ephSP3 = obj.navSP3 ;
        end
        
        % extract the time information from SP3 ephemerides
        function time  = getGNSS_SP3time(obj, idObs)
            if (obj.haveSP3())
                if (nargin == 1)
                    idObs = 0;
                end
                if (idObs == 0)
                    idObs = 1:obj.navSP3.nObs;
                end
                time = obj.navSP3.time(idObs,1);
            else
                time = [];
            end
        end
        
        % extract the position information from SP3 ephemerides
        function pos = getGNSS_SP3coordinates(obj, idObs)
            if (obj.haveSP3())
                if (nargin == 1)
                    idObs = 0;
                end
                if (idObs == 0)
                    idObs = 1:obj.navSP3.nObs;
                end
                pos = obj.navSP3.pos(:,:,idObs);
            else
                pos = [];
            end
        end
        
        % extract the time information from SP3 ephemerides
        function clock = getGNSS_SP3clock(obj, idObs)
            if (obj.haveSP3())
                if (nargin == 1)
                    idObs = 0;
                end
                if (idObs == 0)
                    idObs = 1:obj.navSP3.nObs;
                end
                clock = obj.navSP3.clock(:,idObs);
            else
                clock = [];
            end
        end
        
%         function satObs = getSatObservation(obj, idGNSS, idSat, idObs)
%         end        
    end
    % =========================================================================
    %    PRIVATE METHODS
    % =========================================================================

    methods(Access = 'private')
        % Allocate the memory necessary to store all the observations of
        % the receivers (master station included)
        function allocateMemory(obj, nObs, nRec, nFreqG, nFreqR)
            obj.timeChart = zeros(nObs, nRec+2);                % time table
            obj.dateChart = zeros(nObs,6,nRec+2);               % date
            
            obj.X0  = zeros(3, nRec+1);                  % approximate initial position

            obj.prxG  = zeros(obj.nSat*(nRec+1), nObs, nFreqG); % pseudo-range
            obj.phxG  = zeros(obj.nSat*(nRec+1), nObs, nFreqG); % phase
            obj.dopxG = zeros(obj.nSat*(nRec+1), nObs, nFreqG); % doppler
            obj.snrxG = zeros(obj.nSat*(nRec+1), nObs, nFreqG); % signal to noise ratio
            obj.flagG = zeros(nRec+1, nObs);                    % flags
            
            obj.prxR  = zeros(obj.nSat*(nRec+1), nObs, nFreqR); % pseudo-range
            obj.phxR  = zeros(obj.nSat*(nRec+1), nObs, nFreqR); % phase
            obj.dopxR = zeros(obj.nSat*(nRec+1), nObs, nFreqR); % doppler
            obj.snrxR = zeros(obj.nSat*(nRec+1), nObs, 1);      % signal to noise ratio
            obj.flagR = zeros(nRec+1, nObs);                    % flags
            
            % Up to now, These're not used for the receivers but the memory is allocated anyway for compatibility reason
            obj.dt = zeros(nRec+1,nObs);
            obj.dtDot = zeros(nRec+1,nObs);
        end
        
        % Extract useful info from the ini file
        % and save them in the object
        function err = importIniPar(obj,ini)
            err.val = 0;    % Everything is ok 
            err.msg = '';
            
            % Save the number of receivers
            nR = ini.getData('Receivers','nRec');
            if (isempty(nR))
                err.val = 2;
                err.msg = 'The receiver number has not been specified. ';                
            end
            
            obj.receiverOk = false(nR+1,1);

            % Receivers file
            data_path = ini.getData('Receivers','data_path');             
            if (isempty(data_path)) 
                err.val = 1;
                err.msg = 'The receiver data path is missing.';
                return
            end
            file_names = ini.getData('Receivers','file_name');
            if (isempty(file_names)) 
                err.val = 1;
                err.msg = 'The receivers file names is missing.';
                return
            else
                % If missing, get the number of receiver from the file name list
                if isempty(nR)
                    err.val = 0;
                    err.msg ='';
                    nR = length(file_names);
                end
                
                % Assign the filenames
                for r = nR:-1:1                
                    obj.obsFile(r+1).name = [data_path file_names{r}];
                    if isempty(dir(obj.obsFile(r+1).name))  % Check if the receiver file is available
                        obj.setReceiverStatus(r, 0);
                        err.val = 1;
                        err.msg = 'The receivers file names is missing.';
                        nR = nR - 1;
                    end
                end
            end
            if (nR > 0)
                obj.setNumRec(nR);
            else
                err.val = 1;
                err.msg = 'No rover observations are available, processing interrupted';
                return
            end

            % Master file
            data_path = ini.getData('Master','data_path');
            if (isempty(data_path)) 
                err.val = 1;
                err.msg = 'The master data path is missing.';
                return
            end            
            tmp = ini.getData('Master','file_name');
            obj.obsFile(obj.idM).name = [data_path tmp];            
            if (isempty(tmp) || isempty(dir(obj.obsFile(obj.idM).name))) 
                err.val = 1;
                err.msg = 'The master file name is missing.';
                return
            end
            obj.setReceiverStatus(0,isempty(dir(obj.obsFile(obj.idM).name)));
            
            % Navigation file
            data_path = ini.getData('Navigational','data_path');
            if isempty(data_path) 
                err.val = 1;
                err.msg = 'The navigational data path is missing.';
                return
            end            
            tmp = ini.getData('Navigational','file_name');
            obj.navFile.name = [data_path tmp];
            if (isempty(tmp) || isempty(dir(obj.navFile.name)))
                err.val = 1;
                err.msg = 'The navigational file is missing.';
                return
            end            
            obj.navFile.isSP3 = ini.getData('Navigational','isSP3');
            if isempty(obj.navFile.isSP3)
                obj.navFile.isSP3 = 0;
                err.val = 2;
                err.msg = [err.msg 'The SP3 flag has not been set. '];
            end
            
            % Store receiver positions
            obj.antennasRF.ref = ini.getData('Antennas RF','ref');
            if (isempty(obj.antennasRF.ref))
                err.val = 2;
                err.msg = [err.msg 'The reference antenna has not been set. '];
                obj.antennasRF.ref = 1; % the first receiver is the reference
            end 
            
            obj.antennasRF.pos = zeros(3,nR+1);
            for r = 1:nR
                tmp = ini.getData('Antennas RF',['XYZ_ant' num2str(r)]);
                if (isempty(tmp)) 
                    err.val = 2;
                    err.msg = [err.msg 'The antenna ' num2str(r) ' position has not been set. '];
                    obj.antennasRF.pos(:,r+1) = zeros(3,1);
                else
                    obj.antennasRF.pos(:,r+1) = tmp;
                end
            end        
            
            % Store structure attitude at time T0
            % must be read from INI (ToDo)
            obj.attitude.roll = 0;
            obj.attitude.pitch = 0;
            obj.attitude.yaw = 0;

            % Computation point
            tmp = ini.getData('Antennas RF','XYZ_ev_point');            
            if (isempty(tmp))
                err.val = 2;
                err.msg = [err.msg 'The computation point position has not been set. '];
                obj.antennasRF.pos(:,obj.idM) = zeros(3,1);
            else
                obj.antennasRF.pos(:,obj.idM) = tmp;
            end
        end
        
        % Read RINEX files and fill the object
        function readRINEXs(obj, obsFile, navFile)
            nR = obj.getNumRec()+1; % The first position is for the MASTER            
            
            % GPS
            if (obj.getGNSSstatus(obj.idGPS) == 1)
                nFreqG = obj.getGNSSnFreq(obj.idGPS);

                % Reading Master and Remote observations
                % We use cells to read the receiver files
                prGcell = cell(nR,nFreqG);
                phGcell = cell(nR,nFreqG);
                dopGcell = cell(nR,nFreqG);
                snrGcell = cell(nR,nFreqG);
            else
                nFreqG = 0;
            end
            
            % GLONASS
            if (obj.getGNSSstatus(obj.idGLONASS) == 1)
                nFreqR = obj.getGLONASSnFreq(obj.idGLONASS);
                
                % Reading Master and Remote observations
                % We use cells to read the receiver files
                prRcell = cell(nR,nFreqR);
                phRcell = cell(nR,nFreqR);
                dopRcell = cell(nR,nFreqR);
                snrRcell = cell(nR,1);
            else
                nFreqR = 0;
            end
            
            timeRoundedCell = cell(nR,1);
            timeCell = cell(nR,1);
            dateCell = cell(nR,1);
            posCell = cell(nR,1);
            
            % Reading MASTER observations ------------------------------------------------
            fprintf('Reading RINEX of the master station...\n');
            t0 = tic();
            % It should be better to change load rinex in such a way to
            % already read just the needed observations...
            % I read both GLONASS and GPS            
            if (obj.getGNSSstatus(obj.idGPS) == 1) && (obj.getGNSSstatus(obj.idGLONASS) == 1)
                [prGcell{obj.idM,1}, ~, phGcell{obj.idM,1}, ~, prGcell{obj.idM,2}, ~, phGcell{obj.idM,2}, ~, ...
                    dopGcell{obj.idM,1}, ~, dopGcell{obj.idM,2}, ~, snrGcell{obj.idM,1}, ~, snrGcell{obj.idM,2}, ~, ...
                    prRcell{obj.idM,1}, ~, phRcell{obj.idM,1}, ~, prRcell{obj.idM,2}, ~, phRcell{obj.idM,2}, ~, ...
                    dopRcell{obj.idM,1}, ~, dopRcell{obj.idM,2}, ~, snrRcell{obj.idM,1}, ~, ...
                    timeRoundedCell{obj.idM}, timeCell{obj.idM}, ~, dateCell{obj.idM}, posCell{obj.idM}, ~, obj.ephG, obj.iono, obj.ephR] = ...
                    load_RINEX(navFile.isSP3, obsFile(obj.idM).name, navFile.name);
                % I read only GPS
            elseif (obj.getGNSSstatus(obj.idGPS) == 1)
                [prGcell{obj.idM,1}, ~, phGcell{obj.idM,1}, ~, prGcell{obj.idM,2}, ~, phGcell{obj.idM,2}, ~, ...
                    dopGcell{obj.idM,1}, ~, dopGcell{obj.idM,2}, ~, snrGcell{obj.idM,1}, ~, snrGcell{obj.idM,2}, ~, ...
                    ~, ~, ~, ~, ~, ~, ~, ~, ...
                    ~, ~, ~, ~, ~, ~, ...
                    timeRoundedCell{obj.idM}, timeCell{obj.idM}, ~, dateCell{obj.idM}, posCell{obj.idM}, ~, obj.ephG, obj.iono, ~] = ...
                    load_RINEX(navFile.isSP3, obsFile(obj.idM).name, navFile.name);
                % I read only GLONASS
            elseif (obj.getGNSSstatus(obj.idGLONASS) == 1)
                [~, ~, ~, ~, ~, ~, ~, ~, ...
                    ~, ~, ~, ~, ~, ~, ~, ~, ...
                    prRcell{obj.idM,1}, ~, phRcell{obj.idM,1}, ~, prRcell{obj.idM,2}, ~, phRcell{obj.idM,2}, ~, ...
                    dopRcell{obj.idM,1}, ~, dopRcell{obj.idM,2}, ~, snrRcell{obj.idM,1}, ~, ...
                    timeRoundedCell{obj.idM}, timeCell{obj.idM}, ~, dateCell{obj.idM}, posCell{obj.idM}, ~, ~, obj.iono, obj.ephR] = ...
                    load_RINEX(navFile.isSP3, obsFile(obj.idM).name, navFile.name);
            end
            
            t0 = toc(t0);
            fprintf(['  Data read in ' num2str(t0) ' seconds\n']);
                
            tmp = timeRoundedCell{obj.idM};  % it's not possible to access single array values when they are inside a cell container
            obj.samplingRate(obj.idM) =  1/median(tmp(2:end) - tmp(1:end-1));    % Save master sample rate

            interval = zeros(nR-1,1);
            % Reading receivers ----------------------------------------------------------
            for r=2:nR
                fprintf('Reading RINEX of the receiver %d/%d...\n', r-1, nR-1);
                t0 = tic();

                % I read both GLONASS and GPS
                if (obj.getGNSSstatus(obj.idGPS) == 1) && (obj.getGNSSstatus(obj.idGLONASS) == 1)
                    [prGcell{r,1}, ~, phGcell{r,1}, ~, prGcell{r,2}, ~, phGcell{r,2}, ~, ...
                        dopGcell{r,1}, ~, dopGcell{r,2}, ~, snrGcell{r,1}, ~, snrGcell{r,2}, ~, ...
                        prRcell{r,1}, ~, phRcell{r,1}, ~, prRcell{r,2}, ~, phRcell{r,2}, ~, ...
                        dopRcell{r,1}, ~, dopRcell{r,2}, ~, snrRcell{r,1}, ~, ...
                        timeRoundedCell{r}, timeCell{r}, ~, dateCell{r}, posCell{r}, ~, ~, ~, ~] = ...
                        load_RINEX(navFile.isSP3, obsFile(r).name, navFile.name);
                % I read only GPS
                elseif (obj.getGNSSstatus(obj.idGPS) == 1)
                    [prGcell{r,1}, ~, phGcell{r,1}, ~, prGcell{r,2}, ~, phGcell{r,2}, ~, ...
                        dopGcell{r,1}, ~, dopGcell{r,2}, ~, snrGcell{r,1}, ~, snrGcell{r,2}, ~, ...
                        ~, ~, ~, ~, ~, ~, ~, ~, ...
                        ~, ~, ~, ~, ~, ~, ...
                        timeRoundedCell{r}, timeCell{r}, ~, dateCell{r}, posCell{r}, ~, ~, ~, ~] = ...
                        load_RINEX(navFile.isSP3, obsFile(r).name, navFile.name);
                % I read only GLONASS    
                elseif (obj.getGNSSstatus(obj.idGLONASS) == 1)
                    [~, ~, ~, ~, ~, ~, ~, ~, ...
                        ~, ~, ~, ~, ~, ~, ~, ~, ...
                        prRcell{r,1}, ~, phRcell{r,1}, ~, prRcell{r,2}, ~, phRcell{r,2}, ~, ...
                        dopRcell{r,1}, ~, dopRcell{r,2}, ~, snrRcell{r,1}, ~, ...
                        timeRoundedCell{r}, timeCell{r}, ~, dateCell{r}, posCell{r}, ~, ~, ~, ~] = ...
                        load_RINEX(navFile.isSP3, obsFile(r).name, navFile.name);
                end

                t0 = toc(t0);
               	fprintf(['  Data read in ' num2str(t0) ' seconds\n']);
                
                % Maintain first and last epoch for all the receivers (to check receivers syncronization)
                tmp = timeRoundedCell{r};  % it's not possible to access single array values when they are inside a cell container
                if exist('startTime','var')
                    % First and last epochs with at least one receiver observations
                    epoch_init = min(epoch_init,tmp(1));
                    epoch_end = max(epoch_end,tmp(end));
                    % First and last epochs with all receiver observations available
                    startTime = max(startTime,tmp(1));
                    stopTime = min(stopTime,tmp(end));
                else
                    epoch_init = tmp(1);
                    epoch_end = tmp(end);
                    startTime = tmp(1);
                    stopTime = tmp(end);
                end
                % Sampling rate of the single receiver
                % using median to avoid outliers (# good data must be > # missing data)
                interval(r-1) = median(tmp(2:end) - tmp(1:end-1));
                obj.samplingRate(r) = 1/interval(r-1);    % Save single receiver sampling rate
                clear tmp;
            end
            interval = mean(interval); % mean interval among all the receivers
            time_GPS = (epoch_init:interval:epoch_end)'; % reference GPS time for all the receivers
            
            numObs = length(time_GPS);
            obj.nObs = numObs;
            obj.allocateMemory(numObs, nR-1, nFreqG, nFreqR);
            obj.timeChart(:,1) = time_GPS;  % Set up reference time 
            
            % Cutting MASTER and REMOTEs epochs on reference time
            % Note: normally load RINEX already do this when the MASTER and
            % the REMOTE RINEX are read together (but here I need to read
            % them separately - I have more than one Receivers)
            void = {};
            if (obj.getGNSSstatus(obj.idGPS) == 1) && (obj.getGNSSstatus(obj.idGLONASS) == 1)
                [timeRoundedCell, timeCell, dateCell, prGcell, phGcell, dopGcell, snrGcell, prRcell, phRcell, dopRcell, snrRcell] = obj.recObsCutter(timeRoundedCell, timeCell, dateCell, prGcell, phGcell, dopGcell, snrGcell, prRcell, phRcell, dopRcell, snrRcell, nR, startTime, stopTime);
            elseif (obj.getGNSSstatus(obj.idGPS) == 1)
                [timeRoundedCell, timeCell, dateCell, prGcell, phGcell, dopGcell, snrGcell, ~, ~, ~, ~] = obj.recObsCutter(timeRoundedCell, timeCell, dateCell, prGcell, phGcell, dopGcell, snrGcell, void, void, void, void, nR, startTime, stopTime);
            elseif (obj.getGLONASSstatus(obj.idGPS) == 1)
                [timeRoundedCell, timeCell, dateCell, ~, ~, ~, ~, prRcell, phRcell, dopRcell, snrRcell] = obj.recObsCutter(timeRoundedCell, timeCell, dateCell, void, void, void, void, prRcell, phRcell, dopRcell, snrRcell, nR, startTime, stopTime);
            end            
            
            % Fill object tables
            for r=1:nR
                [ncu ids]= intersect(time_GPS, timeRoundedCell{r});
                obj.timeChart(ids, r+1) = timeCell{r};                                   % time
                obj.dateChart(ids, :, r+1) = dateCell{r};                                       % date
                
                if (obj.getGNSSstatus(obj.idGPS) == 1)
                    for f = 1:nFreqG()
                        obj.prxG((1:obj.nSat)+obj.nSat*(r-1), ids, f) = prGcell{r,f};     % pseudo-range
                        obj.phxG((1:obj.nSat)+obj.nSat*(r-1), ids, f) = phGcell{r,f};     % phase
                        obj.dopxG((1:obj.nSat)+obj.nSat*(r-1), ids, f) = dopGcell{r,f};   % doppler
                        obj.snrxG((1:obj.nSat)+obj.nSat*(r-1), ids, f) = snrGcell{r,f};   % signal to noise ratio
                    end
                    
                    % detect missing epochs (up to now it just checks for fully missing epoch)
                    obj.flagG(r, :) = (obj.timeChart(:, r+1) == 0)';                      % flags
                end
                
                if (obj.getGNSSstatus(obj.idGLONASS) == 1)
                    for f = 1:nFreqR()
                        obj.prxR((1:obj.nSat)+obj.nSat*(r-1), ids, f) = prRcell{r,f};     % pseudo-range
                        obj.phxR((1:obj.nSat)+obj.nSat*(r-1), ids, f) = phRcell{r,f};     % phase
                        obj.dopxR((1:obj.nSat)+obj.nSat*(r-1), ids, f) = dopRcell{r,f};   % doppler                    
                    end                
                    obj.snrxR((1:obj.nSat)+obj.nSat*(r-1), ids, f) = snrRcell{r,1};       % signal to noise ratio
    
                    % detect missing epochs (up to now it just checks for fully missing epoch)
                    obj.flagR(r, :) = (obj.timeChart(:, r+1) == 0)';                      % flags
                end                
                
                obj.X0(:, r) = posCell{r};
            end

            if (obj.getGNSSstatus(obj.idGPS) == 1) 
                clear    prGcell phGcell dopGcell snrGcell;
            end
            if (obj.getGNSSstatus(obj.idGLONASS) == 1)
                clear    prRcell phRcell dopRcell snrRcell;
            end
            clear    timeRoundedCell timeCell dateCell posCell;
            
            %GPS week number
            obj.dateChart(:,:,1) = obj.dateChart(:,:,obj.idM+1);
            obj.dateChart(:,1,1) = obj.dateChart(:,1,1) + 2000;
            obj.week = floor((datenum(obj.dateChart(:,:,1)) - datenum([1980,1,6,0,0,0]))/7);

            
            fprintf('The data is ready!\n')
        end        
        
        % Remove the observations of the satellites that are without ephemerides        
        function cleanNoEphSat(obj)
            nR = obj.getNumRec();
            % When SP3 are available, all the ephemerides for the
            % satellites are available.
            if (~obj.haveSP3())
                % remove satellites without ephemerides (GPS)
                if (obj.getGNSSstatus(obj.idGPS))
                    eph = getGNSSeph(obj, obj.idGPS);
                    rmSat = setdiff(1:obj.nSat,unique(eph(1,:))); % Check satellite idS that are not in the list of svPRN
                    for r = 1:nR
                        obj.prxG(rmSat+obj.nSat*(r-1),:,:) = 0;
                        obj.phxG(rmSat+obj.nSat*(r-1),:,:) = 0;
                        obj.dopxG(rmSat+obj.nSat*(r-1),:,:) = 0;
                        obj.snrxG(rmSat+obj.nSat*(r-1),:,:) = 0;
                    end                    
                end
                % remove satellites without ephemerides (GLONASS)
                if (obj.getGNSSstatus(obj.idGLONASS))
                    eph = getGNSSeph(obj, obj.idGLONASS);
                    rmSat = setdiff(1:obj.nSat,unique(eph(1,:)));
                    for r = 1:nR
                        obj.prxR(rmSat+obj.nSat*(r-1),:,:) = 0;
                        obj.phxR(rmSat+obj.nSat*(r-1),:,:) = 0;
                        obj.dopxR(rmSat+obj.nSat*(r-1),:,:) = 0;
                        obj.snrxR(rmSat+obj.nSat*(r-1),:,:) = 0;
                    end
                end
            end
        end
        
        % Get initial position (idRec == 0 return the master)
        % Commented => it is not currently used
        % function [X0 flagPos] = getX0(obj, idRec)
        %    if (idRec == 0)
        %        [X0 flagPos] = obj.getX0_M;
        %    else
        %        [X0 flagPos] = obj.getX0_R(idRec);
        %    end
        %end

    % =========================================================================
    %    Setter
    % =========================================================================
    
        % Setter of the num of Receiver available
        function setNumRec(obj, nRec)
            obj.nRec = nRec;
        end
        
        % Setter of the num of Frequencies to be used available
        function setGNSSnFreq(obj, idGNSS, nFreq)
            switch idGNSS
                case obj.idGPS
                    obj.nFreqG = nFreq;
                case obj.idGLONASS
                    obj.nFreqR = nFreq;
            end                    
        end
        
        % Setters of the available constellations
        function setGNSS(obj, idGNSS, isAvailable)
            if (nargin == 1)
                isAvailable = 1;
            end
            obj.selectConstellation(idGNSS) = logical(isAvailable);
        end
        
        % cut the observations in a cell array structure according to a
        % start and stop epoch
        function [timeRoundedCell, timeCell, dateCell, prGcell, phGcell, dopGcell, snrGcell, prRcell, phRcell, dopRcell, snrRcell] = recObsCutter(obj, timeRoundedCell, timeCell, dateCell, prGcell, phGcell, dopGcell, snrGcell, prRcell, phRcell, dopRcell, snrRcell, nRec, startTime, stopTime)
            % For each receiver
            for r = 1:nRec
                % Find the indexes of the first and last epoch to keep
                [ncu startId] = min(abs(timeRoundedCell{r}-startTime));
                [ncu stopId] = min(abs(timeRoundedCell{r}-stopTime));
                
                % Cut the array of the data contained in the cell
                timeRoundedCell{r} = goObservation.cutCellItem(timeRoundedCell{r}, startId, stopId);
                timeCell{r} = goObservation.cutCellItem(timeCell{r}, startId, stopId);
                dateCell{r} = goObservation.cutCellItem(dateCell{r}, startId, stopId);
                if (obj.getGNSSstatus(obj.idGPS) == 1)
                    for f = 1:obj.getGNSSnFreq(obj.idGPS)
                        prGcell{r,f} = goObservation.cutCellItem(prGcell{r,f}', startId, stopId)';
                        phGcell{r,f} = goObservation.cutCellItem(phGcell{r,f}', startId, stopId)';
                        dopGcell{r,f} = goObservation.cutCellItem(dopGcell{r,f}', startId, stopId)';
                        snrGcell{r,f} = goObservation.cutCellItem(snrGcell{r,f}', startId, stopId)';
                    end
                end
                if (obj.getGNSSstatus(obj.idGLONASS) == 1)
                    for f = 1:obj.getGNSSnFreq(obj.idGLONASS)
                        prRcell{r,f} = goObservation.cutCellItem(prRcell{r,f}', startId, stopId)';
                        phRcell{r,f} = goObservation.cutCellItem(phRcell{r,f}', startId, stopId)';
                        dopRcell{r,f} = goObservation.cutCellItem(dopRcell{r,f}', startId, stopId)';
                    end
                    snrRcell{r,1} = goObservation.cutCellItem(snrRcell{r,1}', startId, stopId)';
                end
            end
        end
    end
% =========================================================================
%    STATIC PRIVATE METHODS
% =========================================================================

    methods(Static, Access = 'private')
        
        % Extract a subTable from a generic observation table with size [nSat*(nRec+1), nObs, nFreq]
        function subtable = extractData(table, idSat, idRec, idObs, nFreq, nRec)
           nS = (size(table,1)/(nRec+1));
           if (idSat == 0)
               idSat = 1:nS;
           end
           
           if (idRec == 0)
               idRec = (1:nRec); % the receiver having index 1 is the Master
           end
           idRec = idRec + 1; % Remote receiver start from position 2 (the first is occupied by the Master)

           if (idObs == 0)
               idObs = 1:size(table,2);
           end

           if (nFreq == 0)
               nFreq = 1:size(table,3);
           end
           
           if (idRec == 0) % => I want Master observations
               firstIndex = idSat; 
           else
                % Compute the id of the lines containing the selected receivers data
                firstIndex = repmat(nS*(idRec-1),length(idSat),1)+repmat(idSat,length(idRec),1)';
           end
           
           % Extract the table
           subtable = table(firstIndex(:), idObs, nFreq);
        end
        
        % Cut an array stored in a cell of a cell array, according to start
        % and stop indexes
        function tmp = cutCellItem(cellItem, startId, stopId)
            tmp = cellItem(startId:stopId,:);
        end
    end

end
