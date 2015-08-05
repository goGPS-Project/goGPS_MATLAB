%   OBJECT goGNSS
% =========================================================================
%
% DESCRIPTION
%   Object to manage a useful functions / standard parameters of goGPS
%   Under development no methods are used up to now...
%
% EXAMPLE
%   go = goGNSS();
%
% FOR A LIST OF CONSTANTs and METHODS use doc goGNSS
%
%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%----------------------------------------------------------------------------------------------
%
%    Code contributed by Stefano Caldera, Andrea Gatti, Lisa Pertusini
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

classdef goGNSS < handle
    
    % Constant values
    % => to discriminate them from function (in autocompletion) they are
    % written in capital letters
    properties (Constant)
        % GENERIC CONSTANTS -----------------------------------------------
        
        V_LIGHT = 299792458;                  % velocity of light in the void [m/s]
        
        MAX_SAT = 32                          % Maximum number of active satellites in a constellation
        
        % CONSTELLATION REF -----------------------------------------------
        % CRS parameters, according to each GNSS system CRS definition
        % (ICD document in brackets):
        %
        % *_GPS --> WGS-84   (IS-GPS200E)
        % *_GLO --> PZ-90    (GLONASS-ICD 5.1)
        % *_GAL --> GTRF     (Galileo-ICD 1.1)
        % *_BDS --> CGCS2000 (BeiDou-ICD 1.0)
        % *_QZS --> WGS-84   (IS-QZSS 1.5D)
        
        ELL_A_GPS = 6378137;                          % GPS (WGS-84)      Ellipsoid semi-major axis [m]
        ELL_A_GLO = 6378136;                          % GLONASS (PZ-90)   Ellipsoid semi-major axis [m]
        ELL_A_GAL = 6378137;                          % Galileo (GTRF)    Ellipsoid semi-major axis [m]
        ELL_A_BDS = 6378136;                          % BeiDou (CGCS2000) Ellipsoid semi-major axis [m]
        ELL_A_QZS = 6378137;                          % QZSS (WGS-84)     Ellipsoid semi-major axis [m]
        
        ELL_F_GPS = 1/298.257222101;                  % GPS (WGS-84)      Ellipsoid flattening
        ELL_F_GLO = 1/298.257222101;                  % GLONASS (PZ-90)   Ellipsoid flattening
        ELL_F_GAL = 1/298.257222101;                  % Galileo (GTRF)    Ellipsoid flattening
        ELL_F_BDS = 1/298.257222101;                  % BeiDou (CGCS2000) Ellipsoid flattening
        ELL_F_QZS = 1/298.257222101;                  % QZSS (WGS-84)     Ellipsoid flattening
        
        ELL_E_GPS = sqrt(1-(1-goGNSS.ELL_F_GPS)^2);   % GPS (WGS-84)      Eccentricity
        ELL_E_GLO = sqrt(1-(1-goGNSS.ELL_F_GLO)^2);   % GLONASS (PZ-90)   Eccentricity
        ELL_E_GAL = sqrt(1-(1-goGNSS.ELL_F_GAL)^2);   % Galileo (GTRF)    Eccentricity
        ELL_E_BDS = sqrt(1-(1-goGNSS.ELL_F_BDS)^2);   % BeiDou (CGCS2000) Eccentricity
        ELL_E_QZS = sqrt(1-(1-goGNSS.ELL_F_QZS)^2);   % QZSS (WGS-84)     Eccentricity
        
        GM_GPS = 3.986005e14;                     % GPS     Gravitational constant * (mass of Earth) [m^3/s^2]
        GM_GLO = 3.9860044e14;                    % GLONASS Gravitational constant * (mass of Earth) [m^3/s^2]
        GM_GAL = 3.986004418e14;                  % Galileo Gravitational constant * (mass of Earth) [m^3/s^2]
        GM_BDS = 3.986004418e14;                  % BeiDou  Gravitational constant * (mass of Earth) [m^3/s^2]
        GM_QZS = 3.986005e14;                     % QZSS    Gravitational constant * (mass of Earth) [m^3/s^2]
                                                  % (NOTE: these values are not actually called from goGNSS.m
                                                  %        by ecc_anomaly.m for computation time reasons; if
                                                  %        it's needed to change them, please update also
                                                  %        the values in ecc_anomaly.m)
        
        OMEGAE_DOT_GPS = 7.2921151467e-5;             % GPS     Angular velocity of the Earth rotation [rad/s]
        OMEGAE_DOT_GLO = 7.292115e-5;                 % GLONASS Angular velocity of the Earth rotation [rad/s]
        OMEGAE_DOT_GAL = 7.2921151467e-5;             % Galileo Angular velocity of the Earth rotation [rad/s]
        OMEGAE_DOT_BDS = 7.292115e-5;                 % BeiDou  Angular velocity of the Earth rotation [rad/s]
        OMEGAE_DOT_QZS = 7.2921151467e-5;             % QZSS    Angular velocity of the Earth rotation [rad/s]
        
        J2_GLO = 1.0826257e-3;                        % GLONASS second zonal harmonic of the geopotential
        
        PI_ORBIT = 3.1415926535898;                   % pi value used for orbit computation
        CIRCLE_RAD = 2*goGNSS.PI_ORBIT;               % 2 pi (NOTE: this value is not actually called from goGNSS.m
                                                      %             for computation time reasons; if it's needed to
                                                      %             change it, please update also ecc_anomaly.m and
                                                      %             satellite_orbits.m)
        
        % CONSTELLATION SPECIFIC ------------------------------------------
        
        FL1 = 1575.420;  % GPS [MHz]
        FL2 = 1227.600;  %
        FL5 = 1176.450;  %
        
        FR1_base  = 1602.000;  % GLONASS [MHz]
        FR2_base  = 1246.000;  %
        FR1_delta = 0.5625;    %
        FR2_delta = 0.4375;    %
        FR_channels = 6:-1:-7;
        
        FE1  = goGNSS.FL1;     % Galileo [MHz]
        FE5a = goGNSS.FL5;     %
        FE5b = 1207.140;       %
        FE5  = 1191.795;       %
        FE6  = 1278.750;       %
        
        FC1  = 1589.740;       % BeiDou [MHz]
        FC2  = 1561.098;       %
        FC5b = goGNSS.FE5b;    %
        FC6  = 1268.520;       %
        
        FJ1 = goGNSS.FL1;      % QZSS [MHz]
        FJ2 = goGNSS.FL2;      %
        FJ5 = goGNSS.FL5;      %
        FJ6 = goGNSS.FE6;      %
        
        FS1 = goGNSS.FL1;      % SBAS [MHz]
        FS5 = goGNSS.FL5;      %
        
        FG = [goGNSS.FL1 goGNSS.FL2 goGNSS.FL5]*1e6;        % GPS carriers frequencies [Hz]
        LAMBDAG = goGNSS.V_LIGHT ./ goGNSS.FG;              % GPS carriers wavelengths [m]
        
        FR_base  = [goGNSS.FR1_base goGNSS.FR2_base];       % GLONASS carriers base frequencies [Hz]
        FR_delta = [goGNSS.FR1_delta goGNSS.FR2_delta];     % GLONASS carriers delta frequencies [Hz/n]
        FR1 = goGNSS.FR_channels' .* goGNSS.FR_delta(1) + goGNSS.FR_base(1);
        FR2 = goGNSS.FR_channels' .* goGNSS.FR_delta(2) + goGNSS.FR_base(2);
        FR = [goGNSS.FR1 goGNSS.FR2]*1e6;                   % GLONASS carriers frequencies [Hz]
        LAMBDAR = goGNSS.V_LIGHT ./ goGNSS.FR;              % GLONASS carriers wavelengths [m]
        
        FE = [goGNSS.FE1 goGNSS.FE5a goGNSS.FE5b goGNSS.FE5 goGNSS.FE6]*1e6; % Galileo carriers frequencies [Hz]
        LAMBDAE = goGNSS.V_LIGHT ./ goGNSS.FE;                               % Galileo carriers wavelengths [m]
        
        FC = [goGNSS.FC1 goGNSS.FC2 goGNSS.FC5b goGNSS.FC6]*1e6; % BeiDou carriers frequencies [Hz]
        LAMBDAC = goGNSS.V_LIGHT ./ goGNSS.FC;                   % BeiDou carriers wavelengths [m]
        
        FJ = [goGNSS.FJ1 goGNSS.FJ2 goGNSS.FJ5 goGNSS.FJ6]*1e6;  % QZSS carriers frequencies [Hz]
        LAMBDAJ = goGNSS.V_LIGHT ./ goGNSS.FJ;                   % QZSS carriers wavelengths [m]
        
        FS = [goGNSS.FS1 goGNSS.FS5]*1e6;                        % SBAS carriers frequencies [Hz]
        LAMBDAS = goGNSS.V_LIGHT ./ goGNSS.FS;                   % SBAS carriers wavelengths [m]
        
        % CONSTELLATIONS IDs ----------------------------------------------
        
        ID_GPS     = 1 % Id of GPS constellation for goGPS internal use
        ID_GLONASS = 2 % Id of GLONASS constellation for goGPS internal use
        ID_GALILEO = 3 % Id of Galileo constellation for goGPS internal use
        ID_BEIDOU  = 4 % Id of BeiDou constellation for goGPS internal use
        ID_QZSS    = 5 % Id of QZSS constellation for goGPS internal use
        ID_SBAS    = 6 % Id of SBAS constellation for goGPS internal use
        
        % goGPS MODES -----------------------------------------------------
        
        MODE_RT_NAV          = 24;  % Real Time Navigation (Kalman Filter on Code and Phase Double Differences (with/without a constraint)
        MODE_RT_R_MON        = 21;  % Real Time Rover Monitor
        MODE_RT_M_MON        = 22;  % Real Time Master Monitor
        MODE_RT_RM_MON       = 23;  % Real Time Master + Rover Monitor
        
        MODE_PP_LS_C_SA      = 1;   % Post Proc Least Squares on Code Stand Alone
        MODE_PP_LS_CP_SA     = 3;   % Post Proc Least Squares on Code and Phase Stand Alone (BASE FOR FUTURE PPP IMPLEMENTATION)
        MODE_PP_LS_CP_VEL    = 3.1; % Post Proc Least Squares on Code and Phase for Velocity estimation
        MODE_PP_LS_C_DD      = 11;  % Post Proc Least Squares on Code Double Differences
        MODE_PP_LS_CP_DD_L   = 13;  % Post Proc Least Squares on Code and Phase Double Differences with LAMBDA
        MODE_PP_LS_CP_DD_MR  = 16;  % Post Proc Least Squares on Code and Phase Double Differences, Multiple Receivers
        MODE_PP_LS_C_SA_MR   = 17;  % Post Proc Least Squares on Code Stand Alone, Multiple Receivers
        
        MODE_PP_KF_C_SA      = 2;   % Post Proc Kalman Filter on Code Stand Alone
        MODE_PP_KF_C_DD      = 12;  % Post Proc Kalman Filter on Code Double Differencies
        MODE_PP_KF_CP_SA     = 4;   % Post Proc Kalman Filter on Code and Phase Stand Alone
        MODE_PP_KF_CP_DD     = 14;  % Post Proc Kalman Filter on Code and Phase Double Differences
        MODE_PP_KF_CP_DD_MR  = 15;  % Post Proc Kalman Filter on Code and Phase Double Differences, Multiple Receivers
                 
        GMODE_PP = [ goGNSS.MODE_PP_LS_C_SA ...     % Group of post processing modes
            goGNSS.MODE_PP_LS_CP_SA ...
            goGNSS.MODE_PP_LS_C_DD ...
            goGNSS.MODE_PP_LS_CP_DD_L ...
            goGNSS.MODE_PP_LS_CP_VEL ...
            goGNSS.MODE_PP_KF_C_SA ...
            goGNSS.MODE_PP_KF_C_DD ...
            goGNSS.MODE_PP_KF_CP_SA ...
            goGNSS.MODE_PP_KF_CP_DD ...
            goGNSS.MODE_PP_LS_C_SA_MR ...
            goGNSS.MODE_PP_LS_CP_DD_MR ...
            goGNSS.MODE_PP_KF_CP_DD_MR];
        
        GMODE_RT = [ goGNSS.MODE_RT_NAV ...         % Group of real time modes
            goGNSS.MODE_RT_R_MON ...
            goGNSS.MODE_RT_M_MON ...
            goGNSS.MODE_RT_RM_MON];
        
        GMODE_SA = [ goGNSS.MODE_PP_LS_C_SA ...     % Group of stand alone modes
            goGNSS.MODE_PP_LS_CP_SA ...
            goGNSS.MODE_PP_LS_CP_VEL ...
            goGNSS.MODE_PP_KF_C_SA ...
            goGNSS.MODE_PP_KF_CP_SA ...
            goGNSS.MODE_PP_LS_C_SA_MR];
        
        GMODE_DD = [ goGNSS.MODE_PP_LS_C_DD ...     % Group of double differences modes
            goGNSS.MODE_PP_LS_CP_DD_L ...
            goGNSS.MODE_PP_KF_C_DD ...
            goGNSS.MODE_PP_KF_CP_DD ...
            goGNSS.MODE_PP_LS_CP_DD_MR ...
            goGNSS.MODE_PP_KF_CP_DD_MR];
        
        GMODE_MR = [ goGNSS.MODE_PP_LS_C_SA_MR ...  % Group of multi-receiver modes
            goGNSS.MODE_PP_LS_CP_DD_MR ...
            goGNSS.MODE_PP_KF_CP_DD_MR];
       
        GMODE_PH = [ goGNSS.MODE_RT_NAV ...         % Group of modes using Phase
            goGNSS.MODE_RT_R_MON ...
            goGNSS.MODE_RT_M_MON ...
            goGNSS.MODE_RT_RM_MON ...
            goGNSS.MODE_PP_LS_CP_SA ...
            goGNSS.MODE_PP_LS_CP_VEL ...
            goGNSS.MODE_PP_LS_CP_DD_MR ...
            goGNSS.MODE_PP_LS_CP_DD_L ...
            goGNSS.MODE_PP_KF_CP_SA ...
            goGNSS.MODE_PP_KF_CP_DD ...
            goGNSS.MODE_PP_KF_CP_DD_MR];
        
        GMODE_KM = [ goGNSS.MODE_PP_KF_C_SA ...      % Group of modes using Kalman Filter
            goGNSS.MODE_PP_KF_C_DD ... 
            goGNSS.MODE_PP_KF_CP_SA ...
            goGNSS.MODE_PP_KF_CP_DD ...
            goGNSS.MODE_PP_KF_CP_DD_MR];
        
    end
    
    % Creator (empty)
    methods
        % Object to manage a useful functions / standard parameters of goGPS
        function obj = goGNSS()
        end
    end
    
    %   COMPATIBILITY FUNCTION (STATIC)
    % -------------------------------------------------------------------------
    % function to keep compatibility with the past
    methods (Static, Access = 'public')
        function f1 = F1()
            % GPS carriers frequencies F1
            f1 = goGNSS.FG(1);
        end
        function f1 = F2()
            % GPS carriers frequencies F2
            f1 = goGNSS.FG(2);
        end
        
        function lambda1 = LAMBDA1()
            % GPS carriers frequency 1 [Hz]
            lambda1 = goGNSS.LAMBDAG(1);
        end
        function lambda2 = LAMBDA2()
            % GPS carriers frequency 2 [Hz]
            lambda2 = goGNSS.LAMBDAG(2);
            
        end
    end
    
    %   MODE FUNCTION (STATIC)
    % -------------------------------------------------------------------------
    % function to detect a certain kind of processing
    methods (Static, Access = 'public')
        function isPostProcessing = isPP(mode)
            % return whether or not the mode given in use is a Post Processing mode
            isPostProcessing = sum(intersect(mode, goGNSS.GMODE_PP));
        end
        
        function isRealTime = isRT(mode)
            % return whether or not the mode given in use is a Real Time mode
            isRealTime = sum(intersect(mode, goGNSS.GMODE_RT));
        end
        
        function isDoubleDifferences = isDD(mode)
            % return whether or not the mode given in use is a Double Difference mode
            isDoubleDifferences = sum(intersect(mode, goGNSS.GMODE_DD));
        end
        
        function isStandAlone = isSA(mode)
            % return whether or not the mode given in use is a Stand Alone mode
            isStandAlone = sum(intersect(mode, goGNSS.GMODE_SA));
        end
        
        function isMultiReceiver = isMR(mode)
            % return whether or not the mode given in use is a Stand Alone mode
            isMultiReceiver = sum(intersect(mode, goGNSS.GMODE_MR));
        end
        
        function isUsingPhase = isPH(mode)
            % return whether or not the mode given in use uses Phase
            isUsingPhase = sum(intersect(mode, goGNSS.GMODE_PH));
        end
        
         function isKalman = isKM(mode)
            % return whether or not the mode given in use uses Kalman Filter
            isKalman = sum(intersect(mode, goGNSS.GMODE_KM));
        end       
        
    end
    
    
    %   CONSTELLATION FUNCTION (STATIC)
    % -------------------------------------------------------------------------
    methods (Static, Access = 'public')
        
        function [constellations] = initConstellation(GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag)
            
            % SYNTAX:
            %   [constellations] = initConstellation(GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag);
            %
            % INPUT:
            %   GPS_flag = boolean flag for enabling/disabling GPS usage
            %   GLO_flag = boolean flag for enabling/disabling GLONASS usage
            %   GAL_flag = boolean flag for enabling/disabling Galileo usage
            %   BDS_flag = boolean flag for enabling/disabling BeiDou usage
            %   QZS_flag = boolean flag for enabling/disabling QZSS usage
            %   SBS_flag = boolean flag for enabling/disabling SBAS usage (for ranging)
            %
            % OUTPUT:
            %   constellations = struct with multi-constellation settings
            %
            % DESCRIPTION:
            %   Multi-constellation settings and initialization.
            
            GPS_PRN = [1:32];
            GLO_PRN = [1:24];
            GAL_PRN = [1:30];
            BDS_PRN = [1:37];
            QZS_PRN = [193:196];
            SBS_PRN = 0; %SBAS ranging not supported yet
            
            constellations.GPS     = struct('numSat', numel(GPS_PRN), 'enabled', GPS_flag, 'indexes', 0, 'PRN', GPS_PRN);
            constellations.GLONASS = struct('numSat', numel(GLO_PRN), 'enabled', GLO_flag, 'indexes', 0, 'PRN', GLO_PRN);
            constellations.Galileo = struct('numSat', numel(GAL_PRN), 'enabled', GAL_flag, 'indexes', 0, 'PRN', GAL_PRN);
            constellations.BeiDou  = struct('numSat', numel(BDS_PRN), 'enabled', BDS_flag, 'indexes', 0, 'PRN', BDS_PRN);
            constellations.QZSS    = struct('numSat', numel(QZS_PRN), 'enabled', QZS_flag, 'indexes', 0, 'PRN', QZS_PRN);
            constellations.SBAS    = struct('numSat', numel(SBS_PRN), 'enabled', 0,        'indexes', 0, 'PRN', SBS_PRN); %SBAS ranging not supported yet
            
            nSatTot = 0; %total number of satellites used given the enabled constellations
            q = 0;       %counter for enabled constellations
            
            systems = fieldnames(constellations);
            constellations.indexes = [];
            constellations.PRN = [];
            for i = 1 : numel(systems)
                if(constellations.(systems{i}).enabled)
                    nSatTot = nSatTot + constellations.(systems{i}).numSat;
                    q = q + 1;
                    if (q == 1)
                        indexes_tmp = [1 : constellations.(systems{i}).numSat];
                    else
                        indexes_tmp = [indexes_tmp(end) + 1 : indexes_tmp(end) + constellations.(systems{i}).numSat];
                    end
                    constellations.(systems{i}).indexes = indexes_tmp;
                    constellations.indexes = [constellations.indexes, indexes_tmp];
                    constellations.PRN = [constellations.PRN, constellations.(systems{i}).PRN];
                end
            end
            
            constellations.nEnabledSat = nSatTot;
        end
        
        function [lambda] = getGNSSWavelengths(Eph, nSatTot)
            lambda = zeros(nSatTot,2);
            for s = 1 : nSatTot
                pos = find(Eph(30,:) == s,1);
                if (~isempty(pos))
                    switch char(Eph(31,pos))
                        case 'G'
                            lambda(s,1) = goGNSS.getWavelength(goGNSS.ID_GPS, 1);
                            lambda(s,2) = goGNSS.getWavelength(goGNSS.ID_GPS, 2);
                        case 'R'
                            lambda(s,1) = goGNSS.getWavelength(goGNSS.ID_GLONASS, 1, Eph(15,pos));
                            lambda(s,2) = goGNSS.getWavelength(goGNSS.ID_GLONASS, 2, Eph(15,pos));
                        case 'E'
                            lambda(s,1) = goGNSS.getWavelength(goGNSS.ID_GALILEO, 1);
                            lambda(s,2) = goGNSS.getWavelength(goGNSS.ID_GALILEO, 2);
                        case 'C'
                            lambda(s,1) = goGNSS.getWavelength(goGNSS.ID_BEIDOU, 2);
                            lambda(s,2) = goGNSS.getWavelength(goGNSS.ID_BEIDOU, 3);
                        case 'J'
                            lambda(s,1) = goGNSS.getWavelength(goGNSS.ID_QZSS, 1);
                            lambda(s,2) = goGNSS.getWavelength(goGNSS.ID_QZSS, 2);
                        otherwise
                            fprintf('Something went wrong in goGNSS.getGNSSWavelengths()\nUnrecongized Satellite system.\n');
                    end
                end
            end
        end
        
        function [lambda] = getWavelength(idCostellation, freq, GLOFreqNum)
            switch idCostellation
                case goGNSS.ID_GPS
                    lambda = goGNSS.LAMBDAG(freq);
                case goGNSS.ID_GLONASS
                    pos = goGNSS.FR_channels == GLOFreqNum;
                    lambda = goGNSS.LAMBDAR(pos, freq);
                case goGNSS.ID_GALILEO
                    lambda = goGNSS.LAMBDAE(freq);
                case goGNSS.ID_BEIDOU
                    lambda = goGNSS.LAMBDAC(freq);
                case goGNSS.ID_QZSS
                    lambda = goGNSS.LAMBDAJ(freq);
                case goGNSS.ID_SBAS
                    lambda = goGNSS.LAMBDAS(freq);
            end
        end
        
        function [GLOFreqNum] = getGLOFreqNum(satNum, Eph)
            pos1 = find(Eph(30,:) == satNum,1);
            pos2 = find(strcmp(char(Eph(31,:)),'R'));
            pos = interesct(pos1, pos2);
            GLOFreqNum = Eph(15,pos);
        end
        
        function [IonoFactor] = getInterFreqIonoFactor(lambda)
            IonoFactor = (lambda(:,:)./goGNSS.LAMBDAG(1)).^2;
        end
    end
    
    %   USEFUL FUNCTION (STATIC)
    % -------------------------------------------------------------------------
    methods (Static, Access = 'public')
        
        function [XR, dtR] = getBancroftPos(XS, dtS, prR)
            % get a bancroft solution from one receiver
            
            matB = [XS(:,:), prR(:) + goGNSS.V_LIGHT * dtS(:)]; % Bancroft matrix
            b = goGNSS.bancroft(matB);
            XR = b(1:3);
            dtR = b(4)/goGNSS.V_LIGHT;
        end
        
        function RinReader(fileNameObs, constellations)
            if (isempty(constellations)) %then use only GPS as default
                [constellations] = goGNSS.initConstellation(1, 0, 0, 0, 0, 0);
            end
            tic;            
            %number of satellite slots for enabled constellations
            nSatTot = constellations.nEnabledSat;

            %open RINEX observation file (ROVER)
            fid = fopen(fileNameObs,'r');
            txtRin = textscan(fid,'%s','Delimiter','\n','whitespace','');
            fclose(fid);
            txtRin = txtRin{1};
            
            %parse RINEX header
            [version obsTypes knownPos flagFoundTypes interval sysId line] = goGNSS.RinParseHDR(txtRin);
            [obsCols, nType] = obs_type_find(obsTypes, sysId);
            
        end
        
        function [pr1_R, pr1_M, ph1_R, ph1_M, pr2_R, pr2_M, ph2_R, ph2_M, ...
                dop1_R, dop1_M, dop2_R, dop2_M, snr1_R, snr1_M, ...
                snr2_R, snr2_M, time, time_R, time_M, week_R, week_M, ...
                date_R, date_M, pos_R, pos_M, Eph, iono, interval] = ...
                loadRINEX(filename_nav, filename_R_obs, filename_M_obs, constellations, flag_SP3, wait_dlg)
            
            % Check the input arguments
            if (nargin < 6)
                wait_dlg_PresenceFlag = false;
            else
                wait_dlg_PresenceFlag = true;
            end
            if (isempty(filename_M_obs))
                filename_M_obs_PresenceFlag = false;
            else
                filename_M_obs_PresenceFlag = true;
            end
            if (isempty(constellations)) %then use only GPS as default
                [constellations] = goGNSS.initConstellation(1, 0, 0, 0, 0, 0);
            end
            tic;
            
            %number of satellite slots for enabled constellations
            nSatTot = constellations.nEnabledSat;
            
            %fraction of INTERVAL (epoch-to-epoch timespan, as specified in the header)
            %that is allowed as maximum difference between rover and master timings
            %during synchronization
            max_desync_frac = 0.1;
            
            %read navigation files
            if (~flag_SP3)
                if (wait_dlg_PresenceFlag)
                    waitbar(0.5,wait_dlg,'Reading navigation files...')
                end
                
                Eph_G = []; iono_G = zeros(8,1);
                Eph_R = []; iono_R = zeros(8,1);
                Eph_E = []; iono_E = zeros(8,1);
                Eph_C = []; iono_C = zeros(8,1);
                Eph_J = []; iono_J = zeros(8,1);
                
                if (strcmpi(filename_nav(end),'p'))
                    flag_mixed = 1;
                else
                    flag_mixed = 0;
                end
                
                if (constellations.GPS.enabled || flag_mixed)
                    if (exist(filename_nav,'file'))
                        %parse RINEX navigation file (GPS) NOTE: filename expected to
                        %end with 'n' or 'N' (GPS) or with 'p' or 'P' (mixed GNSS)
                        [Eph_G, iono_G] = RINEX_get_nav(filename_nav, constellations);
                    else
                        fprintf('... WARNING: GPS navigation file not found. Disabling GPS positioning. \n');
                        constellations.GPS.enabled = 0;
                    end
                end
                
                if (constellations.GLONASS.enabled)
                    if (exist([filename_nav(1:end-1) 'g'],'file'))
                        %parse RINEX navigation file (GLONASS)
                        [Eph_R, iono_R] = RINEX_get_nav([filename_nav(1:end-1) 'g'], constellations);
                    elseif (~flag_mixed)
                        fprintf('... WARNING: GLONASS navigation file not found. Disabling GLONASS positioning. \n');
                        constellations.GLONASS.enabled = 0;
                    end
                end
                
                if (constellations.Galileo.enabled)
                    if (exist([filename_nav(1:end-1) 'l'],'file'))
                        %parse RINEX navigation file (Galileo)
                        [Eph_E, iono_E] = RINEX_get_nav([filename_nav(1:end-1) 'l'], constellations);
                    elseif (~flag_mixed)
                        fprintf('... WARNING: Galileo navigation file not found. Disabling Galileo positioning. \n');
                        constellations.Galileo.enabled = 0;
                    end
                end
                
                if (constellations.BeiDou.enabled)
                    if (exist([filename_nav(1:end-1) 'b'],'file'))
                        parse RINEX navigation file (BeiDou)
                        [Eph_C, iono_C] = RINEX_get_nav([filename_nav(1:end-1) 'b'], constellations);
                    elseif (~flag_mixed)
                        fprintf('... WARNING: BeiDou navigation file not found. Disabling BeiDou positioning. \n');
                        constellations.BeiDou.enabled = 0;
                    end
                end
                
                if (constellations.QZSS.enabled)
                    if (exist([filename_nav(1:end-1) 'q'],'file'))
                        %parse RINEX navigation file (QZSS)
                        [Eph_J, iono_J] = RINEX_get_nav([filename_nav(1:end-1) 'q'], constellations);
                    elseif (~flag_mixed)
                        fprintf('... WARNING: QZSS navigation file not found. Disabling QZSS positioning. \n');
                        constellations.QZSS.enabled = 0;
                    end
                end
                
                Eph = [Eph_G Eph_R Eph_E Eph_C Eph_J];
                
                if (any(iono_G))
                    iono = iono_G;
                elseif (any(iono_R))
                    iono = iono_R;
                elseif (any(iono_E))
                    iono = iono_E;
                elseif (any(iono_C))
                    iono = iono_C;
                elseif (any(iono_J))
                    iono = iono_J;
                else
                    iono = zeros(8,1);
                    fprintf('... WARNING: ionosphere parameters not found in navigation file(s).\n');
                end
                
                if (wait_dlg_PresenceFlag)
                    waitbar(1,wait_dlg)
                end
            else
                Eph = zeros(33,nSatTot);
                iono = zeros(8,1);
            end
            
            %-------------------------------------------------------------------------------
            
            %open RINEX observation file (ROVER)
            FR_oss = fopen(filename_R_obs,'r');
            txtRin = textscan(FR_oss,'%s','Delimiter','\n','whitespace','');
            txtRin = txtRin{1};
            fclose(FR_oss);
            
            if (filename_M_obs_PresenceFlag)
                %open RINEX observation file (MASTER)
                FM_oss = fopen(filename_M_obs,'r');
            end
            
            %-------------------------------------------------------------------------------
            
            if (wait_dlg_PresenceFlag)
                waitbar(0.5,wait_dlg,'Parsing RINEX headers...')
            end
            
            %parse RINEX header
            [RINEX_version obs_typ_R, pos_R, info_base_R, interval_R, sysId, line] = goGNSS.RinParseHDR(txtRin);
            
            %check RINEX version
            if (isempty(sysId))
                RINEX_version = 2;
            else
                RINEX_version = 3;
            end
            
            %check the availability of basic data to parse the RINEX file (ROVER)
            if (info_base_R == 0)
                error('Basic data is missing in the ROVER RINEX header')
            end
            
            %find observation type columns
            [obs_col_R, nObsTypes_R] = goGNSS.obs_type_find(obs_typ_R, sysId);
            
            %number of lines to be read for each epoch (only for RINEX v2.xx)
            if (RINEX_version == 2)
                nLinesToRead_R = ceil(nObsTypes_R/5);  %maximum of 5 obs per line
            end
            
            if (filename_M_obs_PresenceFlag)
                [obs_typ_M, pos_M, info_base_M, interval_M, sysId] = RINEX_parse_hdr(FM_oss);
                
                %check the availability of basic data to parse the RINEX file (MASTER)
                if (info_base_M == 0)
                    error('Basic data is missing in the ROVER RINEX header')
                end
                
                %find observation type columns
                [obs_col_M, nObsTypes_M] = obs_type_find(obs_typ_M, sysId);
                
                %number of lines to be read for each epoch (only for RINEX v2.xx)
                if (~isstruct(nObsTypes_M))
                    nLinesToRead_M = ceil(nObsTypes_M/5);  %maximum of 5 obs per line
                end
            else
                pos_M = zeros(3,1);
                interval_M = [];
            end
            
            if (wait_dlg_PresenceFlag)
                waitbar(1,wait_dlg)
            end
            
            interval = min([interval_R, interval_M]);
            
            %-------------------------------------------------------------------------------
            
            nEpochs = 10800;
            
            %variable initialization (GPS)
            time_R = zeros(nEpochs,1);
            time_M = zeros(nEpochs,1);
            pr1_R = zeros(nSatTot,nEpochs);
            pr2_R = zeros(nSatTot,nEpochs);
            ph1_R = zeros(nSatTot,nEpochs);
            ph2_R = zeros(nSatTot,nEpochs);
            dop1_R = zeros(nSatTot,nEpochs);
            dop2_R = zeros(nSatTot,nEpochs);
            snr1_R = zeros(nSatTot,nEpochs);
            snr2_R = zeros(nSatTot,nEpochs);
            pr1_M = zeros(nSatTot,nEpochs);
            pr2_M = zeros(nSatTot,nEpochs);
            ph1_M = zeros(nSatTot,nEpochs);
            ph2_M = zeros(nSatTot,nEpochs);
            snr1_M = zeros(nSatTot,nEpochs);
            snr2_M = zeros(nSatTot,nEpochs);
            dop1_M = zeros(nSatTot,nEpochs);
            dop2_M = zeros(nSatTot,nEpochs);
            date_R = zeros(nEpochs,6);
            date_M = zeros(nEpochs,6);
            
            %read data for the first epoch (ROVER)
            [time_R(1), epoch_R, num_sat_R, sat_R, sat_types_R, line] = goGNSS.RinGetEpoch(txtRin, line);
            
            %-------------------------------------------------------------------------------
            
            if (filename_M_obs_PresenceFlag)
                %read data for the first epoch (MASTER)
                [time_M(1), epoch_M, num_sat_M, sat_M, sat_types_M] = RINEX_get_epoch(FM_oss);
            end
            %-------------------------------------------------------------------------------
            
            if (wait_dlg_PresenceFlag)
                waitbar(0.5,wait_dlg,'Parsing RINEX headers...')
            end
            
            if (filename_M_obs_PresenceFlag)
                while ((time_M(1) - time_R(1)) < 0 && abs(time_M(1) - time_R(1)) >= max_desync_frac*interval)
                    
                    %number of lines to be skipped
                    if (RINEX_version == 2)
                        nSkipLines = num_sat_M*nLinesToRead_M;
                    else
                        nSkipLines = num_sat_M;
                    end
                    
                    %skip observations
                    for s = 1 : nSkipLines
                        fgetl(FM_oss);
                    end
                    
                    %read data for the current epoch (MASTER)
                    [time_M(1), epoch_M, num_sat_M, sat_M, sat_types_M] = RINEX_get_epoch(FM_oss);
                end
                
                while ((time_R(1) - time_M(1)) < 0 && abs(time_R(1) - time_M(1)) >= max_desync_frac*interval)
                    
                    %number of lines to be skipped
                    if (RINEX_version == 2)
                        nSkipLines = num_sat_R*nLinesToRead_R;
                    else
                        nSkipLines = num_sat_R;
                    end
                    
                    %skip observations
                    for s = 1 : nSkipLines
                        fgetl(FR_oss);
                    end
                    
                    %read data for the current epoch (ROVER)
                    [time_R(1), epoch_R, num_sat_R, sat_R, sat_types_R] = RINEX_get_epoch(FR_oss);
                end
            end
            
            %read first batch of observations
            %ROVER
            [obs_R, line] = goGNSS.RinGetObs(txtRin, line, num_sat_R, sat_R, sat_types_R, obs_col_R, nObsTypes_R, constellations);
            
            %read ROVER observations
            if (sum(obs_R.P1 ~= 0) == constellations.nEnabledSat)
                pr1_R(:,1) = obs_R.P1;
            else
                pr1_R(:,1) = obs_R.C1;
            end
            pr2_R(:,1) = obs_R.P2;
            ph1_R(:,1) = obs_R.L1;
            ph2_R(:,1) = obs_R.L2;
            dop1_R(:,1) = obs_R.D1;
            dop2_R(:,1) = obs_R.D2;
            snr1_R(:,1) = obs_R.S1;
            snr2_R(:,1) = obs_R.S2;
            
            if (filename_M_obs_PresenceFlag)
                %MASTER
                obs_M = RINEX_get_obs(FM_oss, num_sat_M, sat_M, sat_types_M, obs_col_M, nObsTypes_M, constellations);
                
                %read MASTER observations
                if (sum(obs_M.P1 ~= 0) == constellations.nEnabledSat)
                    pr1_M(:,1) = obs_M.P1;
                else
                    pr1_M(:,1) = obs_M.C1;
                end
                pr2_M(:,1) = obs_M.P2;
                ph1_M(:,1) = obs_M.L1;
                ph2_M(:,1) = obs_M.L2;
                dop1_M(:,1) = obs_M.D1;
                dop2_M(:,1) = obs_M.D2;
                snr1_M(:,1) = obs_M.S1;
                snr2_M(:,1) = obs_M.S2;
            end
            
            if (wait_dlg_PresenceFlag)
                waitbar(1,wait_dlg)
            end
            
            %-------------------------------------------------------------------------------
            
            %define the reference time
            time(1,1) = roundmod(time_R(1),interval);
            date_R(1,:) = epoch_R(1,:);
            if (filename_M_obs_PresenceFlag)
                date_M(1,:) = epoch_M(1,:);
            end
            
            if (wait_dlg_PresenceFlag)
                waitbar(0.5,wait_dlg,'Reading RINEX observations...')
            end
            
            k = 2;
            while (line < length(txtRin))
                
                if (abs((time_R(k-1) - time(k-1))) < max_desync_frac*interval)
                    %read data for the current epoch (ROVER)
                    [time_R(k), epoch_R, num_sat_R, sat_R, sat_types_R, line] = goGNSS.RinGetEpoch(txtRin, line);
                else
                    time_R(k) = time_R(k-1);
                    if (time_R(k-1) ~= 0)
                        fprintf('Missing epoch %f (ROVER)\n', time(k-1));
                    end
                    time_R(k-1) = 0;
                end
                
                if (filename_M_obs_PresenceFlag)
                    if (abs((time_M(k-1) - time(k-1))) < max_desync_frac*interval)
                        %read data for the current epoch (MASTER)
                        [time_M(k), epoch_M, num_sat_M, sat_M, sat_types_M] = RINEX_get_epoch(FM_oss);
                    else
                        time_M(k) = time_M(k-1);
                        if (time_M(k-1) ~= 0)
                            fprintf('Missing epoch %f (MASTER)\n', time(k-1));
                        end
                        time_M(k-1) = 0;
                    end
                end
                
                if (k > nEpochs)
                    %variable initialization (GPS)
                    pr1_R(:,k) = zeros(nSatTot,1);
                    pr2_R(:,k) = zeros(nSatTot,1);
                    ph1_R(:,k) = zeros(nSatTot,1);
                    ph2_R(:,k) = zeros(nSatTot,1);
                    dop1_R(:,k) = zeros(nSatTot,1);
                    dop2_R(:,k) = zeros(nSatTot,1);
                    snr1_R(:,k) = zeros(nSatTot,1);
                    snr2_R(:,k) = zeros(nSatTot,1);
                    pr1_M(:,k) = zeros(nSatTot,1);
                    pr2_M(:,k) = zeros(nSatTot,1);
                    ph1_M(:,k) = zeros(nSatTot,1);
                    ph2_M(:,k) = zeros(nSatTot,1);
                    snr1_M(:,k) = zeros(nSatTot,1);
                    snr2_M(:,k) = zeros(nSatTot,1);
                    dop1_M(:,k) = zeros(nSatTot,1);
                    dop2_M(:,k) = zeros(nSatTot,1);
                    
                    nEpochs = nEpochs  + 1;
                end
                
                date_R(k,:) = epoch_R(1,:);
                if (filename_M_obs_PresenceFlag)
                    date_M(k,:) = epoch_M(1,:);
                end
                
                time(k,1) = time(k-1,1) + interval;
                
                if (abs(time_R(k)-time(k)) < max_desync_frac*interval)
                    
                    %read ROVER observations
                    [obs_R line] = goGNSS.RinGetObs(txtRin, line, num_sat_R, sat_R, sat_types_R, obs_col_R, nObsTypes_R, constellations);
                    
                    %read ROVER observations
                    if (sum(obs_R.P1 ~= 0) == constellations.nEnabledSat)
                        pr1_R(:,k) = obs_R.P1;
                    else
                        pr1_R(:,k) = obs_R.C1;
                    end
                    pr2_R(:,k) = obs_R.P2;
                    ph1_R(:,k) = obs_R.L1;
                    ph2_R(:,k) = obs_R.L2;
                    dop1_R(:,k) = obs_R.D1;
                    dop2_R(:,k) = obs_R.D2;
                    snr1_R(:,k) = obs_R.S1;
                    snr2_R(:,k) = obs_R.S2;
                    %     else
                    %         %number of lines to be skipped
                    %         if (RINEX_version == 2)
                    %             nSkipLines = num_sat_R*nLinesToRead_R;
                    %         else
                    %             nSkipLines = num_sat_R;
                    %         end
                    %
                    %         %skip observations
                    %         for s = 1 : nSkipLines
                    %             fgetl(FR_oss);
                    %         end
                end
                
                if (filename_M_obs_PresenceFlag)
                    
                    if (abs(time_M(k) - time(k)) < max_desync_frac*interval)
                        
                        %read MASTER observations
                        obs_M = RINEX_get_obs(FM_oss, num_sat_M, sat_M, sat_types_M, obs_col_M, nObsTypes_M, constellations);
                        
                        %read MASTER observations
                        if (sum(obs_M.P1 ~= 0) == constellations.nEnabledSat)
                            pr1_M(:,k) = obs_M.P1;
                        else
                            pr1_M(:,k) = obs_M.C1;
                        end
                        pr2_M(:,k) = obs_M.P2;
                        ph1_M(:,k) = obs_M.L1;
                        ph2_M(:,k) = obs_M.L2;
                        dop1_M(:,k) = obs_M.D1;
                        dop2_M(:,k) = obs_M.D2;
                        snr1_M(:,k) = obs_M.S1;
                        snr2_M(:,k) = obs_M.S2;
                        %         else
                        %             %number of lines to be skipped
                        %             if (RINEX_version == 2)
                        %                 nSkipLines = num_sat_M*nLinesToRead_M;
                        %             else
                        %                 nSkipLines = num_sat_M;
                        %             end
                        %
                        %             %skip observations
                        %             for s = 1 : nSkipLines
                        %                 fgetl(FM_oss);
                        %             end
                    end
                end
                
                k = k + 1;
            end
            
            %remove empty slots
            time_R(k:nEpochs) = [];
            time_M(k:nEpochs) = [];
            pr1_R(:,k:nEpochs) = [];
            pr2_R(:,k:nEpochs) = [];
            ph1_R(:,k:nEpochs) = [];
            ph2_R(:,k:nEpochs) = [];
            dop1_R(:,k:nEpochs) = [];
            dop2_R(:,k:nEpochs) = [];
            snr1_R(:,k:nEpochs) = [];
            snr2_R(:,k:nEpochs) = [];
            pr1_M(:,k:nEpochs) = [];
            pr2_M(:,k:nEpochs) = [];
            ph1_M(:,k:nEpochs) = [];
            ph2_M(:,k:nEpochs) = [];
            snr1_M(:,k:nEpochs) = [];
            snr2_M(:,k:nEpochs) = [];
            dop1_M(:,k:nEpochs) = [];
            dop2_M(:,k:nEpochs) = [];
            date_R(k:nEpochs,:) = [];
            date_M(k:nEpochs,:) = [];
            
            %remove rover tail
            if (filename_M_obs_PresenceFlag)
                flag_tail = 1;
                while (flag_tail)
                    if (time_M(end) == 0)
                        date_R(end,:) = [];
                        date_M(end,:) = [];
                        time(end) = [];
                        time_R(end) = [];
                        time_M(end) = [];
                        pr1_R(:,end) = [];
                        pr2_R(:,end) = [];
                        ph1_R(:,end) = [];
                        ph2_R(:,end) = [];
                        dop1_R(:,end) = [];
                        dop2_R(:,end) = [];
                        snr1_R(:,end) = [];
                        snr2_R(:,end) = [];
                        pr1_M(:,end) = [];
                        pr2_M(:,end) = [];
                        ph1_M(:,end) = [];
                        ph2_M(:,end) = [];
                        snr1_M(:,end) = [];
                        snr2_M(:,end) = [];
                        dop1_M(:,end) = [];
                        dop2_M(:,end) = [];
                    else
                        flag_tail = 0;
                    end
                end
            end
            
            if (wait_dlg_PresenceFlag)
                waitbar(1,wait_dlg)
            end
            
            %-------------------------------------------------------------------------------
            
            %close RINEX files
            if (filename_M_obs_PresenceFlag)
                fclose(FM_oss);
            end
            
            %GPS week number
            week_R = date2gps(date_R);
            week_M = date2gps(date_M);
            fprintf('The RINEX file has been read!\n');
            toc
        end
    end
    
    methods (Static, Access = 'private')
        % Bancroft algorithm for the computation of ground coordinates
        % having at least 4 visible satellites.
        function pos = bancroft(matB)
            % SYNTAX:
            %   [pos] = bancroft(B_pass);
            %
            % INPUT:
            %   B_pass = Bancroft matrix
            %
            % OUTPUT:
            %   pos = approximated ground position (X,Y,Z coordinates)
            %
            % DESCRIPTION:
            %   Bancroft algorithm for the computation of ground coordinates
            %   having at least 4 visible satellites.
            %
            % Copyright (C) Kai Borre
            % Kai Borre 04-30-95, improved by C.C. Goad 11-24-96
            %
            % Adapted by Mirko Reguzzoni, Eugenio Realini, 2009
            %----------------------------------------------------------------------------------------------
            
            pos = zeros(4,1);   % Init position of the receiver
            
            nSat = size(matB,1);
            
            for iter = 1:2
                % Compute correctionon XS coordinates
                % due to the rotation of the Earth
                for s = 1:nSat
                    x = matB(s,1);      % X coordinate of the satellite S
                    y = matB(s,2);      % Y coordinate of the satellite S
                    if iter == 1
                        traveltime = 0.072; % [s]
                    else
                        z = matB(s,3);  % Z coordinate of the satellite S
                        % Compute approximate distance^2 between R and S
                        rho = (x-pos(1))^2+(y-pos(2))^2+(z-pos(3))^2;
                        traveltime = sqrt(rho)/goGNSS.V_LIGHT;
                    end
                    angle = traveltime*goGNSS.OMEGAE_DOT;
                    cosa = cos(angle);
                    sina = sin(angle);
                    % Applying correction (rotation)
                    matB(s,1) =	cosa*x + sina*y;
                    matB(s,2) = -sina*x + cosa*y;
                end % i-loop
                
                % Criptical way to implement Bancroft
                if nSat > 4
                    BBB = (matB'*matB)\matB';
                else
                    BBB = inv(matB);
                end
                e = ones(nSat,1);
                alpha = zeros(nSat,1);
                for s = 1:nSat
                    alpha(s) = lorentz(matB(s,:)',matB(s,:)')/2;
                end
                BBBe = BBB*e;
                BBBalpha = BBB*alpha;
                a = lorentz(BBBe,BBBe);
                b = lorentz(BBBe,BBBalpha)-1;
                c = lorentz(BBBalpha,BBBalpha);
                root = sqrt(b*b-a*c);
                r(1) = (-b-root)/a;
                r(2) = (-b+root)/a;
                possible_pos = zeros(4,2);
                for s = 1:2
                    possible_pos(:,s) = r(s)*BBBe+BBBalpha;
                    possible_pos(4,s) = -possible_pos(4,s);
                end
                
                abs_omc = zeros(2,1);
                for s =1:nSat
                    for i = 1:2
                        c_dt = possible_pos(4,i);
                        calc = norm(matB(s,1:3)' -possible_pos(1:3,i))+c_dt;
                        omc = matB(s,4)-calc;
                        abs_omc(i) = abs(omc);
                    end
                end; % j-loop
                
                % discrimination between roots
                if abs_omc(1) > abs_omc(2)
                    pos = possible_pos(:,2);
                else
                    pos = possible_pos(:,1);
                end
            end;
        end
        
        function [version, obsTypes, knownPos, flagFoundTypes, interval, sysId, curLine] = RinParseHDR(txtRin)
            flagFoundTypes = 0;
            obsTypes = cell(0,0);
            sysId = cell(0,0);
            knownPos = [];
            interval = 1; %default to 1 second (1 Hz observations)
            version = 2;
            
            %parse first line
            curLine = 1; txtLine = txtRin{curLine};
            
            %constellation counter for RINEX v3.xx
            c = 1;
            
            %check if the end of the header or the end of the file has been reached
            while isempty(strfind(txtLine,'END OF HEADER'))
                
                answer = strfind(txtLine,'RINEX VERSION'); %RINEX v2.xx
                if ~isempty(answer)
                    version = floor(sscanf(txtLine(1:15),'%f'));
                end    
                
                answer = strfind(txtLine,'# / TYPES OF OBSERV'); %RINEX v2.xx
                if ~isempty(answer)
                    obsTypes{1} = [];
                    nObs = sscanf(txtLine(1:6),'%d');
                    nLinObs = ceil(nObs/9);
                    for i = 1 : nLinObs
                        if (i > 1)
                            curLine = 1; txtLine = txtRin{curLine};
                        end
                        n = min(nObs,9);
                        for k = 1 : n
                            ot = sscanf(txtLine(k*6+1:k*6+6),'%s');
                            obsTypes{1} = [obsTypes{1} ot];
                        end
                        nObs = nObs - 9;
                    end
                    
                    flagFoundTypes = 1;
                end
                
                answer = strfind(txtLine,'SYS / # / OBS TYPES'); %RINEX v3.xx
                if ~isempty(answer)
                    sysId{c} = sscanf(txtLine(1),'%s');
                    nObs = sscanf(txtLine(2:6),'%d');
                    obsTypes.(sysId{c}) = [];
                    nLinObs = ceil(nObs/13);
                    for i = 1 : nLinObs
                        if (i > 1)
                            curLine = 1; txtLine = txtRin{curLine};
                        end
                        n = min(nObs,13);
                        for k = 0 : n-1
                            ot = sscanf(txtLine(6+k*4+1:6+k*4+4),'%s');
                            obsTypes.(sysId{c}) = [obsTypes.(sysId{c}) ot];
                        end
                        nObs = nObs - 13;
                    end
                    
                    c = c + 1;
                    flagFoundTypes = 1;
                end
                
                answer = strfind(txtLine,'APPROX POSITION XYZ');
                if ~isempty(answer)
                    X = sscanf(txtLine(1:14),'%f');
                    Y = sscanf(txtLine(15:28),'%f');
                    Z = sscanf(txtLine(29:42),'%f');
                    knownPos = [X; Y; Z];
                end
                answer = strfind(txtLine,'INTERVAL');
                if ~isempty(answer)
                    interval = sscanf(txtLine(1:10),'%f');
                end
                
                %parse next line
                curLine = curLine +1; txtLine = txtRin{curLine};
            end
            
            %check RINEX version
            if (~isempty(sysId) && (version == 2))
                version = 3;
            end
        end
        
        function [obsCols, nType] = obs_type_find(obsTypes, sysId)
            
            if (isempty(sysId)) %RINEX v2.xx does not have sysId
                
                nType = size(obsTypes{1},2)/2;
                
                %search L1 column
                s1 = strfind(obsTypes{1}, 'L1');
                s2 = strfind(obsTypes{1}, 'LA');
                col_L1 = [s1 s2];
                
                %search L2 column
                s1 = strfind(obsTypes{1}, 'L2');
                s2 = strfind(obsTypes{1}, 'LC');
                col_L2 = [s1 s2];
                
                %search C1 column
                s1 = strfind(obsTypes{1}, 'C1');
                s2 = strfind(obsTypes{1}, 'CA');
                col_C1 = [s1 s2];
                
                %search P1 column
                s1 = strfind(obsTypes{1}, 'P1');
                s2 = strfind(obsTypes{1}, 'CA'); %QZSS does not use P1
                col_P1 = [s1 s2];
                
                %if RINEX v2.12 and GPS/GLONASS P1 observations are not available
                if (length(col_P1) ~= 2 && ~isempty(s2))
                    %keep QZSS CA observations as C1
                    col_P1 = [];
                end
                
                %search P2 column
                s1 = strfind(obsTypes{1}, 'P2');
                s2 = strfind(obsTypes{1}, 'CC');
                col_P2 = [s1 s2];
                
                %search S1 column
                s1 = strfind(obsTypes{1}, 'S1');
                s2 = strfind(obsTypes{1}, 'SA');
                col_S1 = [s1 s2];
                
                %search S2 column
                s1 = strfind(obsTypes{1}, 'S2');
                s2 = strfind(obsTypes{1}, 'SC');
                col_S2 = [s1 s2];
                
                %search D1 column
                s1 = strfind(obsTypes{1}, 'D1');
                s2 = strfind(obsTypes{1}, 'DA');
                col_D1 = [s1 s2];
                
                %search D2 column
                s1 = strfind(obsTypes{1}, 'D2');
                s2 = strfind(obsTypes{1}, 'DC');
                col_D2 = [s1 s2];
                
                obsCols.L1 = (col_L1+1)/2;
                obsCols.L2 = (col_L2+1)/2;
                obsCols.C1 = (col_C1+1)/2;
                obsCols.P1 = (col_P1+1)/2;
                obsCols.P2 = (col_P2+1)/2;
                obsCols.S1 = (col_S1+1)/2;
                obsCols.S2 = (col_S2+1)/2;
                obsCols.D1 = (col_D1+1)/2;
                obsCols.D2 = (col_D2+1)/2;
                
            else %RINEX v3.xx
                for c = 1 : length(sysId)
                    
                    nType.(sysId{c}) = size(obsTypes.(sysId{c}),2)/3;
                    
                    switch sysId{c}
                        case 'G' %GPS
                            idL1 = 'L1C';
                            idL2 = 'L2W';
                            idC1 = 'C1C';
                            idP1 = 'C1P';
                            idP2 = 'C2W';
                            idS1 = 'S1C';
                            idS2 = 'S2W';
                            idD1 = 'D1C';
                            idD2 = 'D2W';
                        case 'R' %GLONASS
                            idL1 = 'L1C';
                            idL2 = 'L2P';
                            idC1 = 'C1C';
                            idP1 = 'C1P';
                            idP2 = 'C2P';
                            idS1 = 'S1C';
                            idS2 = 'S2P';
                            idD1 = 'D1C';
                            idD2 = 'D2P';
                        case 'E' %Galileo
                            idL1 = 'L1X';
                            idL2 = 'L5X';
                            idC1 = 'C1X';
                            idP1 = '...'; % <-- ?
                            idP2 = 'C5X';
                            idS1 = 'S1X';
                            idS2 = 'S5X';
                            idD1 = 'D1X';
                            idD2 = 'D5X';
                        case 'C' %Compass/Beidou
                            idL1 = 'L2I';
                            idL2 = 'L7I';
                            idC1 = 'C2I';
                            idP1 = '...'; % <-- ?
                            idP2 = 'C7I';
                            idS1 = 'S2I';
                            idS2 = 'S7I';
                            idD1 = 'D2I';
                            idD2 = 'D7I';
                        case 'Q' %QZSS
                            idL1 = 'L1C';
                            idL2 = 'L2C';
                            idC1 = 'C1C';
                            idP1 = '...'; %QZSS does not use P1
                            idP2 = 'C2C';
                            idS1 = 'S1C';
                            idS2 = 'S2C';
                            idD1 = 'D1C';
                            idD2 = 'D2C';
                    end
                    
                    %search L1, L2, C1, P1, P2, S1, S2, D1, D2 columns
                    col_L1 = strfind(obsTypes.(sysId{c}), idL1);
                    col_L2 = strfind(obsTypes.(sysId{c}), idL2);
                    col_C1 = strfind(obsTypes.(sysId{c}), idC1);
                    col_P1 = strfind(obsTypes.(sysId{c}), idP1);
                    col_P2 = strfind(obsTypes.(sysId{c}), idP2);
                    col_S1 = strfind(obsTypes.(sysId{c}), idS1);
                    col_S2 = strfind(obsTypes.(sysId{c}), idS2);
                    col_D1 = strfind(obsTypes.(sysId{c}), idD1);
                    col_D2 = strfind(obsTypes.(sysId{c}), idD2);

                    obsCols.(sysId{c}).L1 = (col_L1+2)/3;
                    obsCols.(sysId{c}).L2 = (col_L2+2)/3;
                    obsCols.(sysId{c}).C1 = (col_C1+2)/3;
                    obsCols.(sysId{c}).P1 = (col_P1+2)/3;
                    obsCols.(sysId{c}).P2 = (col_P2+2)/3;
                    obsCols.(sysId{c}).S1 = (col_S1+2)/3;
                    obsCols.(sysId{c}).S2 = (col_S2+2)/3;
                    obsCols.(sysId{c}).D1 = (col_D1+2)/3;
                    obsCols.(sysId{c}).D2 = (col_D2+2)/3;
                end
            end
        end
        
        function [time, datee, num_sat, sat, sat_types, line] = RinGetEpoch(txtRin, line)
            %variable initialization
            time = 0;
            sat = [];
            sat_types = [];
            num_sat = 0;
            datee=[0 0 0 0 0 0]; %Preallocation not useful (see last line of code)
            if (nargout > 3)
                datee_RequestedInOutputFlag = true;
            else
                datee_RequestedInOutputFlag = false;
            end% if
            eoEpoch = 0;
            
            %search data
            while (eoEpoch==0)
                %read the string
                line = line+1; lin = txtRin{line};
                %answer = strfind(lin,'COMMENT');
                keywords = {'COMMENT', 'MARKER NAME', 'MARKER NUMBER', 'APPROX POSITION XYZ', 'ANTENNA: DELTA H/E/N'};
                answer = [];
                s = 1;
                while (s <= length(keywords) && isempty(answer))
                    answer = strfind(lin,keywords{s});
                    s = s + 1;
                end
                %if it is a line that should be skipped read the following one
                while (~isempty(answer) && (line<length(txtRin)))
                    line = line+1; lin = txtRin{line};
                    %check again
                    answer = [];
                    s = 1;
                    while (s <= length(keywords) && isempty(answer))
                        answer = strfind(lin,keywords{s});
                        s = s + 1;
                    end
                end
                %check if the end of file is reached
                if (line==length(txtRin));
                    return
                end
                
                %check RINEX version
                if (~strcmp(lin(1),'>')) %RINEX v2.xx
                    %check if it is a string that should be analyzed
                    if (strcmp(lin(29),'0') || strcmp(lin(29),'1') || strcmp(lin(29),'2'))
                        
                        %save time information
                        data   = textscan(lin(1:26),'%f%f%f%f%f%f');
                        year   = data{1};
                        month  = data{2};
                        day    = data{3};
                        hour   = data{4};
                        minute = data{5};
                        second = data{6};
                        
                        %computation of the GPS time in weeks and seconds of week
                        year = four_digit_year(year);
                        [week, time] = date2gps([year, month, day, hour, minute, second]); %#ok<ASGLU>
                        
                        %number of visible satellites
                        [num_sat] = sscanf(lin(30:32),'%d');
                        
                        %keep just the satellite data
                        lin = ExtractSubstring(lin, 33, 68);
                        
                        %remove 'blank spaces' and unwanted characters at the end of the string
                        lin = RemoveUnwantedTrailingSpaces(lin);
                        
                        %read additional lines, depending on the number of satellites
                        nlines = ceil(num_sat/12);
                        for n = 1 : nlines - 1
                            line = line+1;
                            lin = [lin ExtractSubstring(txtRin{line}, 33, 68)];
                            lin = RemoveUnwantedTrailingSpaces(lin);
                        end
                        
                        pos = 1;
                        sat = zeros(num_sat,1);
                        sat_types = char(32*uint8(ones(num_sat,1))');
                        for i = 1 : num_sat
                            %check if GPS satellites are labeled 'G' or not labeled
                            if (strcmp(lin(pos),' '))
                                type = 'G';
                            else
                                type = lin(pos);
                            end
                            % sat_types = [sat_types; type];
                            sat_types(i) = type;
                            % sat(i) = sscanf(lin(pos+1:pos+2),'%d');
                            sat(i) = mod((lin(pos+1)-48)*10+(lin(pos+2)-48),160);
                            pos = pos + 3;
                        end
                        
                        eoEpoch = 1;
                    end
                    
                else %RINEX v3.xx
                    
                    %check if it is a string that should be analyzed
                    if (strcmp(lin(29),'0') || strcmp(lin(29),'1') || strcmp(lin(29),'2'))
                        
                        %save time information
                        data   = textscan(lin(2:29),'%f%f%f%f%f%f');
                        year   = data{1};
                        month  = data{2};
                        day    = data{3};
                        hour   = data{4};
                        minute = data{5};
                        second = data{6};
                        
                        %computation of the GPS time in weeks and seconds of week
                        [week, time] = date2gps([year, month, day, hour, minute, second]); %#ok<ASGLU>
                        
                        %number of visible satellites
                        [num_sat] = sscanf(lin(33:35),'%d');
                        
                        eoEpoch = 1;
                    end
                end
            end
            
            if datee_RequestedInOutputFlag
                datee = [year month day hour minute second];
            end %if
        end
        
        function [obs_struct, line] = RinGetObs(txtRin, line, nSat, sat, sat_types, obs_col, nObsTypes, constellations)
            %total number of satellites (according to enabled constellations)
            nSatTot = constellations.nEnabledSat;
            
            %array to contain the starting index for each constellation in the total array (with length = nSatTot)
            sat_types_id = zeros(size(sat_types));
            
            %starting index in the total array for the various constellations
            idGPS = constellations.GPS.indexes(1);
            idGLONASS = constellations.GLONASS.indexes(1);
            idGalileo = constellations.Galileo.indexes(1);
            idBeiDou = constellations.BeiDou.indexes(1);
            idQZSS = constellations.QZSS.indexes(1);
            idSBAS = constellations.SBAS.indexes(1);
            
            %output observations structure initialization
            obs_struct = struct('L1', zeros(nSatTot,1), 'L2', zeros(nSatTot,1), ...
                'C1', zeros(nSatTot,1), ...
                'P1', zeros(nSatTot,1), 'P2', zeros(nSatTot,1), ...
                'S1', zeros(nSatTot,1), 'S2', zeros(nSatTot,1), ...
                'D1', zeros(nSatTot,1), 'D2', zeros(nSatTot,1));
            
            obs_tmp = struct('TMP1', zeros(nSatTot,1), 'TMP2', zeros(nSatTot,1));
            
            if (~isempty(sat_types)) %RINEX v2.xx
                
                %convert constellations letter to starting index in the total array
                sat_types_id(sat_types == 'G') = idGPS*constellations.GPS.enabled;
                sat_types_id(sat_types == 'R') = idGLONASS*constellations.GLONASS.enabled;
                sat_types_id(sat_types == 'E') = idGalileo*constellations.Galileo.enabled;
                sat_types_id(sat_types == 'C') = idBeiDou*constellations.BeiDou.enabled;
                sat_types_id(sat_types == 'J') = idQZSS*constellations.QZSS.enabled;
                sat_types_id(sat_types == 'S') = idSBAS*constellations.SBAS.enabled;
                
                %observation types
                nLinesToRead = ceil(nObsTypes/5);  % I read a maximum of 5 obs per line => this is the number of lines to read
                nObsToRead = nLinesToRead * 5;     % Each line contains 5 observations
                
                %data read and assignment
                lin = char(32*uint8(ones(16*nObsToRead,1))'); % preallocate the line to read
                
                % Mask to filter all the possible observations (max 15)
                mask = false(16,nObsToRead);
                mask(2:14,:) = true;
                % preallocate a matrix of 15 strings (of length 14 characters)
                % notice that each observation element has a max length of 13 char,
                % the first character is added as a padding to separate the strings for
                % the command sscanf that now can be launched once for each satellite
                strObs = char(ones(14,nObsToRead)*32);
                
                for s = 1 : nSat
                    
                    %DEBUG
                    if (sat(s) > 32)
                        sat(s) = 32; %this happens only with SBAS; it's already fixed in the multi-constellation version
                    end
                    
                    lin = char(lin*0+32); % clear line -> fill with spaces
                    
                    if (sat_types_id(s) ~= 0)
                        % read all the lines containing the observations needed
                        for l = 1 : (nLinesToRead)
                            line = line+1; linTmp = txtRin{line};
                            linLengthTmp = length(linTmp);
                            lin((80*(l-1))+(1:linLengthTmp)) = linTmp;  %each line has a maximum lenght of 80 characters
                        end
                        linLength = 80*(nLinesToRead-1)+linLengthTmp;
                        
                        % convert the lines read from the RINEX file to a single matrix
                        % containing all the observations
                        strObs(1:13,:) = (reshape(lin(mask(:)),13,nObsToRead));
                        fltObs = sscanf(strObs, '%f'); % read all the observations in the string
                        obsId = 0; % index of the current observation
                        % start parsing the observation string
                        for k = 1 : min(nObsTypes, floor(linLength/16))
                            % check if the element is empty
                            %if (~strcmp(strObs(:,k)','              ')) % if the current val is not missing (full of spaces)
                            if (sum(strObs(:,k))~=448) % if the current val is not missing (full of spaces)
                                obsId = obsId+1;
                                %obs = sscanf(lin(mask(:,k)), '%f');
                                obs = fltObs(obsId);
                                
                                %check and assign the observation type
                                if (any(~(k-obs_col.L1)))
                                    obs_struct.L1(sat_types_id(s)+sat(s)-1) = obs;
                                    if (linLength>=16*k)
                                        %convert signal-to-noise ratio
                                        % faster conversion of a single ASCII character into an int
                                        snr = mod((lin(16*k)-48),16);
                                        obs_tmp.TMP1(sat_types_id(s)+sat(s)-1) = 6 * snr;
                                    end
                                elseif (any(~(k-obs_col.L2)))
                                    obs_struct.L2(sat_types_id(s)+sat(s)-1) = obs;
                                    if (linLength>=16*k)
                                        %convert signal-to-noise ratio
                                        % faster conversion of a single ASCII character into an int
                                        snr = mod((lin(16*k)-48),16);
                                        obs_tmp.TMP2(sat_types_id(s)+sat(s)-1) = 6 * snr;
                                    end
                                elseif (any(~(k-obs_col.C1)))
                                    obs_struct.C1(sat_types_id(s)+sat(s)-1) = obs;
                                elseif (any(~(k-obs_col.P1)))
                                    obs_struct.P1(sat_types_id(s)+sat(s)-1) = obs;
                                elseif (any(~(k-obs_col.P2)))
                                    obs_struct.P2(sat_types_id(s)+sat(s)-1) = obs;
                                elseif (any(~(k-obs_col.D1)))
                                    obs_struct.D1(sat_types_id(s)+sat(s)-1) = obs;
                                elseif (any(~(k-obs_col.D2)))
                                    obs_struct.D2(sat_types_id(s)+sat(s)-1) = obs;
                                elseif (any(~(k-obs_col.S1)))
                                    obs_struct.S1(sat_types_id(s)+sat(s)-1) = obs;
                                elseif (any(~(k-obs_col.S2)))
                                    obs_struct.S2(sat_types_id(s)+sat(s)-1) = obs;
                                end
                            end
                        end
                        if (~obs_struct.S1(sat_types_id(s)+sat(s)-1))
                            obs_struct.S1(sat_types_id(s)+sat(s)-1) = obs_tmp.TMP1(sat_types_id(s)+sat(s)-1);
                        end
                        if (~obs_struct.S2(sat_types_id(s)+sat(s)-1))
                            obs_struct.S2(sat_types_id(s)+sat(s)-1) = obs_tmp.TMP2(sat_types_id(s)+sat(s)-1);
                        end
                    else
                        %skip all the observation lines for the unused satellite
                        for l = 1 : (nLinesToRead)
                            line = line+1; linTmp = txtRin{line};
                        end
                    end
                end
                
            else %RINEX v3.xx
                
                for s = 1 : nSat
                    
                    %read the line for satellite 's'
                    line = line+1; linTmp = txtRin{line};
                    
                    %read the constellation ID
                    sysId = linTmp(1);
                    
                    %read the satellite PRN/slot number
                    satId = str2num(linTmp(2:3));
                    
                    %number of observations to be read on this line
                    nObsToRead = nObsTypes.(sysId);
                    
                    %data read and assignment
                    lin = char(32*uint8(ones(16*nObsToRead,1))'); % preallocate the line to read
                    
                    %keep only the part of 'lin' containing the observations
                    linTmp = linTmp(4:end);
                    linLengthTmp = length(linTmp);
                    
                    %fill in the 'lin' variable with the actual line read (may be shorter than expected)
                    lin(1:linLengthTmp) = linTmp;
                    linLength = length(lin);
                    
                    %compute the index in the total array and check if the current
                    %constellation is required (if not, skip the line)
                    switch (sysId)
                        case 'G'
                            if (constellations.GPS.enabled)
                                index = idGPS + satId - 1;
                            else
                                continue
                            end
                        case 'R'
                            if (constellations.GLONASS.enabled)
                                index = idGLONASS + satId - 1;
                            else
                                continue
                            end
                        case 'E'
                            if (constellations.Galileo.enabled)
                                index = idGalileo + satId - 1;
                            else
                                continue
                            end
                        case 'C'
                            if (constellations.BeiDou.enabled)
                                index = idBeiDou + satId - 1;
                            else
                                continue
                            end
                        case 'J'
                            if (constellations.QZSS.enabled)
                                index = idQZSS + satId - 1;
                            else
                                continue
                            end
                        case 'S'
                            if (constellations.SBAS.enabled)
                                index = idSBAS + satId - 1;
                            else
                                continue
                            end
                    end
                    
                    % Mask to filter all the possible observations
                    mask = false(16,nObsToRead);
                    mask(2:14,:) = true;
                    % preallocate a matrix of n strings (of length 14 characters)
                    % notice that each observation element has a max length of 13 char,
                    % the first character is added as a padding to separate the strings for
                    % the command sscanf that now can be launched once for each satellite
                    strObs = char(ones(14,nObsToRead)*32);
                    
                    % convert the lines read from the RINEX file to a single matrix
                    % containing all the observations
                    strObs(1:13,:) = (reshape(lin(mask(:)),13,nObsToRead));
                    fltObs = sscanf(strObs, '%f'); % read all the observations in the string
                    obsId = 0; % index of the current observation
                    % start parsing the observation string
                    for k = 1 :  min(nObsToRead, floor(linLength/16))
                        % check if the element is empty
                        if (~strcmp(strObs(:,k)','              ')) % if the current val is not missing (full of spaces)
                            obsId = obsId+1;
                            %obs = sscanf(lin(mask(:,k)), '%f');
                            obs = fltObs(obsId);
                            
                            %check and assign the observation type
                            if (any(~(k-obs_col.(sysId).L1)))
                                obs_struct.L1(index) = obs;
                                if (linLength>=16*k)
                                    %convert signal-to-noise ratio
                                    % faster conversion of a single ASCII character into an int
                                    snr = mod((lin(16*k)-48),16);
                                    obs_tmp.TMP1(index) = 6 * snr;
                                end
                            elseif (any(~(k-obs_col.(sysId).L2)))
                                obs_struct.L2(index) = obs;
                                if (linLength>=16*k)
                                    %convert signal-to-noise ratio
                                    % faster conversion of a single ASCII character into an int
                                    snr = mod((lin(16*k)-48),16);
                                    obs_tmp.TMP2(index) = 6 * snr;
                                end
                            elseif (any(~(k-obs_col.(sysId).C1)))
                                obs_struct.C1(index) = obs;
                            elseif (any(~(k-obs_col.(sysId).P1)))
                                obs_struct.P1(index) = obs;
                            elseif (any(~(k-obs_col.(sysId).P2)))
                                obs_struct.P2(index) = obs;
                            elseif (any(~(k-obs_col.(sysId).D1)))
                                obs_struct.D1(index) = obs;
                            elseif (any(~(k-obs_col.(sysId).D2)))
                                obs_struct.D2(index) = obs;
                            elseif (any(~(k-obs_col.(sysId).S1)))
                                obs_struct.S1(index) = obs;
                            elseif (any(~(k-obs_col.(sysId).S2)))
                                obs_struct.S2(index) = obs;
                            end
                        end
                    end
                    if (~obs_struct.S1(index))
                        obs_struct.S1(index) = obs_tmp.TMP1(index);
                    end
                    if (~obs_struct.S2(index))
                        obs_struct.S2(index) = obs_tmp.TMP2(index);
                    end
                end
            end
            
            clear obs_tmp
        end
    end
    
end
