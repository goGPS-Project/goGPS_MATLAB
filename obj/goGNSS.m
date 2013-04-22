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
% Copyright (C) 2009-2013 Mirko Reguzzoni, Eugenio Realini
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
    % => to discriminate them from function (in autocompletition) they are
    % written in capital letters
    properties (Constant)
        % GENERIC CONSTANTS -----------------------------------------------
        
        V_LIGHT = 299792458;                  % velocity of light in the void [m/s]
        
        MAX_SAT = 32                                  % Maximum number of active satellites in a constellation
        
        % CONSTELLATION REF -----------------------------------------------
        % CRS parameters, according to each GNSS system CRS definition
        % (ICD document in brackets):
        %
        % *_GPS --> WGS-84  (IS-GPS200E)
        % *_GLO --> PZ-90   (GLONASS-ICD 5.1)
        % *_GAL --> GTRF    (Galileo-ICD 1.1)
        % *_BDS --> CSG2000 (BeiDou-ICD 1.0)
        % *_QZS --> WGS-84  (IS-QZSS 1.5D)
        
        ELL_A_GPS = 6378137;                          % GPS (WGS-84)     Ellipsoid semi-major axis [m]
        ELL_A_GLO = 6378136;                          % GLONASS (PZ-90)  Ellipsoid semi-major axis [m]
        ELL_A_GAL = 6378137;                          % Galileo (GTRF)   Ellipsoid semi-major axis [m]
        ELL_A_BDS = 6378136;                          % BeiDou (CSG2000) Ellipsoid semi-major axis [m]
        ELL_A_QZS = 6378137;                          % QZSS (WGS-84)    Ellipsoid semi-major axis [m]
        
        ELL_F_GPS = 1/298.257222101;                  % GPS (WGS-84)     Ellipsoid flattening
        ELL_F_GLO = 1/298.257222101;                  % GLONASS (PZ-90)  Ellipsoid flattening
        ELL_F_GAL = 1/298.257222101;                  % Galileo (GTRF)   Ellipsoid flattening
        ELL_F_BDS = 1/298.257222101;                  % BeiDou (CSG2000) Ellipsoid flattening
        ELL_F_QZS = 1/298.257222101;                  % QZSS (WGS-84)    Ellipsoid flattening
        
        ELL_E_GPS = sqrt(1-(1-goGNSS.ELL_F_GPS)^2);   % GPS (WGS-84)     Eccentricity
        ELL_E_GLO = sqrt(1-(1-goGNSS.ELL_F_GLO)^2);   % GLONASS (PZ-90)  Eccentricity
        ELL_E_GAL = sqrt(1-(1-goGNSS.ELL_F_GAL)^2);   % Galileo (GTRF)   Eccentricity
        ELL_E_BDS = sqrt(1-(1-goGNSS.ELL_F_BDS)^2);   % BeiDou (CSG2000) Eccentricity
        ELL_E_QZS = sqrt(1-(1-goGNSS.ELL_F_QZS)^2);   % QZSS (WGS-84)    Eccentricity
        
        GM_GPS = 3.986005e14;                     % GPS     Gravitational constant * (mass of Earth) [m^3/s^2]
        GM_GLO = 3.9860044e14;                    % GLONASS Gravitational constant * (mass of Earth) [m^3/s^2]
        GM_GAL = 3.986004418e14;                  % Galileo Gravitational constant * (mass of Earth) [m^3/s^2]
        GM_BDS = 3.986004418e14;                  % BeiDou  Gravitational constant * (mass of Earth) [m^3/s^2]
        GM_QZS = 3.986005e14;                     % QZSS    Gravitational constant * (mass of Earth) [m^3/s^2]
        
        OMEGAE_DOT_GPS = 7.2921151467e-5;             % GPS     Angular velocity of the Earth rotation [rad/s]
        OMEGAE_DOT_GLO = 7.292115e-5;                 % GLONASS Angular velocity of the Earth rotation [rad/s]
        OMEGAE_DOT_GAL = 7.2921151467e-5;             % Galileo Angular velocity of the Earth rotation [rad/s]
        OMEGAE_DOT_BDS = 7.292115e-5;                 % BeiDou  Angular velocity of the Earth rotation [rad/s]
        OMEGAE_DOT_QZS = 7.2921151467e-5;             % QZSS    Angular velocity of the Earth rotation [rad/s]
        
        J2_GLO = 1.0826257e-3;                        % GLONASS second zonal harmonic of the geopotential
        
        PI_ORBIT = 3.1415926535898;                   % pi value used for orbit computation
        CIRCLE_RAD = 2*goGNSS.PI_ORBIT;               % 2 pi
        
        % CONSTELLATION SPECIFIC ------------------------------------------
        
        FG = [1575.420 1227.600]*1e6;                 % GPS carriers frequencies [Hz]
        LAMBDAG = goGNSS.V_LIGHT ./ goGNSS.FG;        % GPS carriers wavelengths [m]
        
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
        MODE_PP_LS_CP_SA     = 3;   % Post Proc Least Squares on Code Stand Alone (BASE FOR FUTURE PPP IMPLEMENTATION)
        MODE_PP_LS_CP_VEL    = 3.1; % Post Proc Least Squares on Code and Phase for Velocity estimation
        MODE_PP_LS_C_DD      = 11;  % Post Proc Least Squares on Code Double Differences
        MODE_PP_LS_CP_DD_L   = 13;  % Post Proc Least Squares on Code Double Differences with Lambda
        
        MODE_PP_KF_C_SA      = 2;   % Post Proc Kalman Filter on Code Stand Alone
        MODE_PP_KF_C_DD      = 12;  % Post Proc Kalman Filter on Code Double Differencies
        MODE_PP_KF_CP_SA     = 4;   % Post Proc Kalman Filter on Code and Phase Stand Alone
        MODE_PP_KF_CP_DD     = 14;  % Post Proc Kalman Filter on Code and Phase Double Differences
        MODE_PP_KF_CP_DD_MR  = 15;  % Post Proc Kalman Filter on Code and Phase Double Differences Multiple Receivers
        
        GMODE_PP = [ goGNSS.MODE_PP_LS_C_SA ...     % Group of post processing modes
            goGNSS.MODE_PP_LS_CP_SA ...
            goGNSS.MODE_PP_LS_C_DD ...
            goGNSS.MODE_PP_LS_CP_DD_L ...
            goGNSS.MODE_PP_LS_CP_VEL ...
            goGNSS.MODE_PP_KF_C_SA ...
            goGNSS.MODE_PP_KF_C_DD ...
            goGNSS.MODE_PP_KF_CP_SA ...
            goGNSS.MODE_PP_KF_CP_DD ...
            goGNSS.MODE_PP_KF_CP_DD_MR];
        
        GMODE_RT = [ goGNSS.MODE_RT_NAV ...         % Group of real time modes
            goGNSS.MODE_RT_R_MON ...
            goGNSS.MODE_RT_M_MON ...
            goGNSS.MODE_RT_RM_MON];
        
        GMODE_SA = [ goGNSS.MODE_PP_LS_C_SA ...     % Group of stand alone modes
            goGNSS.MODE_PP_LS_CP_SA ...
            goGNSS.MODE_PP_LS_CP_VEL ...
            goGNSS.MODE_PP_KF_C_SA ...
            goGNSS.MODE_PP_KF_CP_SA];
        
        GMODE_DD = [ goGNSS.MODE_PP_LS_C_DD ...     % Group of double differences modes
            goGNSS.MODE_PP_LS_CP_DD_L ...
            goGNSS.MODE_PP_KF_C_DD ...
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
            % return wheather or not the mode given in use is a Post Processing mode
            isPostProcessing = sum(intersect(mode, goGNSS.GMODE_PP));
        end
        
        function isRealTime = isRT(mode)
            % return wheather or not the mode given in use is a Real Time mode
            isRealTime = sum(intersect(mode, goGNSS.GMODE_RT));
        end
        
        function isDoubleDifferences = isDD(mode)
            % return wheather or not the mode given in use is a Double Difference mode
            isDoubleDifferences = sum(intersect(mode, goGNSS.GMODE_DD));
        end
        
        function isStandAlone = isSA(mode)
            % return wheather or not the mode given in use is a Stand Alone mode
            isStandAlone = sum(intersect(mode, goGNSS.GMODE_SA));
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
    end
    
end
