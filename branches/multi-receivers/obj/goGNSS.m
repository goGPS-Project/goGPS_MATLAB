% =========================================================================
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
% LIST of CONSTANT
%
%   V_LIGHT                 velocity of light in the void                 [m/s]
%   F1                      GPS carriers frequencies F1                   [1/s]
%   F2                      GPS carriers frequencies F2                   [1/s]
%   ELL_A                   ellipsoid semi-major axis                     [m]
%   ELL_F                   ellipsoid flattening
%   ELL_E                   eccentricity
%
%   OMEGAE_DOT              angular velocity of the Earth rotation        [rad/s]
%   GM                      gravitational constant (mass of Earth)        [m^3/s^2]
%
%   MAX_SAT                 Maximum number of active satellites in a constellation
%
% CONSTELLATIONS IDs ----------------------------------------------
%
%   ID_GPS                  id of GPS observations (for internal use)
%   ID_GLONASS              id of GLONASS observations (for internal use)
%   ID_GALILEO              id of Galileo observations (for internal use)
%   ID_BEIDOU               id of BeiDou observations (for internal use)
%   ID_QZSS                 id of QZSS observations (for internal use)
%   ID_SBAS                 id of SBAS observations (for internal use)
%
% goGPS MODES -----------------------------------------------------
% 
%   MODE_RT_NAV             Real Time Navigation (Kalman Filter on Code and Phase Double Differences (with/without a constraint)
%   MODE_RT_R_MON           Real Time Rover Monitor
%   MODE_RT_M_MON           Real Time Master Monitor
%   MODE_RT_RM_MON          Real Time Master + Rover Monitor
%    
%   MODE_PP_LS_C_SA         Post Proc Least Squares on Code Stand Alone
%   MODE_PP_LS_CP_SA        Post Proc Least Squares on Code Stand Alone (BASE FOR FUTURE PPP IMPLEMENTATION)
%   MODE_PP_LS_C_DD         Post Proc Least Squares on Code Double Differences
%   MODE_PP_LS_CP_DD_L      Post Proc Least Squares on Code Double Differences with Lambda
%   MODE_PP_LS_CP_VEL       Post Proc Least Squares on Code and Phase for Velocity estimation
%    
%   MODE_PP_KF_C_SA         Post Proc Kalman Filter on Code Stand Alone
%   MODE_PP_KF_C_DD         Post Proc Kalman Filter on Code Double Differencies
%   MODE_PP_KF_CP_SA        Post Proc Kalman Filter on Code and Phase Stand Alone
%   MODE_PP_KF_CP_DD        Post Proc Kalman Filter on Code and Phase Double Differences
%   MODE_PP_KF_CP_DD_MR     Post Proc Kalman Filter on Code and Phase Double Differences Multiple Receivers
%
% LIST of METHODS
%
%   [pos] = bancroft(B_pass)
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

classdef goGNSS < handle

    % Constant values 
    % => to discriminate them from function (in autocompletition) they are
    % written in capital letters
    properties (Constant)
        % GENERIC CONSTANTS -----------------------------------------------
                
        V_LIGHT = 299792458;                  % velocity of light in the void [m/s]
                        
        ELL_A = 6378137;                      % Ellipsoid semi-major axis [m]
        ELL_F = 1/298.257222101;              % Ellipsoid flattening                
        ELL_E = sqrt(1-(1-goGNSS.ELL_F)^2);   % Eccentricity
                 
        OMEGAE_DOT = 7.2921151467e-5;         % Angular velocity of the Earth rotation [rad/s]
        GM = 3.986004418e14;                  % Gravitational constant * (mass of Earth) [m^3/s^2]        
        
        MAX_SAT      = 32                     % Maximum number of active satellites in a constellation

        % CONSTELLATION SPECIFIC ------------------------------------------
        
        F1 = 1575.420*1e6;                    % GPS carriers frequency  [Hz]
        F2 = 1227.600*1e6;                    % GPS carriers frequency  [Hz]                
        LAMBDA1 = goGNSS.V_LIGHT / goGNSS.F1; % GPS carriers wavelength [m]
        LAMBDA2 = goGNSS.V_LIGHT / goGNSS.F2; % GPS carriers wavelength [m]
        
        % CONSTELLATIONS IDs ----------------------------------------------
                
        ID_GPS     = 1
        ID_GLONASS = 2
        ID_GALILEO = 3
        ID_BEIDOU  = 4
        ID_QZSS    = 5
        ID_SBAS    = 6
        
        % goGPS MODES -----------------------------------------------------

        MODE_RT_NAV          = 24;  % Real Time Navigation (Kalman Filter on Code and Phase Double Differences (with/without a constraint)
        MODE_RT_R_MON        = 21;  % Real Time Rover Monitor
        MODE_RT_M_MON        = 22;  % Real Time Master Monitor
        MODE_RT_RM_MON       = 23;  % Real Time Master + Rover Monitor
        
        MODE_PP_LS_C_SA      = 1;   % Post Proc Least Squares on Code Stand Alone
        MODE_PP_LS_CP_SA     = 3;   % Post Proc Least Squares on Code Stand Alone (BASE FOR FUTURE PPP IMPLEMENTATION)
        MODE_PP_LS_C_DD      = 11;  % Post Proc Least Squares on Code Double Differences
        MODE_PP_LS_CP_DD_L   = 13;  % Post Proc Least Squares on Code Double Differences with Lambda
        MODE_PP_LS_CP_VEL    = 3.1; % Post Proc Least Squares on Code and Phase for Velocity estimation

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
        function obj = goGNSS()
        end
    end
    
    %   MODE FUNCTION (STATIC)
    % -------------------------------------------------------------------------
    % function to detect a certain kind of processing
    methods (Static, Access = 'public')
        function isPostProcessing = isPP(mode)
           isPostProcessing = sum(intersect(mode, goGNSS.GMODE_PP));
        end
    
        function isRealTime = isRT(mode)
           isRealTime = sum(intersect(mode, goGNSS.GMODE_RT));
        end   
         
        function isDoubleDifferences = isDD(mode)
            isDoubleDifferences = sum(intersect(mode, goGNSS.GMODE_DD));
        end
        
        function isStandAlone = isSA(mode)
            isStandAlone = sum(intersect(mode, goGNSS.GMODE_SA));
        end
    end
    
    %   USEFUL FUNCTION (STATIC)
    % -------------------------------------------------------------------------
    methods (Static, Access = 'public')   
        
        % Function to get a bancroft solution from one receiver
        % it returns pos + dt => X is a [4 x 1] <=> 3 coord + dt
        %
        % Typical old way of calling the function
        % index = find(no_eph == 0);
        % X = getBancroftPos(XS(index,:), dtS(index), pseudorange(index))
        function X = getBancroftPos(XS, dtS, prR)
            matB = [XS(:,:), prR(:) + goGNSS.V_LIGHT * dtS(:)]; % Bancroft matrix
            X = obj.bancroft(matB);            
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
            %----------------------------------------------------------------------------------------------            pos = zeros(4,1);
            
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
