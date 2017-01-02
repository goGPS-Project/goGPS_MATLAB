%   CLASS GO_settings
% =========================================================================
%
% DESCRIPTION
%   Collector of settings to manage a useful parameters of goGPS
%   This singleton class collects multiple objects containing various
%   parameters
%
% EXAMPLE
%   settings = goGNSS();
%
% FOR A LIST OF CONSTANTs and METHODS use doc goGNSS

%----------------------------------------------------------------------------------------------
%                           goGPS v0.5.9
% Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
% Written by:       Gatti Andrea
% Contributors:     Gatti Andrea, ...
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
classdef GO_settings < handle
    
    properties (Constant)
        V_LIGHT = 299792458;                  % velocity of light in the void [m/s]
        
        MAX_SAT = 32;                         % Maximum number of active satellites in a constellation
        N_SYS = 6;                            % Maximum number of constellations
        
        % Values as defined by standards
        PI_ORBIT = 3.1415926535898;
        CIRCLE_RAD = 6.2831853071796;
        
        % Standard atmosphere - Berg, 1948
        ATM = struct('PRES', 1013.25, ...            % pressure [mbar]
            'STD_TEMP', 291.15, ...         % temperature [K]
            'STD_HUMI', 50.0);              % humidity [%]
        
    end
    
    properties % Public Access
        SS = struct();                       % Satellite System
    end
    
    methods (Access = private)
        % Guard the constructor against external invocation.  We only want
        % to allow a single instance of this class.  See description in
        % Singleton superclass.
        function obj = GO_settings()
            % Initialisation of the variables
            obj.init_SS()
        end
        
        function init_SS_ref(obj)
            % CONSTELLATION REF -----------------------------------------------
            % CRS parameters, according to each GNSS system CRS definition
            % (ICD document in brackets):
            %
            % *_GPS --> WGS-84   (IS-GPS200E)
            % *_GLO --> PZ-90    (GLONASS-ICD 5.1)
            % *_GAL --> GTRF     (Galileo-ICD 1.1)
            % *_BDS --> CGCS2000 (BeiDou-ICD 1.0)
            % *_QZS --> WGS-84   (IS-QZSS 1.5D)
            
            obj.SS.GPS.ell.a = 6378137;                          % GPS (WGS-84)      Ellipsoid semi-major axis [m]
            obj.SS.GLO.ell.a = 6378136;                          % GLONASS (PZ-90)   Ellipsoid semi-major axis [m]
            obj.SS.GAL.ell.a = 6378137;                          % Galileo (GTRF)    Ellipsoid semi-major axis [m]
            obj.SS.BDS.ell.a = 6378136;                          % BeiDou (CGCS2000) Ellipsoid semi-major axis [m]
            obj.SS.QZS.ell.a = 6378137;                          % QZSS (WGS-84)     Ellipsoid semi-major axis [m]
            
            obj.SS.GPS.ell.f = 1/298.257222101;                  % GPS (WGS-84)      Ellipsoid flattening
            obj.SS.GLO.ell.f = 1/298.257222101;                  % GLONASS (PZ-90)   Ellipsoid flattening
            obj.SS.GAL.ell.f = 1/298.257222101;                  % Galileo (GTRF)    Ellipsoid flattening
            obj.SS.BDS.ell.f = 1/298.257222101;                  % BeiDou (CGCS2000) Ellipsoid flattening
            obj.SS.QZS.ell.f = 1/298.257222101;                  % QZSS (WGS-84)     Ellipsoid flattening
            
            % e = sqrt(1 - (1 - f) ^ 2)
            obj.SS.GPS.ell.e = 0.0818191910428158;               % GPS (WGS-84)      Eccentricity
            obj.SS.GLO.ell.e = 0.0818191910428158;               % GLONASS (PZ-90)   Eccentricity
            obj.SS.GAL.ell.e = 0.0818191910428158;               % Galileo (GTRF)    Eccentricity
            obj.SS.BDS.ell.e = 0.0818191910428158;               % BeiDou (CGCS2000) Eccentricity
            obj.SS.QZS.ell.e = 0.0818191910428158;               % QZSS (WGS-84)     Eccentricity
            
            obj.SS.GPS.GM = 3.986005e14;                         % GPS     Gravitational constant * (mass of Earth) [m^3/s^2]
            obj.SS.GLO.GM = 3.9860044e14;                        % GLONASS Gravitational constant * (mass of Earth) [m^3/s^2]
            obj.SS.GAL.GM = 3.986004418e14;                      % Galileo Gravitational constant * (mass of Earth) [m^3/s^2]
            obj.SS.BDS.GM = 3.986004418e14;                      % BeiDou  Gravitational constant * (mass of Earth) [m^3/s^2]
            obj.SS.QZS.GM = 3.986005e14;                         % QZSS    Gravitational constant * (mass of Earth) [m^3/s^2]
            
            obj.SS.GPS.OMEGAE_DOT = 7.2921151467e-5;             % GPS     Angular velocity of the Earth rotation [rad/s]
            obj.SS.GLO.OMEGAE_DOT = 7.292115e-5;                 % GLONASS Angular velocity of the Earth rotation [rad/s]
            obj.SS.GAL.OMEGAE_DOT = 7.2921151467e-5;             % Galileo Angular velocity of the Earth rotation [rad/s]
            obj.SS.BDS.OMEGAE_DOT = 7.292115e-5;                 % BeiDou  Angular velocity of the Earth rotation [rad/s]
            obj.SS.QZS.OMEGAE_DOT = 7.2921151467e-5;             % QZSS    Angular velocity of the Earth rotation [rad/s]
            
            obj.SS.GLO.J2 = 1.0826257e-3;                        % GLONASS second zonal harmonic of the geopotential            
        end
        
        function init_SS_frequencies(obj)
            % GPS Freq [MHz]
            obj.SS.GPS.f.L1 = 1575.420;
            obj.SS.GPS.f.L2 = 1227.600;
            obj.SS.GPS.f.L5 = 1176.450;
            obj.SS.GPS.f_vec = struct2array(obj.SS.GPS.f) * 1e6;                                   % all the frequencies
            obj.SS.GPS.l_vec = obj.V_LIGHT ./ obj.SS.GPS.f_vec;                                    % lambda => wavelengths
            obj.SS.GPS.IF.A1 = obj.SS.GPS.f.L1 ^ 2/ (obj.SS.GPS.f.L1 ^ 2 - obj.SS.GPS.f.L2 ^ 2);   % alpha combination iono-free
            obj.SS.GPS.IF.A2 = obj.SS.GPS.f.L2 ^ 2/ (obj.SS.GPS.f.L1 ^ 2 - obj.SS.GPS.f.L2 ^ 2);   % alpha combination iono-free
            obj.SS.GPS.IF.T = 77; % round(L1/10.23/2)                                              % GPS iono-free parameter
            obj.SS.GPS.IF.N = 60; % round(L2/10.23/2)                                              % GPS iono-free parameter
            
            % GLONASS Freq [MHz]
            obj.SS.GLO.f.base.R1  = 1602.000;
            obj.SS.GLO.f.base.R2  = 1246.000;
            obj.SS.GLO.f.delta.R1 = 0.5625;
            obj.SS.GLO.f.delta.R2 = 0.4375;
            obj.SS.GLO.f.R_channels = 6:-1:-7;
            obj.SS.GLO.f.R1 = obj.SS.GLO.f.R_channels' .* obj.SS.GLO.f.delta.R1 + obj.SS.GLO.f.base.R1;      % all the frequencies
            obj.SS.GLO.f.R2 = obj.SS.GLO.f.R_channels' .* obj.SS.GLO.f.delta.R2 + obj.SS.GLO.f.base.R2;      % all the frequencies
            obj.SS.GLO.f_vec = [obj.SS.GLO.f.R1 obj.SS.GLO.f.R2] * 1e6;                                      % all the frequencies
            obj.SS.GLO.l_vec = obj.V_LIGHT ./ obj.SS.GLO.f_vec;                                          % lambda => wavelengths
            obj.SS.GLO.IF.A1 = obj.SS.GLO.f.R1 .^ 2 ./ (obj.SS.GLO.f.R1 .^ 2 - obj.SS.GLO.f.R2 .^ 2);    % alpha combination iono-free
            obj.SS.GLO.IF.A2 = obj.SS.GLO.f.R2 .^ 2 ./ (obj.SS.GLO.f.R1 .^ 2 - obj.SS.GLO.f.R2 .^ 2);    % alpha combination iono-free
            obj.SS.GLO.IF.T = 9;                                                                         % GLONASS iono-free parameter
            obj.SS.GLO.IF.N = 7;                                                                         % GLONASS iono-free parameter
            
            obj.SS.GAL.f.E1  = 1575.420;                          % Galileo Freq [MHz]
            obj.SS.GAL.f.E5a = 1176.450;
            obj.SS.GAL.f.E5b = 1207.140;
            obj.SS.GAL.f.E5  = 1191.795;
            obj.SS.GAL.f.E6  = 1278.750;
            obj.SS.GAL.f_vec = struct2array(obj.SS.GAL.f);                                             % all the frequencies
            obj.SS.GAL.l_vec = obj.V_LIGHT ./ obj.SS.GAL.f_vec;                                        % lambda => wavelengths
            obj.SS.GAL.IF.A1 = obj.SS.GAL.f.E1 ^ 2/ (obj.SS.GAL.f.E1 ^ 2 - obj.SS.GAL.f.E1 ^ 2);       % alpha combination iono-free
            obj.SS.GAL.IF.A2 = obj.SS.GAL.f.E5a ^ 2/ (obj.SS.GAL.f.E5a ^ 2 - obj.SS.GAL.f.E5a ^ 2);    % alpha combination iono-free
            obj.SS.GAL.IF.T = 154; % round(FE1/10.23);                                                 % GALILEO iono-free parameter
            obj.SS.GAL.IF.N = 115; % round(FE5a/10.23);                                                % GALILEO iono-free parameter
            
            % BeiDou Freq [MHz]
            obj.SS.BDS.f.C2  = 1561.098;
            obj.SS.BDS.f.C5b = 1207.140;
            obj.SS.BDS.f.C6  = 1268.520;
            obj.SS.BDS.f.C1  = 1589.740;
            obj.SS.BDS.f_vec = struct2array(obj.SS.BDS.f);                                           % all the frequencies
            obj.SS.BDS.l_vec = obj.V_LIGHT ./ obj.SS.BDS.f_vec;                                      % lambda => wavelengths
            obj.SS.BDS.IF.A1 = obj.SS.BDS.f.C2 ^ 2/ (obj.SS.BDS.f.C2 ^ 2 - obj.SS.BDS.f.C5b ^ 2);    % alpha combination iono-free
            obj.SS.BDS.IF.A2 = obj.SS.BDS.f.C5b ^ 2/ (obj.SS.BDS.f.C2 ^ 2 - obj.SS.BDS.f.C5b ^ 2);   % alpha combination iono-free
            obj.SS.BDS.IF.T = 763; % round(C2/2.046);                                                % BeiDou iono-free combination parameter
            obj.SS.BDS.IF.N = 590; % round(C5b/2.046);                                               % BeiDou iono-free combination parameter
            
            obj.SS.QZS.f.J1 = 1575.420;                          % QZSS Freq [MHz]
            obj.SS.QZS.f.J2 = 1227.600;
            obj.SS.QZS.f.J5 = 1176.450;
            obj.SS.QZS.f.J6 = 1278.750;
            obj.SS.QZS.f_vec = struct2array(obj.SS.QZS.f);      % all the frequencies
            obj.SS.QZS.l_vec = obj.V_LIGHT ./ obj.SS.QZS.f_vec; % lambda => wavelengths
            obj.SS.QZS.IF.A1 = obj.SS.QZS.f.J1 ^ 2/ (obj.SS.QZS.f.J1 ^ 2 - obj.SS.QZS.f.J2 ^ 2);   % alpha combination iono-free
            obj.SS.QZS.IF.A2 = obj.SS.QZS.f.J2 ^ 2/ (obj.SS.QZS.f.J1 ^ 2 - obj.SS.QZS.f.J2 ^ 2);   % alpha combination iono-free
            obj.SS.QZS.IF.T = 77; % round(J1/10.23/2);                                             % QZSS iono-free combination parameter
            obj.SS.QZS.IF.N = 60; % round(J2/10.23/2);                                             % QZSS iono-free combination parameter
            
            obj.SS.SBS.f.S1 = 1575.420;                         % SBS Freq [MHz]
            obj.SS.SBS.f.S5 = 1176.450;
            obj.SS.SBS.f_vec = struct2array(obj.SS.SBS.f) * 1e6;                                     % all the frequencies
            obj.SS.SBS.l_vec = obj.V_LIGHT ./ obj.SS.SBS.f_vec;                                      % lambda => wavelengths
            obj.SS.SBS.IF.A1 = obj.SS.SBS.f.S1 ^ 2/ (obj.SS.SBS.f.S1 ^ 2 - obj.SS.SBS.f.S5 ^ 2);   % alpha combination iono-free
            obj.SS.SBS.IF.A2 = obj.SS.SBS.f.S5 ^ 2/ (obj.SS.SBS.f.S1 ^ 2 - obj.SS.SBS.f.S5 ^ 2);   % alpha combination iono-free
            obj.SS.SBS.IF.T = 154; % round(S1/10.23)                                                  % SBS iono-free parameter
            obj.SS.SBS.IF.N = 115; % round(S5/10.23)                                                  % SBS iono-free parameter            
        end
        
        function init_SS(obj)             
            % set up satellite system parameters            
            obj.SS.ids = 1 : 6; 
            obj.SS.system_name = {'GPS', 'GLO', 'GAL', 'BDS', 'QZS', 'SBS'};
            obj.SS.sysID = 'GRECJS';
            obj.SS.n_sat = uint8([32 24 30 37 4 0]);
            obj.SS.n_sat_tot = sum(obj.SS.n_sat);
            obj.SS.enabled = logical([1 1 1 1 1 0]);
            
            % duplicate some info in struct (may be useful)
            obj.SS.GPS = struct('id', 1, ...
                                'n_sat', 32, ...
                                'enabled', obj.SS.enabled(1), ...
                                'indexes', [], ...
                                'PRN', 1 : 32, ...
                                'sysID', 'G');
            obj.SS.GLO = struct('id', 2, ...
                                'n_sat', 24, ...
                                'enabled', obj.SS.enabled(2), ...
                                'indexes', [], ...
                                'PRN', [1 : 24], ...
                                'sysID', 'R');
            obj.SS.GAL = struct('id', 3, ...
                                'n_sat', 30, ...
                                'enabled', obj.SS.enabled(3), ...
                                'indexes', [], ...
                                'PRN', 1 : 30, ...
                                'sysID', 'E');
            obj.SS.BDS = struct('id', 4, ...
                                'n_sat', 37, ...
                                'enabled', obj.SS.enabled(4), ...
                                'indexes', [], ...
                                'PRN', 1 : 37, ...
                                'sysID', 'C');
            obj.SS.QZS = struct('id', 5, ...
                                'n_sat', 4, ...
                                'enabled', obj.SS.enabled(5), ...
                                'indexes', [], ...
                                'PRN', 193 : 196, ...
                                'sysID', 'J');
            obj.SS.SBS = struct('id', 6, ...
                                'n_sat', 0, ...
                                'enabled', obj.SS.enabled(6), ...
                                'indexes', [], ...
                                'PRN', [], ...       % not yet enabled for ranging
                                'sysID', 'S');                             
                            
            obj.set_SS(obj.SS.enabled);
            obj.init_SS_ref();
            obj.init_SS_frequencies();
        end
        
    end
    
    methods (Static)
        % Concrete implementation.  See Singleton superclass.
        function obj = get_instance()
            persistent uniqueInstance
            if isempty(uniqueInstance)
                obj = GO_settings();
                uniqueInstance = obj;
            else
                obj = uniqueInstance;
            end
        end
    end
    
    %*** Define your own methods for SingletonImpl.
    methods % Public Access
        function set_SS(obj, enabled_ss)            
            % SYNTAX:
            %   s.set_SS([GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag]);
            %
            % INPUT:
            %   single logical array whose elements are:
            %   GPS_flag = boolean flag for enabling/disabling GPS usage
            %   GLO_flag = boolean flag for enabling/disabling GLONASS usage
            %   GAL_flag = boolean flag for enabling/disabling Galileo usage
            %   BDS_flag = boolean flag for enabling/disabling BeiDou usage
            %   QZS_flag = boolean flag for enabling/disabling QZSS usage
            %   SBS_flag = boolean flag for enabling/disabling SBAS usage (for ranging)
            %
            % OUTPUT:
            %   the results is stored within the object referenced by "s"
            %
            % DESCRIPTION:
            %   Multi-constellation set-up.            
            
            % Check the size of the flag array
            if (nargin == 1)
                enabled_ss = [1 1 1 1 1 0];
            end
            if (numel(enabled_ss) < 6)
                tmp = false(obj.N_SYS, 1);
                tmp(1:numel(enabled_ss)) = enabled_ss;
                enabled_ss = tmp;
                clear tmp;
            end
            obj.SS.enabled = enabled_ss;
            obj.SS.indexes = []; % incremental index of the active satellite system
            obj.SS.PRN = [];     % relative id number in the satellite system
            obj.SS.systems = ''; % char id of the constellation per satellite
                        
            obj.SS.n_sat_tot = 0; % counter for number of satellites
            q = 0;                % counter for enabled constellations
            for s = 1 : numel(obj.SS.system_name)
                if (enabled_ss(s))
                    obj.SS.n_sat_tot = obj.SS.n_sat_tot + obj.SS.(obj.SS.system_name{s}).n_sat;
                    q = q + 1;
                    if (q == 1)
                        ids = [1 : obj.SS.(obj.SS.system_name{s}).n_sat];
                    else
                        ids = [ids(end) + (1 : obj.SS.(obj.SS.system_name{s}).n_sat)];
                    end
                    obj.SS.(obj.SS.system_name{s}).indexes = ids;     
                    obj.SS.(obj.SS.system_name{s}).enabled = true; 
                    obj.SS.indexes = [obj.SS.indexes ids]; % incremental index of the active satellite system
                    obj.SS.PRN = [obj.SS.PRN obj.SS.(obj.SS.system_name{s}).PRN];
                    obj.SS.systems = [obj.SS.systems char(ones(1,obj.SS.(obj.SS.system_name{s}).n_sat) * obj.SS.sysID(s))];
                else
                    obj.SS.(obj.SS.system_name{s}).indexes = []; 
                    obj.SS.(obj.SS.system_name{s}).enabled = false; 
                end
            end
        end                
    end
    
    methods % Public Access (Legacy support)
        function [constellations] = initConstellation(GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag)
            obj.set_SS([GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag])
            constellations = obj.SS;
        end
        
    end
    
end
