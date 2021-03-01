%   CLASS BeiDou
% =========================================================================
%
% DESCRIPTION
%   container of BeiDou Satellite System parameters
%
% REFERENCES
%   CRS parameters, according to each GNSS system CRS definition
%   (ICD document in brackets):
%
%   *_BDS --> CGCS2000 (BeiDou-ICD 1.0)
%   Standard: http://www.beidou.gov.cn/attach/2012/12/27/201212273da29c5eb8274deb8cd2b178228ba2bd.pdf
%
%   Other useful links
%     - http://www.navipedia.net/index.php/BeiDou_Signal_Plan
%     - http://www.navipedia.net/index.php/Reference_Frames_in_GNSS
%     - http://gage6.upc.es/eknot/Professional_Training/PDF/Reference_Systems.pdf
%
%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 Geomatics Research & Development srl (GRed)
%  Written by:        Andrea Gatti, Giulio Tagliaferro ...
%  Contributors:      Andrea Gatti, Giulio Tagliaferro ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
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
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

classdef BeiDou_SS < Satellite_System
    properties (Constant, Access = 'public')
        SYS_EXT_NAME = 'BeiDou';  % full name of the constellation
        SYS_NAME     = 'BDS';     % 3 characters name of the constellation, this "short name" is used as fields of the property list (struct) to identify a constellation
        SYS_C        = 'C';       % Satellite system (ss) character id

        % System frequencies as struct [MHz]
        F = struct('B1', 1561.098, ...
                   'B1C', 1575.42, ...
                   'B2a', 1176.45, ...
                   'B2b', 1207.140, ...
                   'B2ab', 1191.795, ...
                   'B3', 1268.520)

        % Array of supported frequencies [MHz]
        F_VEC = struct2array(BeiDou_SS.F) * 1e6;

        % Array of the corresponding wavelength - lambda => wavelengths
        L_VEC = 299792458 ./ BeiDou_SS.F_VEC;

        N_SAT = 37;       % Maximum number of satellite in the constellation
        PRN = (1 : 37)';  % Satellites id numbers as defined in the constellation

        % CODE2DATA ftp://igs.org/pub/data/format/rinex303.pdf
        CODE_RIN3_ATTRIB  = {'XIQ F' 'XPDZLS' 'XPD' 'XIQZPD F', 'ZPD' 'XIQZPD F'}; % last letter of the observation code
        CODE_RIN3_DEFAULT_ATTRIB  = {'Q' 'P' 'Q' 'Q'}; % last letter of the observation code
        CODE_RIN3_2BAND  = '215786';                % id for the freq as stored in F_VEC
        IONO_FREE_PREF  = ['27';'26';'67'];  % to be evaluated which combination is really better
    end

    properties (Constant, Access = 'private')
        % GPS (WGS84) Ellipsoid semi-major axis [m]
        ELL_A = 6378137;
        % GPS (WGS84) Ellipsoid flattening
        ELL_F = 1/298.257222101;
        % GPS (WGS84) Ellipsoid Eccentricity^2
        ELL_E2 = (1 - (1 - BeiDou_SS.ELL_F) ^ 2);
        % GPS (WGS84) Ellipsoid Eccentricity
        ELL_E = sqrt(BeiDou_SS.ELL_E2);
    end

    properties (Constant, Access = 'public')
        % Structure of orbital parameters (ellipsoid, GM, OMEGA_EARTH_DOT)
        ORBITAL_P = struct('GM', 3.986004418e14, ...               % BeiDou (BeiDou-ICD 1.0) Gravitational constant * (mass of Earth) [m^3/s^2]
                                    'OMEGAE_DOT', 7.2921150e-5, ...         % BeiDou (BeiDou-ICD 1.0) Angular velocity of the Earth rotation [rad/s]
                                    'ELL',struct( ...                       % Ellipsoidal parameters BeiDou (CGCS2000)
                                    'A', BeiDou_SS.ELL_A, ...               % Ellipsoid semi-major axis [m]
                                    'F', BeiDou_SS.ELL_F, ...               % Ellipsoid flattening
                                    'E', BeiDou_SS.ELL_E, ...               % Eccentricity
                                    'E2', BeiDou_SS.ELL_E2));               % Eccentricity^2
        ORBITAL_INC = 55;    % Orbital inclination        
        ORBITAL_RADIUS  = 21528000 + 6378137; % Orbital radius
    end

    methods
        function this = BeiDou_SS(offset)
            % Creator
            % SYNTAX: BeiDou_SS(<offset>);
            if (nargin == 0)
                offset = 0;
            end
            this@Satellite_System(offset);
        end

        function copy = getCopy(this)
            % Get a copy of this
            copy = BeiDou_SS(this.getOffset());
            copy.import(this);
        end
    end
end
