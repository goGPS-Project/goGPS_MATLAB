%   CLASS SBAS_SS
% =========================================================================
%
% DESCRIPTION
%   container of SBAS Satellite System parameters (not yet fully
%   implemented, actually based on GPS)
%
% REFERENCES
%   CRS parameters, according to each GNSS system CRS definition
%       ftp://igs.org/pub/data/format/rinex303.pdf
%   (ICD document in brackets):
%
%   *_SBAS --> WGS-84

%
%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0
%
%--------------------------------------------------------------------------
%  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
%  Written by:        Andrea Gatti, Giulio Tagliaferro ...
%  Contributors:      Andrea Gatti, Giulio Tagliaferro ...
%  A list of all the historical goSBAS contributors is in CREDITS.nfo
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

classdef SBAS_SS < Satellite_System
    properties (Constant, Access = 'public')
        SYS_EXT_NAME = 'SBAS';    % full name of the constellation
        SYS_NAME     = 'SBS';     % 3 characters name of the constellation, this "short name" is used as fields of the property list (struct) to identify a constellation
        SYS_C        = 'S';       % Satellite system (ss) character id

        % System frequencies as struct [MHz]
        F = struct('L1', 1575.420, ...
                   'L5', 1176.450)

        % Array of supported frequencies [MHz]
        F_VEC = struct2array(SBAS_SS.F) * 1e6;

        % Array of the corresponding wavelength - lambda => wavelengths
        L_VEC = 299792458 ./ SBAS_SS.F_VEC;

        N_SAT = 0;           % Maximum number of satellite in the constellation
        PRN = (0 : 0)';      % Satellites id numbers as defined in the constellation

        % CODE2DATA ftp://igs.org/pub/data/format/rinex303.pdf
        CODE_RIN3_ATTRIB  = {'C ' 'XIQ '}; % last letter of the observation code
        CODE_RIN3_DEFAULT_ATTRIB  = {'C' 'Q'}; % last letter of the observation code
        CODE_RIN3_2BAND  = '15';        % id for the freq as stored in F_VEC
        IONO_FREE_PREF  = ['15'];
    end

    properties (Constant, Access = 'private')
        % SBAS (WGS84) Ellipsoid semi-major axis [m]
        ELL_A = 6378137;
        % SBAS (WGS84) Ellipsoid flattening
        ELL_F = 1/298.257223563;
        % SBAS (WGS84) Ellipsoid Eccentricity^2
        ELL_E2 = (1 - (1 - SBAS_SS.ELL_F) ^ 2);
        % SBAS (WGS84) Ellipsoid Eccentricity
        ELL_E = sqrt(SBAS_SS.ELL_E2);
    end

    properties (Constant, Access = 'public')
        % Structure of orbital parameters (ellipsoid, GM, OMEGA_EARTH_DOT)
        ORBITAL_P = struct('GM', 3.986005e14, ...                  % Gravitational constant * (mass of Earth) [m^3/s^2]
                                    'OMEGAE_DOT', 7.2921151467e-5, ...      % Angular velocity of the Earth rotation [rad/s]
                                    'ELL',struct( ...                       % Ellipsoidal parameters SBAS (WGS84)
                                    'A', SBAS_SS.ELL_A, ...              % Ellipsoid semi-major axis [m]
                                    'F', SBAS_SS.ELL_F, ...              % Ellipsoid flattening
                                    'E', SBAS_SS.ELL_E, ...              % Eccentricity
                                    'e2', SBAS_SS.ELL_E2));              % Eccentricity^2
        ORBITAL_INC = 55;    % Orbital inclination        
        ORBITAL_RADIUS  = 20180000 + 6378137; % Orbital radius
    end

    methods
        function this = SBAS_SS(offset)
            % Creator
            % SYNTAX: SBAS_SS(<offset>);
            if (nargin == 0)
                offset = 0;
            end
            this@Satellite_System(offset);
        end

        function copy = getCopy(this)
            % Get a copy of this
            copy = SBAS_SS(this.getOffset());
            copy.import(this);
        end
    end
end
