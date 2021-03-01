%   CLASS QZSS_SS
% =========================================================================
%
% DESCRIPTION
%   container of QZSS Satellite System parameters
%
% REFERENCES
%   CRS parameters, according to each GNSS system CRS definition
%   (ICD document in brackets):
%
%   *_QZS --> GRS80    (IS-QZSS 1.8E)
%   Standard: http://qz-vision.jaxa.jp/USE/is-qzss/DOCS/IS-QZSS_18_E.pdf
%
%   Other useful links
%     - http://www.navipedia.net/index.php/QZSS_Signal_Plan
%     - Ellipsoid: http://www.unoosa.org/pdf/icg/2012/template/WGS_84.pdf
%       note that GM and OMEGAE_DOT are redefined in the standard IS-GPS200H (GPS is not using WGS_84 values)
%     - http://www.navipedia.net/index.php/Reference_Frames_in_GNSS
%     - http://gage6.upc.es/eknot/Professional_Training/PDF/Reference_Systems.pdf

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 Geomatics Research & Development srl (GReD)
%  Written by:        Andrea Gatti
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

classdef QZSS_SS < Satellite_System
    properties (Constant, Access = 'public')
        SYS_EXT_NAME = 'QZSS';    % full name of the constellation
        SYS_NAME     = 'QZS';     % 3 characters name of the constellation, this "short name" is used as fields of the property list (struct) to identify a constellation
        SYS_C        = 'J';       % Satellite system (ss) character id

        % System frequencies as struct [MHz]
        F = struct('L1', 1575.420, ...
                   'L2', 1227.600, ...
                   'L5', 1176.450, ...
                   'LEX6', 1278.750)

        % Array of supported frequencies [MHz]
        F_VEC = struct2array(QZSS_SS.F) * 1e6;

        % Array of the corresponding wavelength - lambda => wavelengths
        L_VEC = 299792458 ./ QZSS_SS.F_VEC;

        N_SAT = 9;           % Maximum number of satellite in the constellation
        PRN = (1 : 9)';      % Satellites id numbers as defined in the constellation

        % CODE2DATA ftp://igs.org/pub/data/format/rinex303.pdf
        CODE_RIN3_ATTRIB  = {'ZXLSC F' 'XLS F', 'XIQ F', 'XLS F'}; % last letter of the observation code
        CODE_RIN3_DEFAULT_ATTRIB  = {'C' 'S' 'Q' 'S'}; % last letter of the observation code
        CODE_RIN3_2BAND  = '1256';                        % id for the freq as stored in F_VEC
        IONO_FREE_PREF  = ['12';'15';'16';'25';'26';'56']; % to be evaluated which combination is really better
    end

    properties (Constant, Access = 'private')
        % QZSS (GRS80) Ellipsoid semi-major axis [m]
        ELL_A = 6378137;
        % QZSS (GRS80) Ellipsoid flattening
        ELL_F = 1/298.257222101;
        % QZSS (GRS80) Ellipsoid Eccentricity^2
        ELL_E2 = (1 - (1 - QZSS_SS.ELL_F) ^ 2);
        % QZSS (GRS80) Ellipsoid Eccentricity
        ELL_E = sqrt(QZSS_SS.ELL_E2);
    end

    properties (Constant, Access = 'public')
        % Structure of orbital parameters (ellipsoid, GM, OMEGA_EARTH_DOT)
        ORBITAL_P = struct('GM', 3.986005e14, ...                  % Gravitational constant * (mass of Earth) [m^3/s^2]
                                    'OMEGAE_DOT', 7.2921151467e-5, ...      % Angular velocity of the Earth rotation [rad/s]
                                    'ELL',struct( ...                       % Ellipsoidal parameters QZSS (GRS80)
                                    'A', QZSS_SS.ELL_A, ...             % Ellipsoid semi-major axis [m]
                                    'F', QZSS_SS.ELL_F, ...             % Ellipsoid flattening
                                    'E', QZSS_SS.ELL_E, ...             % Eccentricity
                                    'E2', QZSS_SS.ELL_E2));             % Eccentricity^2
        ORBITAL_INC = 48;    % Orbital inclination        
        ORBITAL_RADIUS  = 42164000 + 6378137; % Orbital radius
    end

    methods
        function this = QZSS_SS(offset)
            % Creator
            % SYNTAX: QZSS_SS(<offset>)
            if (nargin == 0)
                offset = 0;
            end
            this@Satellite_System(offset);
        end

        function copy = getCopy(this)
            % Get a copy of this
            copy = QZSS_SS(this.getOffset());
            copy.import(this);
        end
    end
end
