%   CLASS IRNSS_SS
% =========================================================================
%
% DESCRIPTION
%   container of IRNSS Satellite System parameters
%
% REFERENCES
%   CRS parameters, according to each GNSS system CRS definition
%   (ICD document in brackets):
%
%   *_QZS --> GRS80    (IS-IRNSS 1.8E)
%   Standard: https://www.google.it/url?sa=t&rct=j&q=&esrc=s&source=web&cd=2&ved=0ahUKEwjF5LXHyOPWAhUJLZoKHenOD4MQFggvMAE&url=https%3A%2F%2Fforum.nasaspaceflight.com%2Findex.php%3Faction%3Ddlattach%3Btopic%3D36710.0%3Battach%3D634319&usg=AOvVaw0sOE3-CqbC2o1p2tjjO8IM
%
%   Other useful links
%     - http://www.navipedia.net/index.php/IRNSS_Signal_Plan
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

classdef IRNSS_SS < Satellite_System
    properties (Constant, Access = 'public')
        SYS_EXT_NAME = 'IRNSS';    % full name of the constellation
        SYS_NAME     = 'IRN';     % 3 characters name of the constellation, this "short name" is used as fields of the property list (struct) to identify a constellation
        SYS_C        = 'I';       % Satellite system (ss) character id

        % System frequencies as struct [MHz]
        F = struct('L5', 1176.450, ...
                   'S',  2492.028)
        
        % Array of supported frequencies [MHz]
        F_VEC = struct2array(IRNSS_SS.F) * 1e6;

        % Array of the corresponding wavelength - lambda => wavelengths
        L_VEC = 299792458 ./ IRNSS_SS.F_VEC;

        N_SAT = 9;                    % Maximum number of satellite in the constellation
        PRN = (1 : IRNSS_SS.N_SAT)';  % Satellites id numbers as defined in the constellation
        
        % CODE2DATA ftp://igs.org/pub/data/format/rinex303.pdf
        CODE_RIN3_ATTRIB  = {'XCBA F', 'XCBA F'}; % last letter of the observation code e.g. C5A - C5B - C5C - C5X
        CODE_RIN3_DEFAULT_ATTRIB  = {'C' 'C'}; % last letter of the observation code
        CODE_RIN3_2BAND  = '59';             % id for the freq as stored in F_VEC e.g. L5 -> C5A, S -> C9A
        IONO_FREE_PREF  = ['59'];  % to be evaluated which combination is really better
end

    properties (Constant, Access = 'private')
        % IRNSS (GRS80) Ellipsoid semi-major axis [m]
        ELL_A = 6378137;
        % IRNSS (GRS80) Ellipsoid flattening
        ELL_F = 1/298.257222101;
        % IRNSS (GRS80) Ellipsoid Eccentricity^2
        ELL_E2 = (1 - (1 - IRNSS_SS.ELL_F) ^ 2);
        % IRNSS (GRS80) Ellipsoid Eccentricity
        ELL_E = sqrt(IRNSS_SS.ELL_E2);
    end

    properties (Constant, Access = 'public')
        % Structure of orbital parameters (ellipsoid, GM, OMEGA_EARTH_DOT)
        ORBITAL_P = struct('GM', 3.986005e14, ...                  % Gravitational constant * (mass of Earth) [m^3/s^2]
                           'OMEGAE_DOT', 7.2921151467e-5, ...      % Angular velocity of the Earth rotation [rad/s]
                           'ELL',struct( ...                       % Ellipsoidal parameters IRNSS (GRS80)
                           'A', IRNSS_SS.ELL_A, ...             % Ellipsoid semi-major axis [m]
                           'F', IRNSS_SS.ELL_F, ...             % Ellipsoid flattening
                           'E', IRNSS_SS.ELL_E, ...             % Eccentricity
                           'E2', IRNSS_SS.ELL_E2));             % Eccentricity^2
        ORBITAL_INC = 31;    % Orbital inclination        
        ORBITAL_RADIUS  = 35786000 + 6378137; % Orbital radius
    end

    methods
        function this = IRNSS_SS(offset)
            % Creator
            % SYNTAX: IRNSS_SS(<offset>)
            if (nargin == 0)
                offset = 0;
            end
            this@Satellite_System(offset);
        end

        function copy = getCopy(this)
            % Get a copy of this
            copy = IRNSS_SS(this.getOffset());
            copy.import(this);
        end
    end
end
