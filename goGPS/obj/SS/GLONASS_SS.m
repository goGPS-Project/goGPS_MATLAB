%   CLASS GLONASS_SS
% =========================================================================
%
% DESCRIPTION
%   container of GLONASS Satellite System parameters
%
% REFERENCES
%   CRS parameters, according to each GNSS system CRS definition
%   (ICD document in brackets):
%
%   *_GLO --> PZ-90.02 (GLONASS-ICD 5.1) (GLONASS-ICD 2016)
%   Standard: http://gauss.gge.unb.ca/GLONASS.ICD-98.pdf
%    newer -> http://kb.unavco.org/kb/assets/727/ikd51en.pdf
%
%   Specifications
%   L1   specification: http://russianspacesystems.ru/wp-content/uploads/2016/08/ICD-GLONASS-CDMA-L1.-Edition-1.0-2016.pdf
%   L2   specification: http://russianspacesystems.ru/wp-content/uploads/2016/08/ICD-GLONASS-CDMA-L2.-Edition-1.0-2016.pdf
%   L3   specification: http://russianspacesystems.ru/wp-content/uploads/2016/08/ICD-GLONASS-CDMA-L3.-Edition-1.0-2016.pdf
%   CDMA specification: http://russianspacesystems.ru/wp-content/uploads/2016/08/ICD-GLONASS-CDMA-General.-Edition-1.0-2016.pdf
%
%   Other useful links
%     - http://www.navipedia.net/index.php/GLONASS_Signal_Plan
%     - http://www.navipedia.net/index.php/Reference_Frames_in_GNSS
%     - http://gage6.upc.es/eknot/Professional_Training/PDF/Reference_Systems.pdf
%     - http://gpsworld.com/wp-content/uploads/2017/08/Almanac-GPSWorld-Aug2017.pdf


%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti, Giulio Tagliaferro
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

classdef GLONASS_SS < Satellite_System
    properties (Constant, Access = 'private')
        % System frequencies as struct [MHz]
        f0 = struct('BASE', struct('G1', 1602.000, 'G2', 1246.000, 'G3', 1202.025), ...
                    'DELTA',struct('G1', 0.5625, 'G2', 0.4375, 'G3', 0), ...
                    'G1_CHANNELS', -7:1:6, ...
                    'G2_CHANNELS', -7:1:6, ...
                    'G3_CHANNELS', -7:1:8);

        % GLONASS (PZ-90.02) Ellipsoid semi-major axis [m]
        ELL_A = 6378136;
        % GLONASS (PZ-90.02) Ellipsoid flattening
        ELL_F = 1/298.25784;
        % GLONASS (PZ-90.02) Ellipsoid Eccentricity^2
        ELL_E2 = (1 - (1 - GLONASS_SS.ELL_F) ^ 2);
        % GLONASS (PZ-90.02) Ellipsoid Eccentricity
        ELL_E = sqrt(GLONASS_SS.ELL_E2);
    end

    properties (Constant, Access = 'public')
        SYS_EXT_NAME = 'GLONASS'; % full name of the constellation
        SYS_NAME     = 'GLO';     % 3 characters name of the constellation, this "short name" is used as fields of the property list (struct) to identify a constellation
        SYS_C        = 'R';       % Satellite system (ss) character id

        % System frequencies as struct [MHz]
        F = struct('BASE', GLONASS_SS.f0.BASE, ...
                   'DELTA', GLONASS_SS.f0.DELTA, ...
                   'G1_CHANNELS', GLONASS_SS.f0.G1_CHANNELS, ...
                   'G2_CHANNELS', GLONASS_SS.f0.G2_CHANNELS, ...
                   'G3_CHANNELS', GLONASS_SS.f0.G3_CHANNELS, ...
                   'G1', GLONASS_SS.f0.G1_CHANNELS' .* GLONASS_SS.f0.DELTA.G1 + GLONASS_SS.f0.BASE.G1, ...
                   'G2', GLONASS_SS.f0.G2_CHANNELS' .* GLONASS_SS.f0.DELTA.G2 + GLONASS_SS.f0.BASE.G2, ...
                   'G3', GLONASS_SS.f0.G3_CHANNELS' .* GLONASS_SS.f0.DELTA.G3 + GLONASS_SS.f0.BASE.G3);

        % Array of supported frequencies [MHz]
        F_VEC = [[GLONASS_SS.F.G1; 0; 0], [GLONASS_SS.F.G2; 0; 0], GLONASS_SS.F.G3] * 1e6;

        % Array of the corresponding wavelength - lambda => wavelengths
        L_VEC = 299792458 ./ GLONASS_SS.F_VEC;

        N_SAT = 24;       % Maximum number of satellite in the constellation
        PRN = (1 : 24)';  % Satellites id numbers as defined in the constellation

        % http://gpsworld.com/wp-content/uploads/2017/08/Almanac-GPSWorld-Aug2017.pdf
        % prn to glonass channel
        PRN2CH = [1 -4 5 6 1 -4 5 6 -2 -7 0 -1 -2 -7 0 -1 4 -3 3 2 4 -3 3 2];
        % Same as PRN2IDCH but return the id of the frequency for a certain prn
        PRN2IDCH = GLONASS_SS.PRN2CH + 8;
        
        % CODE2DATA ftp://igs.org/pub/data/format/rinex303.pdf
        CODE_RIN3_ATTRIB  = {'PC F' 'PC F' 'XIQ F'}; % last letter of the observation code
        CODE_RIN3_DEFAULT_ATTRIB  = {'C' 'C' 'Q'}; % last letter of the observation code
        CODE_RIN3_2BAND  = '123';             % id for the freq as stored in F_VEC
        IONO_FREE_PREF  = ['12';'13';'23'];  % to be evaluated which combination is really better
    end

    properties (Constant, Access = 'public')
        % Structure of orbital parameters (ellipsoid, GM, OMEGA_EARTH_DOT)
        ORBITAL_P = struct('GM', 3.986004418e14, ...               % GLONASS (PZ-90.02) Gravitational constant * (mass of Earth) [m^3/s^2]
                                    'OMEGAE_DOT', 7.292115e-5, ...          % GLONASS (PZ-90.02) Angular velocity of the Earth rotation [rad/s]
                                    'ELL',struct( ...                       % Ellipsoidal parameters GLONASS (PZ-90.02)
                                    'A', GLONASS_SS.ELL_A, ...          % Ellipsoid semi-major axis [m]
                                    'F', GLONASS_SS.ELL_F, ...          % Ellipsoid flattening
                                    'E', GLONASS_SS.ELL_E, ...          % Eccentricity
                                    'E2', GLONASS_SS.ELL_E2));          % Eccentricity^2

        % GLONASS second zonal harmonic of the geopotential
        J2 = 1.0826257e-3;
    end

    methods
        function this = GLONASS_SS(offset)
            % Creator
            if (nargin == 0)
                offset = 0;
            end
            this@Satellite_System(offset);
        end

        function copy = getCopy(this)
            % Get a copy of this
            copy = GLONASS_SS(this.getOffset());
        end
    end

    methods (Access = public)
        function iono_free = getIonoFree(this, f_id, channel)
            % Init iono parameters: iono-free combination is computed with the first two carriers in F_VEC (use f_id to change the frequencies to use)
            % SYNTAX: iono_free = glanass_obj.getIonoFree(<f_id>, <channel>)
            % Structure of the iono free combination parameters (alpha 1, alpha 2, T ed N
            % iono_free
            %   .alpha1
            %   .alpha2
            %   .T
            %   .N
            if nargin < 2
                f_id = [1 2];
            end
            if nargin < 3
                channel = 6;
            end
            if isempty(f_id)
                f_id = [1 2];
            end
            if numel(channel) == 1
                channel = [channel channel];
            end

            iono_free.alpha1 = this.F_VEC(channel(1), f_id(1)) .^ 2 ./ (this.F_VEC(channel(1), f_id(1)) .^ 2 - this.F_VEC(channel(2), f_id(2)) .^ 2);
            iono_free.alpha2 = this.F_VEC(channel(2), f_id(2)) .^ 2 ./ (this.F_VEC(channel(1), f_id(1)) .^ 2 - this.F_VEC(channel(2), f_id(2)) .^ 2);
            gcd_f = gcd(this.F_VEC(channel(1), f_id(1)),this.F_VEC(channel(2), f_id(2)));
            iono_free.T = this.F_VEC(channel(1), f_id(1))/gcd_f;
            iono_free.N = this.F_VEC(channel(2), f_id(2))/gcd_f;
        end

        function name = getFreqName(this)
            % Get the name of the frequencies for the current constellation
            % SYNTAX: name = this.getFreqNames();
            name = fieldnames(this.F.BASE);
        end

    end
end
