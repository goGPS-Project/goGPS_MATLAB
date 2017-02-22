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
%   The new specifications are available only in Russian
%   L1 specification: http://russianspacesystems.ru/wp-content/uploads/2016/08/IKD-L1-s-kod.-razd.-Red-1.0-2016.pdf
%   L2 specification: http://russianspacesystems.ru/wp-content/uploads/2016/08/IKD-L2-s-kod.-razd.-Red-1.0-2016.pdf
%   L3 specification: http://russianspacesystems.ru/wp-content/uploads/2016/08/IKD-L3-s-kod.-razd.-Red-1.0-2016.pdf
%
%   Other useful links
%     - http://www.navipedia.net/index.php/GLONASS_Signal_Plan            
%     - http://www.navipedia.net/index.php/Reference_Frames_in_GNSS
%     - http://gage6.upc.es/eknot/Professional_Training/PDF/Reference_Systems.pdf


%--------------------------------------------------------------------------
%               ___ ___ ___ 
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.9.1
% 
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Gatti Andrea
%  Contributors:     Gatti Andrea, ...
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
        f0 = struct('base', struct('R1', 1602.000, 'R2', 1246.000, 'R3', 1201.000), ...
                    'delta',struct('R1', 0.5625, 'R2', 0.4375, 'R3', 0.4375), ...
                    'R1_channels', -7:1:6, ...
                    'R2_channels', -7:1:6, ...
                    'R3_channels', -7:1:8);
                
        % GLONASS (PZ-90.02) Ellipsoid semi-major axis [m]
        ell_a = 6378136;
        % GLONASS (PZ-90.02) Ellipsoid flattening
        ell_f = 1/298.25784;
        % GLONASS (PZ-90.02) Ellipsoid Eccentricity^2
        ell_e2 = (1 - (1 - GLONASS_SS.ell_f) ^ 2);
        % GLONASS (PZ-90.02) Ellipsoid Eccentricity
        ell_e = sqrt(GLONASS_SS.ell_e2);
    end
    
    properties (Constant, Access = 'public')
        % System frequencies as struct [MHz] 
        f = struct('base', GLONASS_SS.f0.base, ...
                   'delta', GLONASS_SS.f0.delta, ...
                   'R1_channels', GLONASS_SS.f0.R1_channels, ...
                   'R2_channels', GLONASS_SS.f0.R2_channels, ...
                   'R3_channels', GLONASS_SS.f0.R3_channels, ...
                   'R1', GLONASS_SS.f0.R1_channels' .* GLONASS_SS.f0.delta.R1 + GLONASS_SS.f0.base.R1, ...
                   'R2', GLONASS_SS.f0.R2_channels' .* GLONASS_SS.f0.delta.R2 + GLONASS_SS.f0.base.R2, ...
                   'R3', GLONASS_SS.f0.R3_channels' .* GLONASS_SS.f0.delta.R3 + GLONASS_SS.f0.base.R3);
        
        % Array of supported frequencies [MHz]
        f_vec = [[GLONASS_SS.f.R1; 0; 0], [GLONASS_SS.f.R2; 0; 0], GLONASS_SS.f.R3] * 1e6;
        
        % Array of the corresponding wavelength - lambda => wavelengths
        l_vec = 299792458 ./ GLONASS_SS.f_vec;   
        
        char_id = 'R'     % Satellite system (ss) character id
        n_sat = 24;       % Maximum number of satellite in the constellation
        prn = (1 : 24)';  % Satellites id numbers as defined in the constellation
    end
        
    properties (Constant, Access = 'public')
        % Structure of orbital parameters (ellipsoid, GM, OMEGA_EARTH_DOT)
        orbital_parameters = struct('GM', 3.986004418e14, ...               % GLONASS (PZ-90.02) Gravitational constant * (mass of Earth) [m^3/s^2]
                                    'OMEGAE_DOT', 7.292115e-5, ...          % GLONASS (PZ-90.02) Angular velocity of the Earth rotation [rad/s]
                                    'ell',struct( ...                       % Ellipsoidal parameters GLONASS (PZ-90.02)
                                        'a', GLONASS_SS.ell_a, ...          % Ellipsoid semi-major axis [m]
                                        'f', GLONASS_SS.ell_f, ...          % Ellipsoid flattening
                                        'e', GLONASS_SS.ell_e, ...          % Eccentricity
                                        'e2', GLONASS_SS.ell_e2));          % Eccentricity^2
                                    
        % GLONASS second zonal harmonic of the geopotential
        J2 = 1.0826257e-3; 
    end
                
    methods
        function this = GLONASS_SS(offset)
            % Creator            
            if (nargin == 0)
                offset = 0;
            end            
            this.updateGoIds(offset);
        end
    end
    
    methods (Access = public)
        function iono_free = getIonoFree(this, f_id, channel)
            % Init iono parameters: iono-free combination is computed with the first two carriers in f_vec (use f_id to change the frequencies to use)
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
            
            iono_free.alpha1 = this.f_vec(channel(1), f_id(1)) .^ 2 ./ (this.f_vec(channel(1), f_id(1)) .^ 2 - this.f_vec(channel(2), f_id(2)) .^ 2);
            iono_free.alpha2 = this.f_vec(channel(2), f_id(2)) .^ 2 ./ (this.f_vec(channel(1), f_id(1)) .^ 2 - this.f_vec(channel(2), f_id(2)) .^ 2);
            gcd_f = gcd(this.f_vec(channel(1), f_id(1)),this.f_vec(channel(2), f_id(2)));
            iono_free.T = this.f_vec(channel(1), f_id(1))/gcd_f;
            iono_free.N = this.f_vec(channel(2), f_id(2))/gcd_f;
        end        
    end
end
