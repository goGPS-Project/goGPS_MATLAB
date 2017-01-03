%   CLASS GLONASS_SS
% =========================================================================
%
% DESCRIPTION
%   container of GLONASS Satellite System parameters
%
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
classdef GLONASS_SS < Satellite_System
    properties
        indexes       % Satellites unique id numbers in goGPS
        J2            % GLONASS second zonal harmonic of the geopotential
    end
            
    methods
        function obj = GLONASS_SS(offset)
            % Creator
            % CONSTELLATION REF -----------------------------------------------
            % CRS parameters, according to each GNSS system CRS definition
            % (ICD document in brackets):
            %
            % *_GLO --> PZ-90.02 (GLONASS-ICD 5.1)
            % Standard: http://gauss.gge.unb.ca/GLONASS.ICD-98.pdf
            %  newer -> http://kb.unavco.org/kb/assets/727/ikd51en.pdf
            % The new specifications are available only in Russian
            % L1 specification: http://russianspacesystems.ru/wp-content/uploads/2016/08/IKD-L1-s-kod.-razd.-Red-1.0-2016.pdf
            % L2 specification: http://russianspacesystems.ru/wp-content/uploads/2016/08/IKD-L2-s-kod.-razd.-Red-1.0-2016.pdf
            % L3 specification: http://russianspacesystems.ru/wp-content/uploads/2016/08/IKD-L3-s-kod.-razd.-Red-1.0-2016.pdf
            
            if (nargin == 0)
                offset = 0;
            end
            
            obj.char_id = 'R';
            obj.n_sat = 24;
            obj.prn = (1 : 24)';

            % http://www.navipedia.net/index.php/GLONASS_Signal_Plan            
            obj.f.base.R1  = 1602.000;                              % GLONASS (GLONASS-ICD 5.1) Freq [MHz]
            obj.f.base.R2  = 1246.000;                              % GLONASS (GLONASS-ICD 5.1) Freq [MHz]
            obj.f.base.R3  = 1201.000;                              % GLONASS (GLONASS-ICD 2016) Freq [MHz]
            obj.f.delta.R1 = 0.5625;                                % GLONASS (GLONASS-ICD 5.1) Freq [MHz]
            obj.f.delta.R2 = 0.4375;                                % GLONASS (GLONASS-ICD 5.1) Freq [MHz]
            obj.f.delta.R3 = 0.4375;                                % GLONASS (GLONASS-ICD 2016) Freq [MHz]
            obj.f.R_channels = 6:-1:-7;
            obj.f.R1 = obj.f.R_channels' .* obj.f.delta.R1 + obj.f.base.R1;
            obj.f.R2 = obj.f.R_channels' .* obj.f.delta.R1 + obj.f.base.R1;
            obj.f_vec = [obj.f.R1 obj.f.R2] * 1e6;                  % all the frequencies
            obj.l_vec = 299792458 ./ obj.f_vec;                     % lambda => wavelengths

            obj.init_iono(9, 7);

            obj.orbital_parameters.GM = 3.986004418e14;            % GLONASS (PZ-90.02) Gravitational constant * (mass of Earth) [m^3/s^2]
            obj.orbital_parameters.OMEGAE_DOT = 7.292115e-5;       % GLONASS (PZ-90.02) Angular velocity of the Earth rotation [rad/s]

            % Other useful links (the reference is the standard):
            % http://www.navipedia.net/index.php/Reference_Frames_in_GNSS
            % http://gage6.upc.es/eknot/Professional_Training/PDF/Reference_Systems.pdf
            obj.orbital_parameters.ell.a = 6378136.0;               % GLONASS (PZ-90.02) Ellipsoid semi-major axis [m]
            obj.orbital_parameters.ell.f = 1/298.25784;             % GLONASS (PZ-90.02) Ellipsoid flattening
            obj.orbital_parameters.ell.e = sqrt(1 - (1 - obj.orbital_parameters.ell.f) ^ 2); % GLONASS (PZ-90.02) Eccentricity
            obj.J2 = 1.0826257e-3;                                  % GLONASS second zonal harmonic of the geopotential            

            obj.update_indexes(offset);
        end
    end
end