%   CLASS GPS_SS
% =========================================================================
%
% DESCRIPTION
%   container of GPS Satellite System parameters
%

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

classdef GPS_SS < Satellite_System
    properties
        go_ids       % Satellites unique id numbers in goGPS
    end
            
    methods
        function obj = GPS_SS(offset)            
            % Creator
            % CONSTELLATION REF -----------------------------------------------
            % CRS parameters, according to each GNSS system CRS definition
            % (ICD document in brackets):
            %
            % *_GPS --> WGS-84   (IS-GPS200E)
            % Standard IS-GPS-200H: http://www.gps.gov/technical/icwg/IS-GPS-200H.pdf
            
            if (nargin == 0)
                offset = 0;
            end
            
            obj.char_id = 'G';
            obj.n_sat = 32;
            obj.prn = (1 : 32)';

            % http://www.navipedia.net/index.php/GPS_Signal_Plan            
            obj.f.L1 = 1575.420;                                    % GPS (IS-GPS200H) Freq [MHz]
            obj.f.L2 = 1227.600;                                    % GPS (IS-GPS200H) Freq [MHz]
            obj.f.L5 = 1176.450;                                    % GPS (IS-GPS200H) Freq [MHz]
            obj.f_vec = struct2array(obj.f) * 1e6;                  % all the frequencies
            obj.l_vec = 299792458 ./ obj.f_vec;                     % lambda => wavelengths

            obj.initIono(77, 60);

            obj.orbital_parameters.GM = 3.986005e14;                % GPS (IS-GPS200H) Gravitational constant * (mass of Earth) [m^3/s^2]
            obj.orbital_parameters.OMEGAE_DOT = 7.2921151467e-5;    % GPS (IS-GPS200H) Angular velocity of the Earth rotation [rad/s]

            % Other useful links
            % Ellipsoid: http://www.unoosa.org/pdf/icg/2012/template/WGS_84.pdf
            % note that GM and OMEGAE_DOT are redefined in the standard IS-GPS200H (GPS is not using WGS_84 values)
            % http://www.navipedia.net/index.php/Reference_Frames_in_GNSS
            % http://gage6.upc.es/eknot/Professional_Training/PDF/Reference_Systems.pdf
            obj.orbital_parameters.ell.a = 6378137;                 % GPS (WGS84) Ellipsoid semi-major axis [m]
            obj.orbital_parameters.ell.f = 1/298.257223563;         % GPS (WGS84) Ellipsoid flattening
            obj.orbital_parameters.ell.e = sqrt(1 - (1 - obj.orbital_parameters.ell.f) ^ 2); % GPS (WGS84)    Eccentricity

            obj.updateGoIds(offset);
        end
    end
end
