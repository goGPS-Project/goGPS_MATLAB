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

classdef QZSS_SS < Satellite_System
    properties (Constant, Access = 'public')
        % CONSTELLATION REF -----------------------------------------------
        
        % System frequencies as struct [MHz]
        f = struct('J1', 1575.420, ... 
                   'J2', 1227.600, ...
                   'J5', 1176.450, ...
                   'J6', 1278.750)
        
        % Array of supported frequencies [MHz]
        f_vec = struct2array(QZSS_SS.f) * 1e6;  
        
        % Array of the corresponding wavelength - lambda => wavelengths
        l_vec = 299792458 ./ GPS_SS.f_vec;   
        
        char_id = 'J'        % Satellite system (ss) character id
        n_sat = 4;           % Maximum number of satellite in the constellation
        prn = (193 : 196)';  % Satellites id numbers as defined in the constellation
    end
    
    properties (Constant, Access = 'private')
        % QZSS (GRS80) Ellipsoid semi-major axis [m]
        ell_a = 6378137;
        % QZSS (GRS80) Ellipsoid flattening
        ell_f = 1/298.257222101;
        % QZSS (GRS80) Ellipsoid Eccentricity^2
        ell_e2 = (1 - (1 - QZSS_SS.ell_f) ^ 2);
        % QZSS (GRS80) Ellipsoid Eccentricity
        ell_e = sqrt(QZSS_SS.ell_e2);
    end
    
    properties (Constant, Access = 'public')
        % Structure of orbital parameters (ellipsoid, GM, OMEGA_EARTH_DOT)
        orbital_parameters = struct('GM', 3.986005e14, ...                  % Gravitational constant * (mass of Earth) [m^3/s^2]
                                    'OMEGAE_DOT', 7.2921151467e-5, ...      % Angular velocity of the Earth rotation [rad/s]
                                    'ell',struct( ...                       % Ellipsoidal parameters QZSS (GRS80)
                                        'a', QZSS_SS.ell_a, ...             % Ellipsoid semi-major axis [m]
                                        'f', QZSS_SS.ell_f, ...             % Ellipsoid flattening
                                        'e', QZSS_SS.ell_e, ...             % Eccentricity
                                        'e2', QZSS_SS.ell_e2));             % Eccentricity^2
    end
    
    methods
        function this = QZSS_SS(offset)
            % Creator                        
            if (nargin == 0)
                offset = 0;
            end
            this.updateGoIds(offset);
        end
    end
end
