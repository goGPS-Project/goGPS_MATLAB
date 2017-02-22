%   CLASS Galileo_SS
% =========================================================================
%
% DESCRIPTION
%   container of Galileo Satellite System parameters
%
% REFERENCES
%   CRS parameters, according to each GNSS system CRS definition
%   (ICD document in brackets):
%
%   *_GAL --> GTRF     (Galileo-ICD 1.1)
%   Standard: https://www.gsc-europa.eu/system/files/galileo_documents/Galileo_OS_SIS_ICD.pdf
%
%   Other useful links
%     - http://www.navipedia.net/index.php/Galileo_Signal_Plan
%     - http://www.navipedia.net/index.php/Reference_Frames_in_GNSS
%    Ellipsoid definition is actually coming from this presentation:
%     - http://gage6.upc.es/eknot/Professional_Training/PDF/Reference_Systems.pdf
%    at the moment of writing the Galileo Geodetic Reference Service Provider (GRSP) website (http://www.ggsp.eu/)
%    is actually offline, and the GTRF16v01 (or any other RF) cannot be found online.
%    Since GTRF seems to be based on ITRF (that does not define a reference ellispoid) http://itrf.ign.fr/faq.php?type=answer
%    The ellipsoid found in the presentation is used as reference
%   
%   Useful data for the future (PPP):
%     - https://www.gsc-europa.eu/support-to-developers/galileo-iov-satellite-metadata#3.2


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

classdef Galileo_SS < Satellite_System            
    properties (Constant, Access = 'public')        
        % System frequencies as struct [MHz]
        f = struct('E1', 1575.420, ...
                   'E5a', 1176.450, ...
                   'E5b', 1207.140, ...
                   'E5', 1191.795, ...
                   'E6', 1278.750) 
        
        % Array of supported frequencies [MHz]
        f_vec = struct2array(Galileo_SS.f) * 1e6;  
        
        % Array of the corresponding wavelength - lambda => wavelengths
        l_vec = 299792458 ./ Galileo_SS.f_vec;   
        
        char_id = 'E'     % Satellite system (ss) character id
        n_sat = 30;       % Maximum number of satellite in the constellation
        prn = (1 : 30)';  % Satellites id numbers as defined in the constellation
    end
    
    properties (Constant, Access = 'private')
        % GPS (WGS84) Ellipsoid semi-major axis [m]
        ell_a = 6378137;
        % GPS (WGS84) Ellipsoid flattening
        ell_f = 1/298.257222101;
        % GPS (WGS84) Ellipsoid Eccentricity^2
        ell_e2 = (1 - (1 - Galileo_SS.ell_f) ^ 2);
        % GPS (WGS84) Ellipsoid Eccentricity
        ell_e = sqrt(Galileo_SS.ell_e2);
    end
    
    properties (Constant, Access = 'public')
        % Structure of orbital parameters (ellipsoid, GM, OMEGA_EARTH_DOT)
        orbital_parameters = struct('GM', 3.986004418e14, ...               % Galileo (Galileo-ICD) Gravitational constant * (mass of Earth) [m^3/s^2]
                                    'OMEGAE_DOT', 7.2921151467e-5, ...      % Galileo (Galileo-ICD) Angular velocity of the Earth rotation [rad/s]
                                    'ell',struct( ...                       % Ellipsoidal parameters Galileo (GTRF)
                                        'a', Galileo_SS.ell_a, ...          % Ellipsoid semi-major axis [m]
                                        'f', Galileo_SS.ell_f, ...          % Ellipsoid flattening
                                        'e', Galileo_SS.ell_e, ...          % Eccentricity
                                        'e2', Galileo_SS.ell_e2));          % Eccentricity^2
    end
    
    methods
        function this = Galileo_SS(offset)
            % Creator            
            if (nargin == 0)
                offset = 0;
            end
            this.updateGoIds(offset);
        end
    end
end
