%   CLASS Galileo_SS
% =========================================================================
%
% DESCRIPTION
%   container of Galileo Satellite System parameters
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
classdef Galileo_SS < Satellite_System
    properties
        indexes       % Satellites unique id numbers in goGPS
    end
            
    methods
        function obj = Galileo_SS(offset)
            % Creator
            % CONSTELLATION REF -----------------------------------------------
            % CRS parameters, according to each GNSS system CRS definition
            % (ICD document in brackets):
            %
            % *_GAL --> GTRF     (Galileo-ICD 1.1)
            % Standard: https://www.gsc-europa.eu/system/files/galileo_documents/Galileo_OS_SIS_ICD.pdf
            
            if (nargin == 0)
                offset = 0;
            end
            
            obj.char_id = 'E';
            obj.n_sat = 30;
            obj.prn = (1 : 30)';

            % http://www.navipedia.net/index.php/Galileo_Signal_Plan
            obj.f.E1  = 1575.420;                                   % Galileo Freq [MHz]
            obj.f.E5a = 1176.450;                                   % Galileo Freq [MHz]
            obj.f.E5b = 1207.140;                                   % Galileo Freq [MHz]
            obj.f.E5  = 1191.795;                                   % Galileo Freq [MHz]
            obj.f.E6  = 1278.750;                                   % Galileo Freq [MHz]
            obj.f_vec = struct2array(obj.f) * 1e6;                  % all the frequencies
            obj.l_vec = 299792458 ./ obj.f_vec;                     % lambda => wavelengths

            obj.init_iono(154, 115);
            
            obj.orbital_parameters.GM = 3.986004418e14;             % Galileo (Galileo-ICD) Gravitational constant * (mass of Earth) [m^3/s^2]
            obj.orbital_parameters.OMEGAE_DOT = 7.2921151467e-5;    % Galileo (Galileo-ICD) Angular velocity of the Earth rotation [rad/s]

            % Other useful links (the reference is the standard):
            % http://www.navipedia.net/index.php/Reference_Frames_in_GNSS
            % Ellipsoid definition is actually coming from this presentation:
            % http://gage6.upc.es/eknot/Professional_Training/PDF/Reference_Systems.pdf
            % at the moment of writing the Galileo Geodetic Reference Service Provider (GRSP) website (http://www.ggsp.eu/)
            % is actually offline, and the GTRF16v01 (or any other RF) cannot be found online.
            % Since GTRF seems to be based on ITRF (that does not define a reference ellispoid) http://itrf.ign.fr/faq.php?type=answer
            % The ellipsoid found in the presentation is used as reference
            obj.orbital_parameters.ell.a = 6378137;                 % Galileo (GTRF) Ellipsoid semi-major axis [m]
            obj.orbital_parameters.ell.f = 1/298.257222101;         % Galileo (GTRF) Ellipsoid flattening
            obj.orbital_parameters.ell.e = sqrt(1 - (1 - obj.orbital_parameters.ell.f) ^ 2);      % Galileo (GTRF)    Eccentricity
            
            % Useful data for the future (PPP):
            % https://www.gsc-europa.eu/support-to-developers/galileo-iov-satellite-metadata#3.2

            obj.update_indexes(offset);
        end
    end
end