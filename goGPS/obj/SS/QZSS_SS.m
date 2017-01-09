%   CLASS QZSS_SS
% =========================================================================
%
% DESCRIPTION
%   container of QZSS Satellite System parameters
%
%----------------------------------------------------------------------------------------------
%                           goGPS v0.9.1
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
classdef QZSS_SS < Satellite_System
    properties
        go_ids       % Satellites unique id numbers in goGPS
    end
            
    methods
        function obj = QZSS_SS(offset)
            % Creator
            % CONSTELLATION REF -----------------------------------------------
            % CRS parameters, according to each GNSS system CRS definition
            % (ICD document in brackets):
            %
            % *_QZS --> GRS80    (IS-QZSS 1.8E)
            % Standard: http://qz-vision.jaxa.jp/USE/is-qzss/DOCS/IS-QZSS_18_E.pdf
                        
            if (nargin == 0)
                offset = 0;
            end
            
            obj.char_id = 'J';
            obj.n_sat = 4;
            obj.prn = (193 : 196)';

            % http://www.navipedia.net/index.php/QZSS_Signal_Plan
            obj.f.J1  = 1575.420;                                   % QZSS Freq [MHz]
            obj.f.J2  = 1227.600;                                   % QZSS Freq [MHz]
            obj.f.J5  = 1176.450;                                   % QZSS Freq [MHz]
            obj.f.J6  = 1278.750;                                   % QZSS Freq [MHz]
            obj.f_vec = struct2array(obj.f) * 1e6;                  % all the frequencies
            obj.l_vec = 299792458 ./ obj.f_vec;                     % lambda => wavelengths

            obj.initIono(77, 60);
            
            obj.orbital_parameters.GM = 3.986005e14;                % QZSS (IS-QZSS 1.8E) Gravitational constant * (mass of Earth) [m^3/s^2]
            obj.orbital_parameters.OMEGAE_DOT = 7.2921151467e-5;    % QZSS (IS-QZSS 1.8E) Angular velocity of the Earth rotation [rad/s]

            % Other useful links (the reference is the standard):
            % http://www.navipedia.net/index.php/Reference_Frames_in_GNSS
            % http://gage6.upc.es/eknot/Professional_Training/PDF/Reference_Systems.pdf
            obj.orbital_parameters.ell.a = 6378137;                 % QZSS (GRS80) Ellipsoid semi-major axis [m]
            obj.orbital_parameters.ell.f = 1/298.257222101;         % QZSS (GRS80) Ellipsoid flattening
            obj.orbital_parameters.ell.e = sqrt(1 - (1 - obj.orbital_parameters.ell.f) ^ 2);      % QZSS (GRS80) Eccentricity            

            obj.updateGoIds(offset);
        end
    end
end
