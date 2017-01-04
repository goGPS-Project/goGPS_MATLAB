%   CLASS BeiDou
% =========================================================================
%
% DESCRIPTION
%   container of BeiDou Satellite System parameters
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
classdef BeiDou_SS < Satellite_System
    properties
        go_ids       % Satellites unique id numbers in goGPS
    end
            
    methods
        function obj = BeiDou_SS(offset)
            % Creator
            % CONSTELLATION REF -----------------------------------------------
            % CRS parameters, according to each GNSS system CRS definition
            % (ICD document in brackets):
            %
            % *_BDS --> CGCS2000 (BeiDou-ICD 1.0)
            % Standard: http://www.beidou.gov.cn/attach/2012/12/27/201212273da29c5eb8274deb8cd2b178228ba2bd.pdf
            
            if (nargin == 0)
                offset = 0;
            end
            
            obj.char_id = 'C';
            obj.n_sat = 37;
            obj.prn = (1 : 37)';

            % http://www.navipedia.net/index.php/BeiDou_Signal_Plan            
            obj.f.C2  = 1561.098;                                   % BeiDou (BeiDou-ICD 1.0) Freq [MHz]
            obj.f.C5b = 1207.140;                                   % BeiDou (BeiDou-ICD 1.0) Freq [MHz]
            obj.f.C6  = 1268.520;                                   % BeiDou (BeiDou-ICD 1.0) Freq [MHz]
            obj.f.C1  = 1589.740;                                   % BeiDou (BeiDou-ICD 1.0) Freq [MHz]
            obj.f_vec = struct2array(obj.f) * 1e6;                  % all the frequencies
            obj.l_vec = 299792458 ./ obj.f_vec;                     % lambda => wavelengths

            obj.initIono(763, 590);

            obj.orbital_parameters.GM = 3.986004418e14;             % BeiDou (BeiDou-ICD 1.0) Gravitational constant * (mass of Earth) [m^3/s^2]
            obj.orbital_parameters.OMEGAE_DOT = 7.2921150e-5;       % BeiDou (BeiDou-ICD 1.0) Angular velocity of the Earth rotation [rad/s]

            % Other useful links (the reference is the standard):
            % http://www.navipedia.net/index.php/Reference_Frames_in_GNSS
            % http://gage6.upc.es/eknot/Professional_Training/PDF/Reference_Systems.pdf
            obj.orbital_parameters.ell.a = 6378137.0;               % BeiDou (CGCS2000) Ellipsoid semi-major axis [m]
            obj.orbital_parameters.ell.f = 1/298.257222101;         % BeiDou (CGCS2000) Ellipsoid flattening
            obj.orbital_parameters.ell.e = sqrt(1 - (1 - obj.orbital_parameters.ell.f) ^ 2); % BeiDou (CGCS2000)    Eccentricity

            obj.updateGoIds(offset);
        end
    end
end
