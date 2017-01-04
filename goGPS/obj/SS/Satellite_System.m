%   CLASS Satellite_System
% =========================================================================
%
% DESCRIPTION
%   Abstract class with basic properties and methods that a satellite
%   system must have
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
classdef Satellite_System < handle
    properties (SetAccess = protected, GetAccess = public)
        f             % System frequencies
        f_vec         % Array of supported frequencies
        l_vec         % Array of the corresponding wavelength (pre-computed for performance reason)
        
        char_id       % Satellite system (ss) character id
        n_sat         % Maximum number of satellite in the constellation
        prn           % Satellites id numbers as defined in the constellation
        
        % Structure of orbital parameters (ellipsoid, GM, OMEGA_EARTH_DOT)
        orbital_parameters 
        %   .ell.a      semi-major axis
        %   .ell.f      flattening
        %   .ell.e      eccentricity
        %   .GM
        %   .OMEGAE_DOT
        %   .J2 (only for GLONASS)
        
        % Structure of the iono free combination parameters (alpha 1, alpha 2, T ed N
        iono_free
        %   .alpha1
        %   .alpha2
        %   .T
        %   .N
    end
    
    properties (Abstract)
        go_ids       % Satellites unique id numbers in goGPS
    end
    
    methods (Access = protected)
        function initIono(obj, T, N)
            % Init iono parameters: iono-free combination is computed with
            % the first two carriers in f_vec
            obj.iono_free.alpha1 = obj.f_vec(1) .^ 2 ./ (obj.f_vec(1) .^ 2 - obj.f_vec(2) .^ 2);
            obj.iono_free.alpha2 = obj.f_vec(2) .^ 2 ./ (obj.f_vec(1) .^ 2 - obj.f_vec(2) .^ 2);
            obj.iono_free.T = T;
            obj.iono_free.N = N;
        end
        
    end
    
    methods
        function updateGoIds(obj, offset)
            %  Update the satellites unique id numbers in goGPS
            if nargin == 1
                offset = 0;
            end
            obj.go_ids = offset + (1 : obj.n_sat);
        end        
    end
end
