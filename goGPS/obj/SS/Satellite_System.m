%   CLASS Satellite_System
% =========================================================================
%
% DESCRIPTION
%   Abstract class with basic properties and methods that a satellite
%   system must have
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

classdef Satellite_System < handle
    properties (Constant, Abstract)
        f             % System frequencies as struct
        f_vec         % Array of supported frequencies
        l_vec         % Array of the corresponding wavelength (pre-computed for performance reason)
        
        char_id       % Satellite system (ss) character id
        n_sat         % Maximum number of satellite in the constellation
        prn           % Satellites id numbers as defined in the constellation
        
        % Structure of orbital parameters (ellipsoid, GM, OMEGA_EARTH_DOT)
        orbital_parameters 
        %   .GM
        %   .OMEGAE_DOT
        %   .J2 (only for GLONASS)        
        %   .ell.a      semi-major axis
        %   .ell.f      flattening
        %   .ell.e      eccentricity
        %   .ell.e2     eccentricity^2
    end
    
    properties
        go_ids       % Satellites unique id numbers in goGPS
    end
        
    methods (Access = public)        
        function iono_free = getIonoFree(this, f_id)
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
            iono_free.alpha1 = this.f_vec(f_id(1)) .^ 2 ./ (this.f_vec(f_id(1)) .^ 2 - this.f_vec(f_id(2)) .^ 2);
            iono_free.alpha2 = this.f_vec(f_id(2)) .^ 2 ./ (this.f_vec(f_id(1)) .^ 2 - this.f_vec(f_id(2)) .^ 2);
            gcd_f = gcd(this.f_vec(f_id(1)),this.f_vec(f_id(2)));
            iono_free.T = this.f_vec(f_id(1))/gcd_f;
            iono_free.N = this.f_vec(f_id(2))/gcd_f;
        end        
    end
    
    methods
        function updateGoIds(this, offset)
            %  Update the satellites unique id numbers in goGPS
            if nargin == 1
                offset = 0;
            end
            this.go_ids = offset + (1 : this.n_sat);
        end        
    end
end
