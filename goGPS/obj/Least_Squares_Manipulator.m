%   CLASS Least_Square_Manipulator
% =========================================================================
%
% DESCRIPTION
%   Manipulate a least square system
%
% EXAMPLE
%   LSM = Least_Square_Manipulator();
%
% SEE ALSO
%   - Least_Square
% FOR A LIST OF CONSTANTs and METHODS use doc Main_Settings

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.2 beta 1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Giulio Tagliaferro
%  Contributors:     
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
classdef Least_Squares_Manipulator < Least_Squares
    properties
        y0_epoch; % [n_obs x 1] epoch of the observation
    end
    methods
        function this = Least_Square_Manipulator()
        end
        function sortSystemByEpoch(this)
            [ep_sort, idx] = sort(this.y0_epoch);
            this.A  = this.A(idx, :);
            this.y0 = this.y0(idx);
            this.b  = this.b(idx);
            this.Q = this.Q(idx, :);
            this.Q = this.Q(:, idx);
            if not(isempty(this.res))
                this.res = this.res(idx);
            end
            this.y0_epoch = ep_sort;
        end
        function splitParamEpochwise(this, col, rate)
            
        end
        
    end
    
end