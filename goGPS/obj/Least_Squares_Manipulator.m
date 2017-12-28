%   CLASS Least_Square_Manipulator
% =========================================================================
%
% DESCRIPTION
%   Efficently manipulate sparse least squares system 
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
%    |___/                    v 0.6.0 alpha 1 - nightly
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
classdef Least_Squares_Manipulator < handle
    properties
        A % Values of the normal matrices [n_obs x n_param_per_epoch]
        A_idx % index of the paramter [n_obs x n_param_per_epoch]
        N %
        y % observations  [ n_obs x 1]
        epoch % epoch of the obseravtions and of the A lines [ n_obs x 1]
        param_class % [n_param x 1] each paramter can be part of a class [ 1 : x , 2 : y , 3 : z, 4: clock, 5: inter channel/frequency/system biases,
        %                                                                  6 : tropo, 7: tropo inclination north, 8: tropo inclination east ]
    end
    properties (Access = private)
        Ncc % part of the normal matrix with costant paramters
        Nee % diagonal part of the normal matrix with epoch wise or multi epoch wise paramters
        Nce % cross element between constant and epoch varying paramters 
    end
    methods
        function this = Least_Square_Manipulator()
        end
        function [x, res, s02, Cxx] = solve(this)
        end
        function regularize(this, reg_opt)
        end
        function setUpPPP(this, rec, ppp_opt)
            % get double frequency iono_free for all the systems
            obs_set = Observation_Set();
            for s = rec.cc.sys_c
                obs_set.merge(rec.getPrefIonoFree('L',s));
            end
            % set up number of parametrs requires
            synt_obs = rec.getSyntTwin(obs_set);
            diff_obs = obs_set.obs - synt_obs;
            n_epochs = rec.time.length;
            n_coo = 3;
            n_clocks = n_epochs;
            n_iob = []; 
        end
    end
end