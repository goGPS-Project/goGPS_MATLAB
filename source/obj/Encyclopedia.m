%   SINGLETON CLASS Encyclopedia
% =========================================================================
%
% DESCRIPTION
%   This class allows the declaration of the definitions
%
% EXAMPLE
%   log = Encyclopedia.getInstance();
%
%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0
%
%--------------------------------------------------------------------------
%  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
%  Written by:        Giulio Tagliaferro
%  Contributors:      Giulio Tagliaferro
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

classdef Encyclopedia < handle
    properties
        observation_group = {};
    end
    methods (Static)
        function this = getInstance(std_out)
            % Concrete implementation.  See Singleton superclass.
            persistent unique_instance_encyclopedia__
            if isempty(unique_instance_encyclopedia__)
                this = Logger();
                unique_instance_encyclopedia__ = this;
            else
                this = unique_instance_encyclopedia__;
            end
        end
    end
    
    methods 
        function addObsGroup(this, obs_group)
            this.observation_group{end+1} = obs_group;
        end
        
        function id = getObsGroupId(this, obs_group)
                id = Core_Utils.findAinB(obs_group, this.observation_group);
                if id == 0
                    this.addObsGroup(obs_group);
                    id = length(this.observation_group);
                end
        end
    end
end
