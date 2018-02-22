%   CLASS Core_SEID
% =========================================================================
%
% DESCRIPTION
%   class to compute a SEID processing
%
% EXAMPLE
%   seid = Core_SEID();
%
% FOR A LIST OF CONSTANTs and METHODS use doc goGNSS


%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.6.0 alpha 1 - nightly
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

classdef Core_SEID < handle

    properties % Public Access
        log;
        state;
    end

    methods (Static)
        function this = Core_SEID()
            % Core object creator
            this.log = Logger.getInstance();
            this.state = Global_Configuration.getCurrentSettings();
        end

        % Concrete implementation.  See Singleton superclass.
        function this = getInstance()
            % Get the persistent instance of the class
            persistent unique_instance_core_seid_
            unique_instance_core_seid_ = [];

            if isempty(unique_instance_core_seid_)
                this = Core();
                unique_instance_core_seid_ = this;
            else
                this = unique_instance_core_seid_;
            end
        end
    end

    % =========================================================================
    %  CONSTELLATION MANAGEMENT
    % =========================================================================

    methods % Public Access
    end

    methods % Public Access (Legacy support)
    end

end
