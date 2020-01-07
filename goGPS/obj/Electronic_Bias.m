
%   Electronic_Bias
% =========================================================================
%
% DESCRIPTION
%   Class to manage fixing
%
% EXAMPLE
%   EB = Electronic_Bias();
%
% SEE ALSO
% FOR A LIST OF CONSTANTs and METHODS use doc Network

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:                    Giulio Tagliaferro
%  Contributors:                  Elec
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
classdef Electronic_Bias < handle
    properties (Constant)
        CONST = 1;
        SPLINE_LIN = 2;
        SPLINE_CUBE = 3;
    end
    properties (GetAccess = public, SetAccess = public)
        o_code
        type
        data
        time_data %% WARNING must be constant rate
    end
    
    methods (Access = public)
        function this = Electronic_Bias(o_code, data, time_data, type)
            % constructor
            % SYNTAX
            % this = Electronic_Bias(data, <time_data>, <type>)
            this.o_code = o_code;
            this.data = data;
            if nargin == 2
                this.type = this.CONST;
            elseif nargin == 3
                this.type = this.SPLINE_LIN;
            elseif nargin == 4
                this.type = type;
                this.time_data = time_data;
            end
        end
        function bias = getBias(this,time)
            % get the electronic bias
            %
            % SYTANX
            % bias = getBias(this,<time>)
            if nargin == 1
                if this.type == this.CONST
                    bias = this.data;
                else
                    log  = Core.getLogger;
                    log.addError('Type not const but no time requested')
                end
            end
            if this.type == this.CONST
                bias = this.data;
            elseif this.type == this.SPLINE_LIN
                 spl_rate = this.time_data.getRate;
                 int_time = (time - this.time_data.first);
                 ep_id = floor(int_time/spl_rate);
                 spline_v = Core_Utils.spline(rem(int_time,spl_rate)/spl_rate,1);
                 bias = spline_v(:,1) .* this.data(ep_id +1) + spline_v(:,2) .* this.data(min(ep_id +2,length(this.data)));
            elseif this.type == this.SPLINE_CUBE
                 spl_rate = this.time_data.getRate;
                 int_time = (time - this.time_data.first);
                 ep_id = floor(int_time/spl_rate);
                 spline_v = Core_Utils.spline(rem(int_time,spl_rate)/spl_rate,3);
                 bias = spline_v(:,1) .* this.data(ep_id +1) + spline_v(:,2) .* this.data(min(ep_id +2,length(this.data))) + spline_v(:,3) .* this.data(min(ep_id +3,length(this.data))) + spline_v(:,4) .* this.data(min(ep_id +4,length(this.data)));
            end
            
        end
    end
end