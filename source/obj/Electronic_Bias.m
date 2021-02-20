
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
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 (GReD srl) Giulio Tagliaferro
%  Written by:        Giulio Tagliaferro
%  Contributors:      Giulio Tagliaferro, ...
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
        SPLINE_CUB = 6;
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
                this.type = LS_Parametrization.CONST;
            elseif nargin == 3
                this.type = LS_Parametrization.SPLINE_LIN;
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
            log  = Core.getLogger;

            if nargin == 1
                if this.type == this.CONST
                    bias = this.data;
                else
                    log.addError('Type not const but no time requested')
                end
            end
            if this.type == LS_Parametrization.CONST
                bias = this.data;
            elseif this.type == LS_Parametrization.SPLINE_LIN
                 spl_rate = this.time_data.getRate;
                 int_time = (time - this.time_data.first);
                 ep_id = floor(int_time/spl_rate);
                 spline_v = Core_Utils.spline(rem(int_time,spl_rate)/spl_rate,1);
                 bias = spline_v(:,1) .* this.data(ep_id +1) + spline_v(:,2) .* this.data(min(ep_id +2,length(this.data)));
            elseif this.type == LS_Parametrization.SPLINE_CUB
                 spl_rate = this.time_data.getRate;
                 int_time = (time - this.time_data.first);
                 ep_id = floor(int_time/spl_rate);
                 interp_ep = ep_id >= 0 & ep_id < length(this.data);
                 bias = zeros(this.time_data.length,1);
                 if any(interp_ep)
                     spline_v = Core_Utils.spline(rem(int_time(interp_ep),spl_rate)/spl_rate,3);
                     bias(interp_ep) = spline_v(:,1) .* this.data(ep_id(interp_ep) +1)' + spline_v(:,2) .* this.data(min(ep_id(interp_ep) +2,length(this.data)))' + spline_v(:,3) .* this.data(min(ep_id(interp_ep) +3,length(this.data)))' + spline_v(:,4) .* this.data(min(ep_id(interp_ep) +4,length(this.data)))';
                 end
                 if any(~interp_ep)
                     bf_ep = ep_id < 0;
                     bias(bf_ep) = bias(find(interp_ep,1,'first'));
                     aft_ep = ep_id >= length(this.data);
                     bias(aft_ep) =bias(find(interp_ep,1,'last'));
                     %log.addWarning('Epoch requested out of time boundaries, extrapolating');
                 end
            end
            
        end
    end
end