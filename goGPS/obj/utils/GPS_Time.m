%   CLASS GPS_Time
% =========================================================================
%
% DESCRIPTION
%   Class to manage times and dates in various format (GPS / UTC/ ...)
%
% EXAMPLE
%   settings = GPS_Time();
%
% FOR A LIST OF CONSTANTs and METHODS use doc goGNSS

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
classdef GPS_Time < handle
    
    properties (Constant)
        DEFAULT_DATE_FORMAT = 'YYYY/MM/DD hh:mm:ss.ssss'; % String representing the format of visualization of the time
    end
    
    properties (SetAccess = private, GetAccess = public)
    end
    
    methods
        function obj = Receiver()
        end
    end
    
    methods 
    end    
    
    methods (Static, Access = 'protected')
        function setTimeFormat(time_format)
            persistent time_format__;
            time_format__ = time_format;
        end
        
        function time_format = getTimeFormat()
            persistent time_format__;
            if isempty(time_format__)
                time_format__ = DEFAULT_DATE_FORMAT;
            end
            time_format = time_format__;
        end
    end
end
