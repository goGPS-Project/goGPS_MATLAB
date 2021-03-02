%   CLASS Reflectometry
% =========================================================================
%
% DESCRIPTION
%   Class to manage reflectometry functions
%
% EXAMPLE
%   rfx = Reflectometry();
%
% CONSTRUCTOR SYNTAX
%   rfx = Reflectometry();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Reflectometry
%
% COMMENTS
% The class stores arrays of time, not just a single element,
% it has been designed this way because MATLAB works best on arrays
%

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 Geomatics Research & Development srl (GReD)
%  Written by:        Andrea Gatti
%  Contributors:      Andrea Gatti, ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------


classdef Reflectometry < handle
    
    properties (Constant, GetAccess = public)
    end
    
    properties (SetAccess = public, GetAccess = public) % set permission have been changed from private to public (Giulio)
    end
        
    % =========================================================================
    %    
    % =========================================================================
    
    methods
        function this = Reflectometry()
            % Constructor
            %
            % SYNTAX
            %   t = Coordinates();            
        end         
    end
    
    % =========================================================================
    %    TESTS
    % =========================================================================
    
    methods (Static, Access = 'public')
        function getAvgNyquist(rec, el_lim)
            % Show fresnel zones for a receiver
            %
            % INPUT
            %   rec     receiver
            %   el_lim  elevation limits [deg] (e.g. [5 15])
            %
            % SYNTAX
            %   Reflectometry.getAvgNyquist(rec, getAvgNyquist)
            log = Core.getLogger();
            cc = Core.getConstellationCollector;
            
            diff_el_crit = abs(diff(el_lim)) - 2;
            for r = 1 : numel(rec)
                work = rec(r).work;
                log.addMessage('Displaying fresnel zones for the receiver rec %s', rec(r).getMarkerName);
                [snr, id_snr] =  = work.getSNR();
                id_snr = find(id_snr);
                [az, el] = work.getAzEl(work.go_id(id_snr));
                azi = 0;
                [lat, lon, h_ellips] = work.getMedianPosGeodetic();
                freq_list = freq_list = unique(work.obs_code(id_snr,2));                
                for f = 1 : numel(freq_list)
                    id_f = find(work.obs_code(id_snr,2) == freq_list(f));
                    for s = 1 : numel(id_f)
                        sat_name = cc.getSatName(work.go_id(id_snr(id_f(s))));
                        % find data within 20 degrees of azimuth and using elevation angle limits
                        id_ok = el(:,id_f(s)) >= el_lim(1) & el(:,id_f(s)) >= el_lim(2);
                        % sin of elevation
                        sin_el = sind(el(id_ok, id_f(s)));
                        % elevation limits
                        diff_el = max(el(id_ok, id_f(s))) - min(el(id_ok, id_f(s)));
                        diff_t = 2 * (max(sin_el) - min(sin_el)) / work.wl(id_snr(id_f(s)));
                        if diff_el > diff_el_crit
                            fprintf(1,'%s Azim. %6.2f  Nyquist %7.2f (m) %s ElevAngles: %3.0f %3.0f \n', ...
              sat_name, azi, avg_nyq, riseORset(riseset), minElevAngle, maxElevAngle);
                            avg_nyq = sum(id_ok) / (2 * diff_t);  
                        end
                    end
                end
            end
        end
    end
    
end
