%   CLASS Core_Antenna
% =========================================================================
%
% DESCRIPTION
%   Collector of antenna models
%   Container of the ATX file()
%
% EXAMPLE
%   atx = Core_Antenna();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Core_Antenna


%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 2
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti, ...
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

classdef Core_Antenna < handle
    
    %% PROPERTIES CONSTANTS
    % ==================================================================================================================================================
    properties (Constant)
    end
    
    %% PROPERTIES MISC
    % ==================================================================================================================================================
    properties (GetAccess = private, SetAccess = private) % Public Access
        creation_time = GPS_Time(now);
    end
    
    %% PROPERTIES ANTENNA
    % ==================================================================================================================================================
    properties
        list        % Array of Antenna objects
        type        % List of antenna names present in ant_list (mirrored copy) (for look-up)
        serial      % List of Antennas serial number (for look-up)
        start  % Validity interval (MATLAB TIME)
        stop   % Validity interval (MATLAB TIME)
    end
    
    %% METHOD CREATOR
    % ==================================================================================================================================================
    methods (Static, Access = public)
        % Concrete implementation.  See Singleton superclass.
        function this = Core_Antenna()
            % Core object creator
        end
    end
    
    %% METHODS INIT & STATIC GETTERS & SETTERS
    % ==================================================================================================================================================
    methods (Static, Access = public)
        function this = fromAntex(file_name)
            this = Core_Antenna();
            this.init();
            this.importAntex(file_name);
        end
    end
    
    %% METHODS INIT
    % ==================================================================================================================================================
    methods
        function init(this, atx_file_name)
            this.type = '';
            this.serial = '';
            this.list = [];
            if nargin == 2
                this.importAntex(atx_file_name);
            end
        end
        
        function importAntex(this, file_name)
            if nargin == 0
                file_name = Core.getState.getAtxFile();
            end
            log = Core.getLogger();
            
            % open RINEX observation file
            fid = fopen(file_name,'r');
            txt = fread(fid,'*char')';
            % try to see if carriage return is present in the file (Windows stupid standard)
            % On Windows file lines ends with char(13) char(10)
            % instead of just using char(10)
            if ~isempty(find(txt(1:min(1000,numel(txt))) == 13, 1, 'first'))
                has_cr = true;  % The file has carriage return - I hate you Bill!
            else
                has_cr = false;  % The file is UNIX standard
            end
            % txt = txt(txt ~= 13);  % remove carriage return - I hate you Bill!
            fclose(fid);
            
            % get new line separators
            nl = regexp(txt, '\n')';
            if nl(end) <  (numel(txt) - double(has_cr))
                nl = [nl; numel(txt)];
            end
            lim = [[1; nl(1 : end - 1) + 1] (nl - 1 - double(has_cr))];
            lim = [lim (lim(:,2) - lim(:,1) + 1)];
            lim(lim(:,3) < 64, :) = []; % Do not consider any line shorter than 65 characters
            lim((txt((lim(:,1) + 60)) == 'C'), :) = []; % detect COMMENT and ignore them
            
            id_start = find((txt((lim(:,1) + 60)) == 'S') & (txt((lim(:,1) + 69)) == 'A'));  % detect START OF ANTENNA
            id_stop = find((txt((lim(:,1) + 60)) == 'E') & (txt((lim(:,1) + 69)) == 'T'));   % detect END OF ANTENNA
            id_ant = find((txt((lim(:,1) + 60)) == 'T'));
            n_ant = numel(id_start);
            if n_ant ~= numel(id_stop) || n_ant ~= numel(id_ant)
                log.addError('Corrupted antex file!\n antennas are not properly opened or closed by (START / END)');
            else
                ant_type = txt(repmat(lim(id_ant, 1), 1, 20) + repmat(0:19, n_ant, 1));
                ant_serial = txt(repmat(lim(id_ant, 1), 1, 20) + repmat(20:39, n_ant, 1)); % Used for Satellites antenna
                ant_list = Antenna();
                for a = 1 : n_ant
                    ant_list(a) = Antenna.fromText(txt(lim(id_start(a) + 1, 1) : lim(id_stop(a) - 1, 2)), lim((id_start(a) + 1) : (id_stop(a) - 1), :), has_cr);
                end
            end
            
            % Append imported antennas
            if isempty(this.type)
                this.type = ant_type;
                this.serial = ant_serial;
                this.list = ant_list;
                tmp = [ant_list.start];
                this.start = tmp.getMatlabTime;
                tmp = [ant_list.stop];
                this.stop = tmp.getMatlabTime;
            else
                this.type = [this.type; ant_type];
                this.serial = [this.serial; ant_serial];
                this.list = [this.list; ant_list];
                tmp = [ant_list.start];
                this.start = [this.start; tmp.getMatlabTime];
                tmp = [ant_list.stop];
                this.stop = [this.stop tmp.getMatlabTime];
            end
        end
    end
    
    %% METHODS UTILITIES
    % ==================================================================================================================================================
    methods
        function sat_type = getSatAntennaType(this, serial, time)
            % get the satellite type given a 3ch antenna (e.g. G01, E12, ...)
            %
            % INPUT
            %   serial      3ch serial value
            %   time        GPS_Time reference time for antenna validity
            %
            % SYNTAX
            %   sat_type = this.getSatAntennaType(serial, time)
            if ~iscell(serial)
                serial = {serial};
            end
            time = time.getMatlabTime;
            ss = Core_Utils.code3Char2Num(this.serial(:, 1:3)); % stored serials expressed as a number
            sat_type = cell(size(serial(:), 1), 1);
            for a = 1 : size(serial(:), 1)
                id_sat = find(ss == Core_Utils.code3Char2Num(serial{a}(1:3)) & ...
                    this.start <= time & ...
                    this.stop >= time, 1, 'first');
                sat_type{a} = this.type(id_sat, :);
            end
        end
        
        function ant = getAntenna(this, type, serial, time, log_lev)
            % get the satellite type given a 3ch antenna (e.g. G01, E12, ...)
            %
            % INPUT
            %   serial      3ch serial value
            %   time        GPS_Time reference time for antenna validity
            %
            % SYNTAX
            %   sat_type = this.getAntenna(type, serial, time)
            if ~isempty(type) && ~iscell(type)
                type = {type};
            end
            if ~isempty(serial) && ~iscell(serial)
                serial = {serial};
            end
            
            log = Core.getLogger;
            
            time = time.getMatlabTime;
            % time check
            t_check = find(this.start <= time & this.stop >= time);
            
            n_ant = max(size(serial(:), 1), size(type(:), 1));
            none_code = Core_Utils.code4Char2Num('NONE');
            ant(n_ant + 1) = Antenna(); ant(end) = []; % init output, define empty array of Antennas
            for a = 1 : n_ant
                if isempty(type) || isempty(type{a})
                    cmp_str = this.serial;
                    ant_name = serial{a};
                elseif isempty(serial) || isempty(serial{a})
                    cmp_str = this.type;
                    ant_name = type{a};
                else
                    cmp_str = [this.type this.serial];
                    ant_name = [type{a} serial{a}];
                end
                
                if (numel(ant_name) >= 4) && Core_Utils.code4Char2Num(ant_name(1:4)) == none_code
                    % ignoring antenna type == 'NONE'
                else
                    % Get the antenna id containing the searched antenna file
                    id_ant = ((strfind(serialize(cmp_str')', ant_name) - 1) / size(cmp_str, 2) + 1);
                    id_ant(rem(id_ant,1) > 0) = []; % remove wrong matches
                    if isempty(id_ant) && ~(isempty(type) || isempty(type{a})) && (Core_Utils.code4Char2Num(ant_name(17:20)) ~= none_code) % search for same antenna with no radome
                        % Check with no radome
                        if isempty(serial) || isempty(serial{a})
                            cmp_str = this.type;
                            ant_name = type{a};
                            ant_name(17:20) = 'NONE';
                        else
                            cmp_str = [this.type this.serial];
                            ant_name = [type{a} serial{a}];
                            ant_name(17:20) = 'NONE';
                        end
                        id_ant = ((strfind(serialize(cmp_str')', ant_name) - 1) / size(cmp_str, 2) + 1);
                        id_ant(rem(id_ant,1) > 0) = []; % remove wrong matches
                        
                        id_ant = intersect(id_ant, t_check);
                        if ~isempty(id_ant)
                            if (nargin > 4)
                                log.addWarning(sprintf('"%s" antenna PCO/PCV found but with no radome', ant_name), log_lev);
                            else
                                log.addWarning(sprintf('"%s" antenna PCO/PCV found but with no radome', ant_name));
                            end
                        end
                    else
                        if isempty(id_ant)
                            if (nargin > 4)
                                log.addWarning(sprintf('"%s" antenna PCO/PCV not found', ant_name), log_lev);
                            else
                                log.addWarning(sprintf('"%s" antenna PCO/PCV not found', ant_name));
                            end
                        else
                            id_ant = intersect(id_ant, t_check);
                            if isempty(id_ant)
                                if (nargin > 4)
                                    log.addWarning(sprintf('"%s" antenna PCO/PCV found but it is out of validity range,\n      not using corrections', ant_name), log_lev);
                                else
                                    log.addWarning(sprintf('"%s" antenna PCO/PCV found but it is out of validity range,\n      not using corrections', ant_name));
                                end
                            else
                                id_ant = id_ant(1);
                            end
                        end
                        if ~isempty(id_ant)
                            ant(a) = this.list(id_ant);
                        end
                    end
                end
            end
        end
    end
    
    
end
