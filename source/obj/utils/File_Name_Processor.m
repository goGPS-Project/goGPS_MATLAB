%   CLASS File_Name_Processor
% =========================================================================
%
% DESCRIPTION
%   Class to process file names with keywords
%
% EXAMPLE
%   fnp = File_Name_Processor();
%   fnp.dateKeyRep('${WWWW}/igs${WWWWD}_${6H}.clk_30s', GPS_Time(datenum('2017/02/05 23:12:11')))
%
% FOR A LIST OF CONSTANTs and METHODS use doc FTP_Server
%
% REQUIRES:
%   goGPS settings;
%
% COMMENTS
%   Nothing to report

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0
%
%--------------------------------------------------------------------------
%  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
%  Written by:        Andrea Gatti
%  Contributors:      Andrea Gatti, Giulio Tagliaferro
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

classdef File_Name_Processor < handle

    properties (Constant)
        GPS_WEEK = '${WWWW}';
        GPS_WD = '${WWWWD}';
        GPS_DOW = '${D}';
        GPS_3H = '${3H}';
        GPS_6H = '${6H}';
        GPS_HH = '${HH}';
        GPS_RHH = '${rHH}'; % roundd hour
        GPS_YY = '${YY}';
        GPS_YYDOY = '${YYDOY}';
        GPS_YYYY = '${YYYY}';
        GPS_DOY = '${DOY}';
        GPS_FIRST_DOY_TRIMESTER = '${DOYT}';
        GPS_SESSION = '${S}';
        GPS_MM = '${MM}';
        GPS_DD = '${DD}';
        GPS_1D = '${1D}';
        GPS_QQ = '${QQ}';
        GPS_5M = '${5M}';
        GPS_MIN = '${mm}';
        GPS_SEC = '${ss}';
        VMF_RES = '${VMFR}';
        VMF_SOURCE = '${VMFS}';
    end

    properties (SetAccess = private, GetAccess = public)
        log = Logger.getInstance();
    end

    properties (SetAccess = protected, GetAccess = protected)
        addr;        % ip address (or alias)
        port = 21;   % connection port
        path = '/';  % path to files
    end

    methods
        function this = File_Name_Processor()
            % Constructor
        end

        function file_name_out = dateKeyRep(this, file_name, date, session, vmf_res, vmf_source)
            % substitute time placeholder with the proper format
            %
            % SYNTAX: 
            %     file_name = this.dateKeyRep(file_name, date, session)
            
            narginchk(3,6)
            if (nargin < 4)
                session = '';
            end
            if (nargin < 5)
                vmf_res = '';
            end
            if (nargin < 6)
                vmf_source = '';
            end
            file_name_out = file_name;
            if any(file_name_out == '$') && ~isempty(regexp(file_name_out, '\$\{Y|(DOY)', 'once'))
                [year, doy] = date.getDOY();
                year = sprintf('%04d', year);
                doy  = sprintf('%03d', doy);
                file_name_out = strrep(file_name_out, this.GPS_DOY, doy);
                file_name_out = strrep(file_name_out, this.GPS_FIRST_DOY_TRIMESTER, sprintf('%03d', GPS_Time.doy2fdoytr(doy,year)));
                
                file_name_out = strrep(file_name_out, this.GPS_YY, year(3:4));
                file_name_out = strrep(file_name_out, this.GPS_YYYY, year);
                file_name_out = strrep(file_name_out, this.GPS_YYDOY, [year(3:4) doy]);
            end
            if any(file_name_out == '$') && ~isempty(regexp(file_name_out, '\$\{(WWWW)|(WD)|(DOW)|(3H)|(6H)|(HH)|(QQ)|(5M)|(mm)|(ss)\}', 'once'))
                [gps_week, gps_sow, gps_dow] = date.getGpsWeek();
                file_name_out = strrep(file_name_out, this.GPS_WEEK, sprintf('%04d', gps_week(1)));
                if any(file_name_out == '$')
                    file_name_out = strrep(file_name_out, this.GPS_DOW, sprintf('%01d', gps_dow(1)));
                end
                if any(file_name_out == '$')
                    file_name_out = strrep(file_name_out, this.GPS_WD, sprintf('%04d%01d', gps_week(1), gps_dow(1)));
                    file_name_out = strrep(file_name_out, this.GPS_3H, sprintf('%02d', fix((gps_sow(1) - double(gps_dow(1)) * 86400)/(3*3600))*3));
                    file_name_out = strrep(file_name_out, this.GPS_6H, sprintf('%02d', fix((gps_sow(1) - double(gps_dow(1)) * 86400)/(6*3600))*6));
                    file_name_out = strrep(file_name_out, this.GPS_HH, sprintf('%02d', fix((gps_sow(1) - double(gps_dow(1)) * 86400)/(3600))));
                    file_name_out = strrep(file_name_out, this.GPS_RHH, sprintf('%02d', round((gps_sow(1) - double(gps_dow(1)) * 86400)/(3600))));
                    file_name_out = strrep(file_name_out, this.GPS_QQ, sprintf('%02d', mod(15 * fix((gps_sow(1) - double(gps_dow(1)) * 86400)/(900)), 60)));
                    file_name_out = strrep(file_name_out, this.GPS_5M, sprintf('%02d', mod(5 * fix((gps_sow(1) - double(gps_dow(1)) * 86400)/(300)), 60)));
                end
                if any(file_name_out == '$')
                    file_name_out = strrep(file_name_out, this.GPS_MIN, sprintf('%02d', mod(fix((gps_sow(1) - double(gps_dow(1)) * 86400)/(60)), 60)));
                    file_name_out = strrep(file_name_out, this.GPS_SEC, sprintf('%02d', mod(fix((gps_sow(1) - double(gps_dow(1)) * 86400)), 60)));
                end
            end
            if any(file_name_out == '$') && ~isempty(regexp(file_name_out, '\$\{(MM)|DD)|(1D)\}', 'once'))
                [date] = datevec(date.getMatlabTime());
                file_name_out = strrep(file_name_out, this.GPS_MM, sprintf('%02d', date(2)));
                file_name_out = strrep(file_name_out, this.GPS_DD, sprintf('%02d', date(3)));
                file_name_out = strrep(file_name_out, this.GPS_1D, sprintf('%d', date(3)));
            end            
            if any(file_name_out == '$') && ~isempty(session)
                file_name_out = strrep(file_name_out, this.GPS_SESSION, sprintf('%01d', session));
            end
            if any(file_name_out == '$') && ~isempty(vmf_res)
                file_name_out = strrep(file_name_out, this.VMF_RES, char(vmf_res));
            end
            if any(file_name_out == '$') && ~isempty(vmf_source)
                file_name_out = strrep(file_name_out, this.VMF_SOURCE, char(vmf_source));
            end
        end

        function step_sec = getStepSec(this, file_name)
            % get the shorter time keyword present in file_name [seconds]
            % SYNTEX: step_sec = this.getStepSec(file_name);
            % Check for GPS time placeholders
            step_sec = 0;
            if ~isempty(strfind(file_name, this.GPS_SEC))
                step_sec = 1;
            elseif ~isempty(strfind(file_name, this.GPS_MIN))
                step_sec = 60;
            elseif ~isempty(strfind(file_name, this.GPS_5M))
                step_sec = 300;
            elseif ~isempty(strfind(file_name, this.GPS_QQ))
                step_sec = 900;
            elseif ~isempty(strfind(file_name, this.GPS_HH))
                step_sec = 3600;
            elseif ~isempty(strfind(file_name, this.GPS_3H))
                step_sec = 3 * 3600;
            elseif ~isempty(strfind(file_name, this.GPS_6H))
                step_sec = 6 * 3600;
            elseif ~isempty((strfind(file_name, this.GPS_DOW))) || ~isempty(strfind(file_name, this.GPS_DD)) || ~isempty(strfind(file_name, this.GPS_WD)) || ~isempty(strfind(file_name, this.GPS_DOY)) || ~isempty(strfind(file_name, this.GPS_1D)) || ~isempty(strfind(file_name, this.GPS_YYDOY))
                step_sec = 24 * 3600;
            elseif ~isempty(strfind(file_name, this.GPS_WEEK))
                step_sec = 24 * 3600 * 7;
            elseif ~isempty(strfind(file_name, this.GPS_FIRST_DOY_TRIMESTER))
                step_sec = 24 * 3600 * 91;
            elseif ~isempty(strfind(file_name, this.GPS_YYYY)) || ~isempty(strfind(file_name, this.GPS_YY))
                step_sec = 24 * 3600 * 365;
            end
        end

        function [file_name_lst, date_list] = dateKeyRepBatch(this, file_name, date_start, date_stop, session_list, session_start, session_stop, vmf_res, vmf_source)
            % substitute time placeholder with the proper format
            % SYNTAX: 
            %   file_name = this.dateKeyRepBatch(file_name, date_start, date_stop)
            % NOTE: I consider only two possible file formats:
            %         - dependend on UTC time (year, doy)
            %         - dependent on GPS time (week, day of week, h of the day)

            if nargin < 5
                session_list = '0';
                session_start = '0';
                session_stop = '0';
            end
            if nargin < 8
                vmf_res = '';
                vmf_source = '';
            end
            session = ~isempty(strfind(file_name, this.GPS_SESSION));
            sss_start = strfind(session_list, session_start);
            sss_stop = strfind(session_list, session_stop);

            % Check for GPS time placeholders
            step_sec = this.getStepSec(file_name);

            if (step_sec > 0)
                file_name_lst = {};
                date_list = date_start.getCopy();
                if not(abs(step_sec - (24 * 3600 * 91)) < 1e-3)
                    date0 = GPS_Time((floor(((date_start.getNominalTime(0.5).getMatlabTime() - GPS_Time.GPS_ZERO) * 86400) / step_sec) * step_sec) / 86400 + GPS_Time.GPS_ZERO);
                else
                    date0 = date_start.getCopy();
                end
                date0 = date0.getNominalTime(60);
                date1 = date_stop.round(1);
                
                date_list.toUnixTime(); % keep an higher precision
                
                i = 1;
                % Find all the file in the interval of dates
                while (date0.getMatlabTime() <= date1.getMatlabTime())
                    tmp = this.dateKeyRep(file_name, date0,'',vmf_res,vmf_source);
                    if session
                        % run over session
                        for s = sss_start : length(session_list)
                            file_name_lst{i} = this.keyRep(tmp, this.GPS_SESSION, session_list(s)); %#ok<AGROW>
                            i = i + 1;
                        end
                        sss_start = 1;
                    else
                        file_name_lst{i} = tmp; %#ok<AGROW>
                        i = i + 1;
                    end
                    date0.addIntSeconds(step_sec);
                    date0 = date0.getNominalTime(60);
                    if (date0.getMatlabTime() <= date1.getMatlabTime())
                        date_list.append(date0);
                    end
                end
            else
                date_list = GPS_Time();
                i = 1;
                % run over session
                if session
                    for s = sss_start : length(session_list)
                        file_name_lst{i} = this.keyRep(file_name, this.GPS_SESSION, session_list(s)); %#ok<AGROW>
                        i = i + 1;
                    end
                else
                    if iscell(file_name)
                        file_name_lst = file_name;
                    else
                        file_name_lst = {file_name};
                    end
                end
            end

            if session
                sss_clean = length(file_name_lst) - (length(session_list)-sss_stop) + 1 : length(file_name_lst);
                file_name_lst(sss_clean) = [];
            end
            if iscell(file_name_lst)
                file_name_lst = file_name_lst';
            end
        end
    end

    methods (Static)
        function dir_path_out = getInnerDir(dir_path_list)
            % Get recursive dir list
            %
            % SYNTAX
            %   dir_path_out = getInnerDir(dir_path_list)
            
            dir_path_out = cell(0);
            % make base_dir a cell
            if ~iscell(dir_path_list)
                dir_path_list = {dir_path_list};
            end
            
            % loop on base_dirs
            for d = 1 : numel(dir_path_list)
                dir_path_out = [dir_path_out dir_path_list(d)];
                dir_path = dir_path_list{d};
                fprintf(' Scanning %s for subdirs\n', dir_path);
                % Search fo subdirs
                tmp = dir(dir_path);
                tmp = tmp([tmp.isdir]);
                dir_list = {};
                for i = numel(tmp) : -1 : 1
                    % remove hidden or special dirs (log)
                    if not(tmp(i).name(1) == '.') && not(strcmp(tmp(i).name, 'log'))
                        dir_list = [dir_list {fullfile(dir_path, tmp(i).name)}];
                    end
                end
                               
                if ~isempty(dir_list)                                    
                    for c = 1 : numel(dir_list)
                        cur_dir = dir_list(c);
                        dir_path_out = [dir_path_out File_Name_Processor.getInnerDir(cur_dir)];
                    end
                end
            end
        end
        
        function file_list_out = findFiles(dir_path_list, prefix)
            % Get recursive file list
            %
            % SYNTAX
            %   file_list_out = findFiles(dir_path_list, prefix)                        
            file_list_out = cell(0);
            % make base_dir a cell
            if ~iscell(dir_path_list)
                dir_path_list = {dir_path_list};
            end
            
            % loop on base_dirs
            for d = 1 : numel(dir_path_list)
                dir_path = dir_path_list{d};
                fprintf(' Scanning %s for subdirs\n', dir_path);
                % Search fo subdirs
                tmp = dir(dir_path);
                tmp_dir_list = tmp([tmp.isdir]);
                tmp_file_list = tmp(~[tmp.isdir]);
                dir_list = {};
                for i = 1 : numel(tmp_dir_list)
                    % remove hidden or special dirs (log)
                    if not(tmp_dir_list(i).name(1) == '.') && not(strcmp(tmp_dir_list(i).name, 'log'))
                        dir_list = [dir_list {fullfile(dir_path, tmp_dir_list(i).name)}];
                    end
                end
                               
                if ~isempty(dir_list)                                    
                    for c = 1 : numel(dir_list)
                        cur_dir = dir_list(c);
                        file_list_out = [file_list_out File_Name_Processor.findFiles(cur_dir, prefix)];
                    end
                end
                
                for i = 1 : numel(tmp_file_list)
                    if length(tmp_file_list(i).name) > length(prefix) &&  strcmp(tmp_file_list(i).name(1:length(prefix)), prefix)
                        file_list_out = [file_list_out {fullfile(dir_path, tmp_file_list(i).name)}];
                    end
                end
            end
        end
        
        function dir_path = getFullDirPath(dir_path, dir_base, dir_fallback, empty_fallback)
            % Get the full path given the relative one and the relative dir_base
            % It changes the folder with dir_fallback only if the changed folder exist
            %
            % SYNTAX: 
            %   dir_path = getFullDirPath(dir_path, <dir_base default = pwd>, fallback_dir_base);

            if nargin == 1
                dir_base = pwd;
            end
            if (nargin < 4)
                empty_fallback = dir_base;
            end
            if isempty(dir_base)
                dir_base = pwd;
            end
                       
            dir_path_bk = dir_path;

            fnp = File_Name_Processor;
            dir_path = fnp.checkPath(dir_path);
            if isunix
                if ~isempty(dir_path)
                    if (dir_path(1) == '/')
                        dir_path = fnp.checkPath(dir_path);
                    else
                        if ~isempty(dir_base)
                            if (dir_base(1) ~= '/')
                                dir_base = fnp.getFullDirPath(dir_base, pwd);
                            end
                        end
                        dir_path = fnp.checkPath([dir_base filesep dir_path]);
                    end
                else
                    dir_path = fnp.checkPath(empty_fallback);
                end
            else
                if length(dir_path) > 1
                    if (dir_path(2) == ':')
                        dir_path = fnp.checkPath(dir_path);
                    else
                        if length(dir_base) > 1
                            if (dir_base(2) == ':')
                                dir_base = fnp.getFullDirPath(dir_base, pwd);
                            end
                        end
                        if ~(numel(dir_path) >= 2 && strcmp(dir_path(1:2), '\\'))
                            dir_path = fnp.checkPath([dir_base filesep dir_path]);
                        end
                    end
                else
                    dir_path = fnp.checkPath(empty_fallback);
                end
            end
            
            dir_path = File_Name_Processor.remDirModifiers(dir_path);
                
            % Fallback if not exist
            if (nargin >= 3) && ~isempty(dir_fallback)
                if ~exist(dir_path, 'file')
                    prefix = '';
                    if numel(dir_fallback) > 2
                        prefix = iif(strcmp(dir_fallback(1:2),'\\'), '\\', '');
                    end
                    tmp = fnp.getFullDirPath(dir_path_bk, dir_fallback);
                    % change dir_path only if the new path exist
                    if exist(tmp, 'file')
                        dir_path = tmp;
                    end
                    dir_path = [prefix dir_path];
                end
            end
        end
        
        function path = remDirModifiers(path)
            % Utility to remove from a path dir changes with "." and ".."
            %
            % SYNTAX
            %   path = File_Name_Processor.remDirModifiers(path)
            
            if ~isempty(path)
                % remove './'
                path = strrep(path, [filesep '.' filesep], filesep);
                
                prefix = '';
                if numel(path) > 2
                    prefix = iif(strcmp(path(1:2),'\\'), '\\', '');
                end
                
                % extract sub folder names
                list = regexp(path,['[^' iif(filesep == '\', '\\', filesep) ']*'],'match');
                
                % search for "../"
                dir_up = find(strcmp(list,'..'));
                offset = 0;
                for i = dir_up
                    list(max(1, (i - 1 : i) - offset)) = [];
                    offset = offset + 2;
                end
                
                % restore full path start
                if isunix
                    path = [prefix filesep strCell2Str(list, filesep)];
                else
                    path = [prefix strCell2Str(list, filesep)];
                end
            end
        end

        function dir_path = getRelDirPath(dir_path, dir_base)
            % Get the full path given the relative one and the relative dir_base
            % SYNTAX: dir_path = getRelativeDirPath(dir_path, <dir_base default = pwd>);

            if nargin == 1
                dir_base = pwd;
            end

            fnp = File_Name_Processor;
            dir_base = fnp.getFullDirPath(dir_base, pwd);
            list_base = regexp(dir_base, ['[^' iif(filesep == '\', '\\', filesep) ']*'], 'match');
            n_dir_base = numel(list_base);

            if iscell(dir_path)
                for j = 1 : numel(dir_path)
                    dir_path{j} = fnp.getFullDirPath(dir_path{j}, dir_base);
                    drive_path = regexp(dir_path{j},'.(?=\:)', 'match', 'once');
                    drive_base = regexp(dir_base,'.(?=\:)', 'match', 'once');
                    if (isempty(drive_path) || drive_path == drive_base)
                        list_path = regexp(dir_path{j}, ['[^' iif(filesep == '\', '\\', filesep) ']*'], 'match');
                        n_dir_path = numel(list_path);
                        i = 0;
                        while (i < n_dir_base) && (i < n_dir_path) && strcmp(list_base{i+1}, list_path{i+1})
                            i = i + 1;
                        end
                        if (n_dir_base -i) > 0
                            list_path = [{repmat(['..' filesep],1, n_dir_base - i)} list_path(i+1:end)];
                        else
                            list_path = list_path(i+1:end);
                        end
                        
                        dir_path{j} = strrep(strCell2Str(list_path, filesep), [filesep filesep], filesep);
                    end
                end
            else
                % if dir_path is a network address && dir_base is not
                prefix = '';
                if (numel(dir_path) > 2 && strcmp(dir_path(1:2),'\\')) && ~((numel(dir_base) > 2 && strcmp(dir_base(1:2),'\\')))
                    n_dir_base = 0;
                    list_base = {};
                    prefix = '\\';
                end
                dir_path = fnp.getFullDirPath(dir_path, dir_base);
                drive_path = regexp(dir_path,'.(?=\:\\)', 'match', 'once');
                drive_base = regexp(dir_base,'.(?=\:\\)', 'match', 'once');
                if (isempty(drive_path) || drive_path == drive_base)
                    list_path = regexp(dir_path, ['[^' iif(filesep == '\', '\\', filesep) ']*'], 'match');
                    n_dir_path = numel(list_path);
                    i = 0;
                    while (i < n_dir_base) && (i < n_dir_path) && strcmp(list_base{i+1}, list_path{i+1})
                        i = i + 1;
                    end
                    if (n_dir_base -i) > 0
                        if i > 0
                            list_path = [{repmat(['..' filesep],1, n_dir_base - i)} list_path(i+1:end)];                            
                        elseif isempty(prefix)
                            prefix = filesep;
                        end
                        
                    else
                        list_path = list_path(i+1:end);
                    end
                    
                    dir_path = [prefix strrep(strCell2Str(list_path, filesep), [filesep filesep], filesep)];
                end
            end
        end

        function file_name = getFileName(file_name)
            % Get only the file name of a full path
            [~, name, extension] = fileparts(file_name);
            file_name = [name extension];
        end
        
        function path = getPath(file_name)
            % Get only the path up to the folder name from a full path
            [path, ~, ~] = fileparts(file_name);
        end
        
        function file_name = keyRep(file_name, key, substitution)
            % Substitute a key in the file_name with another value
            file_name = strrep(file_name, key, substitution);
        end

        function [str_cell] = toIniStringComment(str_cell)
            % Get the list of accepted substitution
            if (nargin == 0)
                str_cell = {};
            end
            str_cell = Ini_Manager.toIniStringComment('Special Keywords that can be used in file names:', str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf(' - %s    4 char GPS week', File_Name_Processor.GPS_WEEK), str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf(' - %s 4+1 char GPS week + day of the week', File_Name_Processor.GPS_WD), str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf(' - %s       1 char day of the week', File_Name_Processor.GPS_DOW), str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf(' - %s      2 char GPS hour (00, 03, 06, 09, 12, 15, 18, 21)', File_Name_Processor.GPS_3H), str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf(' - %s      2 char GPS hour (00, 06, 12, 18)', File_Name_Processor.GPS_6H), str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf(' - %s      2 char GPS hour', File_Name_Processor.GPS_HH), str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf(' - %s      2 char GPS quarter of hour (00, 15, 30, 45)', File_Name_Processor.GPS_QQ), str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf(' - %s      2 char GPS five minutes (05, 10, ... , 55)', File_Name_Processor.GPS_5M), str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf(' - %s 2+3 char GPS year + day of year', File_Name_Processor.GPS_YYDOY), str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf(' - %s    4 char GPS year', File_Name_Processor.GPS_YYYY), str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf(' - %s      2 char GPS year', File_Name_Processor.GPS_YY), str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf(' - %s      2 char GPS month', File_Name_Processor.GPS_MM), str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf(' - %s      2 char GPS day', File_Name_Processor.GPS_DD), str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf(' - %s     3 char GPS day of the year', File_Name_Processor.GPS_DOY), str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf(' - %s       1 char session', File_Name_Processor.GPS_SESSION), str_cell);
        end

        function full_path = getFullPath(dir_path, file_path)
            % Get the list of files combining dir_path + file_path
            % It works when file_path is a cell or a simple string
            % SYNTAX: full_path = IO_Settings.getFullPath(dir_path, file_path);

            fnp = File_Name_Processor();

            if iscell(file_path)
                full_path = cell(numel(file_path),1);
                for f = 1 : numel(file_path)
                    full_path{f} = fnp.checkPath(fullfile(dir_path, file_path{f}));
                end
            else
                full_path = fnp.checkPath(fullfile(dir_path, file_path));
            end
        end

        function [universal_path, is_valid] = checkPath(path, prj_home)
            % Conversion of path OS specific to a universal one: "/" or "\" are converted to filesep
            %
            % SYNTAX:
            %   universal_path = checkPath(path)
            %
            % INPUT:
            %   path
            %
            % OUTPUT:
            %   universal_path
            %   < is_valid >            optional, contains the status of existence
            %
            % DESCRIPTION:
            %   Conversion of path OS specific to a universal one: "/" or "\" are converted to filesep
            %   if the second parameter is present is_valid contains the status of existence
            %       2 => is a file
            %       7 => is a folder

            if instr(path, '${PRJ_HOME}')
                if nargin < 2 
                %    state = Core.getState();
                %    prj_home = state.getHomeDir;
                    prj_home = '';
                end
                path = strrep(path, '${PRJ_HOME}', prj_home);
            end
            
            prefix = '';
            if numel(path) > 2
                prefix = iif(strcmp(path(1:2),'\\'), '\', '');
            end
            
            if not(isempty(path))
                if (iscell(path))
                    % for each line of the cell
                    universal_path = cell(size(path));
                    for c = 1 : length(path)
                        check_drive_win = strfind(path{c},':');
                        if (~isempty(check_drive_win))
                            path{c} = path{c}(max(1, check_drive_win - 1) : end);
                        end
                        universal_path{c} = regexprep(path{c}, '(\\(?![ ]))|(\/)', filesep);
                        universal_path{c} = regexprep(universal_path{c}, ['\' filesep '\' filesep], filesep);
                    end
                    if (nargout == 2)
                        is_valid = zeros(size(path));
                        for c = 1 : length(path)
                            is_valid(c) = exist(universal_path{c}, 'file'); % if it is a file is_valid contains 2, if it is a dir it contains 7
                        end
                    end
                else
                    if isstring(path)
                        path = char(path);
                    end
                    if ischar(path)
                        check_drive_win = strfind(path,':');
                        if (~isempty(check_drive_win))
                            path = path(max(1, check_drive_win - 1) : end);
                        end
                        universal_path = regexprep(path, '(\\(?![ ]))|(\/)', filesep);
                        universal_path = regexprep(universal_path, ['\' filesep '\' filesep], filesep);
                        if (nargout == 2)
                            is_valid = exist(path, 'file'); % if it is a file is_valid contains 2, if it is a dir it contains 7
                        end
                    else
                        universal_path = '';
                        is_valid = 0;
                    end
                end
            else
                universal_path = '';
                is_valid = 0;
            end
            universal_path = [prefix universal_path];
        end

    end
end