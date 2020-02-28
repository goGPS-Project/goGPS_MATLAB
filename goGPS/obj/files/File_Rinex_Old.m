%   CLASS File_Rinex_Old
% =========================================================================
%
% DESCRIPTION
%   Class to store file_paths for RINEX files
%
% EXAMPLE
%   fr = File_Rinex_Old(file_name, verbosity_lev);
%
% FOR A LIST OF CONSTANTS and METHODS use doc File_Rinex_Old

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b7
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

classdef File_Rinex_Old < Exportable

    properties (SetAccess = protected, GetAccess = protected)
        log = Logger.getInstance(); % Handler to the log object
    end

    properties (SetAccess = protected, GetAccess = public)
        is_valid = false;                            % flag, if true it means that the object contains at least one valid rinex file
        base_dir = {'../data/project/default_DD/RINEX/'};           % directory containing all the files
        file_name_list = {'yamatogawa_rover', 'yamatogawa_master'}; % file names (they can be multiple files for different days)
        ext = {'.obs', '.obs'};                                     % file names extension (they can be multiple files for different days)
        is_valid_list = false(1, 2);                 % for each element of file_name_list check the validity of the file

        is_composed = false;                         % when this flag is set, it means that the file_name depends on variables such as DOY DOW YYYY SSSS MM ecc...

        marker_name = {};                            % marker name of the files
        coo = Coordinates.fromXYZ([0 0 0])           % receiver coordinates
        first_epoch = GPS_Time();                    % first epoch stored in the RINEX (updated after checkValidity)
        last_epoch = GPS_Time();                     % last epoch stored in the RINEX (updated after checkValidity)

        verbosity_lev = 0;                           % Level of verbosity  0 = all messages, see Logger for more information (if negative) no messages at all

        eoh = 0;                                     % end of header line - store the last line of the header
    end

    properties (SetAccess = private, GetAccess = private)
        id_date = 2:28;                              % Last character containing the 6 fields of the date (in case of 4digits year), it the pends on the type of rinex files (28 -> OBS RINEX 3)
    end

    methods
        function this = File_Rinex_Old(file_name, verbosity_lev, flag_header_only)
            % Creator of File_Rinex_Old simple parser            
            %
            % INPUT 
            %   file_name        file path (can be a cell array)
            %   verbosity_lev    verbosity level for logging purposes
            %                    (default = 100)
            %   flag_header_only if true and first/last epochs are not found 
            %                    in header do not search them in the file
            %                    (default = false)
            %                    this option may leave the obj corrupted (use with care)
            %
            % SYNTAX
            %   File_Rinex_Old (file_name , verbosity_level, flag_header_only)
            %
            if nargin < 3 || isempty(flag_header_only)
                flag_header_only = false;
            end
            
            if nargin == 0
                % Empty File_Rinex_Old;
            else
                % fill the path with the imported file names
                if ~iscellstr(file_name)
                    file_name = {file_name};
                end
                this.file_name_list = {};
                this.ext = {};
                if nargin == 1 || isempty(verbosity_lev)
                    verbosity_lev = 100;
                end
                
                if verbosity_lev < 50 && numel(file_name) > 20
                    w_bar = Go_Wait_Bar.getInstance(numel(file_name), 'Checking rinex files');
                    w_bar.createNewBar
                    for f = 1 : numel(file_name)
                        [this.base_dir{f}, this.file_name_list{f}, this.ext{f}] = fileparts(checkPath(file_name{f}));
                        w_bar.go();
                    end
                    w_bar.close();
                else
                    for f = 1 : numel(file_name)
                        [this.base_dir{f}, this.file_name_list{f}, this.ext{f}] = fileparts(checkPath(file_name{f}));
                    end
                end
                
                if nargin >= 2
                    this.verbosity_lev = verbosity_lev;
                end
                
                this.checkValidity(flag_header_only);
            end
        end
        
        function copy = getCopy(this)
            for r = numel(this):-1:1
                copy(r) = File_Rinex_Old();
                copy(r).copyFrom(this(r));
            end
        end
        
        function copyFrom(this, File_Rinex_Old)
            % Copy from an object of the same type
            %
            % SYNTAX
            %   this.copyFrom(time)
            this.is_valid       = File_Rinex_Old.is_valid;
            this.base_dir       = File_Rinex_Old.base_dir;
            this.file_name_list = File_Rinex_Old.file_name_list;
            this.ext            = File_Rinex_Old.ext;
            this.is_valid_list  = File_Rinex_Old.is_valid_list ;
            this.is_composed    = File_Rinex_Old.is_composed;
            this.first_epoch    = File_Rinex_Old.first_epoch;
            this.last_epoch     = File_Rinex_Old.last_epoch;
            this.verbosity_lev  = File_Rinex_Old.verbosity_lev;
            this.eoh            = File_Rinex_Old.eoh;
        end
    end

    methods
        function checkValidity(this, flag_header_only)
            % Update the status of validity of the files here pointed
            %
            % SYNTAX
            %   this.checkValidity();

            % pre allocate
            this.is_valid_list = false(1, numel(this.file_name_list));
            this.eoh = zeros(1, numel(this.file_name_list));
            this.first_epoch = GPS_Time();
            this.last_epoch = GPS_Time();
            this.coo = Coordinates;           % receiver coordinates
            % for each file present in the list
            for f = 1 : numel(this.file_name_list)
                marker_name = '';
                full_path = fullfile(this.base_dir{f}, [this.file_name_list{f} this.ext{f}]);

                % check the existence
                this.is_valid_list(f) = exist(full_path, 'file') == 2;
                if this.is_valid_list(f)
                    % try to find the first and the last epoch stored in the file
                    try
                        %%
                        fid = fopen(fullfile(this.base_dir{f}, [this.file_name_list{f} this.ext{f}]));
                        l = 1;
                        line = fgetl(fid);
                        date_start = '';
                        date_stop = '';
                        coo = '';
                        while isempty(strfind(line,'END OF HEADER')) && ischar(line) %#ok<*STREMP>
                            l = l + 1;
                            line = fgetl(fid);
                            if numel(line) > 70 && isempty(marker_name) && line(64) == 'K' % read marker_name
                                marker_name = strtrim(regexp(line, '.*(?=MARKER NAME)', 'match', 'once'));
                            end
                            if numel(line) > 76
                                if strcmp(line(61:64), 'TIME')
                                    tmp = regexp(line, '.*(?=GPS         TIME OF FIRST OBS)', 'match', 'once');
                                    if ~isempty(tmp)
                                        date_start = strtrim(tmp);
                                    end
                                    tmp = regexp(line, '.*(?=GPS         TIME OF LAST OBS)', 'match', 'once');
                                    if ~isempty(tmp)
                                        date_stop = strtrim(tmp);
                                    end
                                else                                    
                                    tmp = regexp(line, '.*(?=APPROX POSITION XYZ)', 'match', 'once');
                                    if ~isempty(tmp)
                                        coo = strtrim(tmp);
                                    end                                    
                                end                            
                            end
                        end
                        
                        if ~isempty(coo)
                            this.coo.append(Coordinates.fromStringXYZ(coo));
                        else
                            this.coo.append(Coordinates.fromXYZ([0 0 0]));
                        end
                        this.eoh(f) = l;
                        
                        % If I did not found date of start and date of end in the file header
                        % Try to read them
                        if ~isempty(date_start)
                            this.first_epoch.addEpoch(date_start, [], true);
                        else
                            epoch_line = fgetl(fid);
                            if epoch_line == -1
                                error('Empty file');
                            end
                            % try to guess the time format
                            [id_start, id_stop] = regexp(epoch_line, '[.0-9]*');
                            this.id_date = id_start(1) : id_stop(6); % save first and last char limits of the date in the line -> suppose it composed by 6 fields
                            
                            this.first_epoch.addEpoch(epoch_line(this.id_date), [], true);
                        end
                        this.log.addStatusOk(['"' this.file_name_list{f} this.ext{f} '" appears to be a valid RINEX'], this.verbosity_lev);
                        this.log.addMessage(sprintf('        first epoch found at: %s', this.first_epoch.last.toString()), this.verbosity_lev);
                        
                        if ~isempty(date_stop)
                            this.last_epoch.addEpoch(date_stop, [], true);
                            this.log.addMessage(sprintf('        last  epoch found at: %s', this.last_epoch.last.toString()), this.verbosity_lev);
                        else
                            if flag_header_only
                                this.log.addWarning('Last epoch not found in header, search in file is not enabled\nThe last epoch has not been saved within the  object', this.verbosity_lev);
                            else
                                % go to the end of the file to search for the last epoch
                                % to be sure to find at least one line containing a valid epoch, go to the end of the file minus 5000 characters
                                fseek(fid, -10000, 'eof');
                                fgetl(fid); % Probably i'm not at the beginning of a line -> disregard the first reading
                                % Start searching for a valid epoch
                                line = fgetl(fid);
                                time = [];
                                loop_n = 1;
                                while isempty(time) && ischar(line)
                                    while ischar(line)
                                        % An epoch line has the second character containing the year of the observation
                                        % e.g. RINEX 2:     " 15  8 23  0  0  0.0000000  0  8G05G07G28G02G06G09G30G13"
                                        %      RINEX 2 NAV: "G01 2006 10 01 00 00 00 0.798045657575E-04 0.227373675443E-11 0.000000000000E+00"
                                        %      RINEX 3:     "> 2016  7 18  0  1  0.0000000  0 36"
                                        %      RINEX 3 NAV: " 3 98  2 15  0 15  0.0 0.163525342941D-03 0.363797880709D-11 0.108000000000D+05"
                                        %      RINEX 3 MET: " 1996  4  1  0  0 15  987.1   10.6   89.5"
                                        % this check could not work when comment are present after the header
                                        if (numel(line) > 20) && ~isempty(regexp(line(1:15),'( [0-9]{1,4} [ 0-9]{1}[0-9]{1} [ 0-9]{1}[0-9]{1})', 'once'))
                                            % Check that the read epoch is within 30 days from the first epoch
                                            % (if it's further in time it's probably a misleading false epoch line)
                                            time = GPS_Time(line(this.id_date), [], true);
                                            time_diff = ((time - this.first_epoch.last())/86400);
                                            if (time_diff < 30) && (time_diff > 0)
                                                epoch_line = line;
                                            end
                                        end
                                        line = fgetl(fid);
                                    end
                                    loop_n = loop_n + 1;
                                    fseek(fid, loop_n * -10000, 'eof'); % If no valid time have been found try to go back more...
                                    line = fgetl(fid);
                                end
                                this.last_epoch.addEpoch(epoch_line(this.id_date), [], true);
                                this.log.addMessage(sprintf('        last  epoch found at: %s', this.last_epoch.last.toString()), this.verbosity_lev);
                            end
                        end
                        fclose(fid);
                        this.is_valid_list(f) = true;                        
                    catch ex
                        this.log.addWarning(['"' this.file_name_list{f} this.ext{f} '" appears to be a corrupted RINEX file'], this.verbosity_lev);
                        this.is_valid_list(f) = false;
                    end
                else
                    this.log.addError(['"' this.file_name_list{f} this.ext{f} '" appears to be missing'], this.verbosity_lev);
                    this.is_valid_list(f) = false;
                end
                % store marker_name
                if isempty(marker_name)
                    this.marker_name{f} = upper(this.file_name_list{f}(1:min(4, numel(this.file_name_list{f}))));
                else
                    this.marker_name{f} = marker_name;
                end
            end
            this.is_valid = all(this.is_valid_list);
            if (~this.is_valid)
                this.log.addWarning('Some or all the RINEX files are corrupted or missing!!!', this.verbosity_lev);
            end
        end

        function file_name = getFileName(this, file_number)
            % Get the full path of a file (if the object contains a list of files, the id can be specified)
            % SYNTAX: file_name = this.getFileName(<file_number = 1>)
            if nargin == 1
                file_number = 1;
            end
            file_name = fullfile(this.base_dir{file_number}, [this.file_name_list{file_number} this.ext{file_number}]);
        end
        
        function first_epoch = getFirstEpoch(this, session)
            % Get the first epoch of a session
            %
            % SYNTAX
            %   first_epoch = this.getFirstEpoch(session)
            if nargin == 1
                if any(this.is_valid_list)
                    first_epoch = this.first_epoch.getCopy;
                else
                    first_epoch = GPS_Time();
                end
            else
                if Core.getState.isRinexSession
                    if this.is_valid_list(session)
                        id_ss = cumsum(this.is_valid_list);
                        first_epoch = this.first_epoch.getEpoch(id_ss(session));
                    else
                        first_epoch = GPS_Time();
                    end
                else
                    first_epoch = this.state.getSessionLimits(session).first;
                end
            end
        end

        function last_epoch = getLastEpoch(this, session)
            % Get the last epoch of a session
            %
            % SYNTAX
            %   first_epoch = this.getFirstEpoch(session)
            if nargin == 1
                if any(this.is_valid_list)
                    last_epoch = this.last_epoch.getCopy;
                else
                    last_epoch = GPS_Time();
                end
            else
                if Core.getState.isRinexSession
                    if this.is_valid_list(session)
                        id_ss = cumsum(this.is_valid_list);
                        last_epoch = this.last_epoch.getEpoch(id_ss(session));
                    else
                        last_epoch = GPS_Time();
                    end
                else
                    last_epoch = this.state.getSessionLimits(session).last;
                end
            end
        end
        
        function line_num = getEOH(this, file_number)
            % Get the end of header line of a RINEX file (if the object contains a list of files, the id can be specified)
            % SYNTAX: line_num = this.getEOH(<file_number = 1>)
            if nargin == 1
                file_number = 1;
            end
            line_num = this.eoh(file_number);
        end

        function validity = isValid(this, file_number)
            % Get the validity of a RINEX file or the object (if the object contains a list of files, the id can be specified)
            % SYNTAX: validity = isValid(<file_number>)
            validity = false;
            for r = 1 : numel(this)
                
                if (nargin == 1)
                    validity = validity || any([this(r).is_valid_list]);
                else
                    validity = validity || this(r).is_valid_list(file_number);
                end
            end
        end

        function keepFiles(this, date_start, date_stop)
            % keep only files contained in the two dates
            %
            % SYNTAX
            % this.keepFiles(date_start, date_stop)
            for r = 1 : numel(this)
                [in_bound] = Core_Utils.timeIntersect(this(r).first_epoch, this(r).last_epoch, date_start, date_stop);
                to_keep = this(r).is_valid_list;
                to_keep(to_keep) = in_bound;
                this(r).keep(to_keep);
            end
        end
        
        function has_obs = hasObsInSession(this, date_start, date_stop)
            % has observation in the two dates
            %
            % SYNTAX
            % this.keepFiles(date_start, date_stop)
            has_obs = false(numel(this),1);
            for r = 1 : numel(this)
                if this(r).isValid
                    in_bound = Core_Utils.timeIntersect(this(r).first_epoch, this(r).last_epoch, date_start, date_stop);
                    has_obs(r) = any(in_bound);
                end
            end
        end
        
        function keep(this,id)
            % keep only element specified by id
            %
            % SYNTAX
            % this.keep(id)
            this.base_dir = this.base_dir(id);
            this.file_name_list = this.file_name_list(id);
            this.ext = this.ext(id);
            this.eoh = this.eoh(id);
            this.first_epoch = this.first_epoch.getEpoch(id(this.is_valid_list));
            this.last_epoch  = this.last_epoch.getEpoch(id(this.is_valid_list));
            this.is_valid_list = this.is_valid_list(id);
        end
    end
end
