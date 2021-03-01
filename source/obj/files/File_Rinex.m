%   CLASS File_Rinex
% =========================================================================
%
% DESCRIPTION
%   Class to store file_paths for RINEX files
%
% EXAMPLE
%   fr = File_Rinex(file_name, verbosity_lev);
%
% FOR A LIST OF CONSTANTS and METHODS use doc File_Rinex

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 Geomatics Research & Development srl (GReD)
%  Written by:        Andrea Gatti
%  Contributors:      Andrea Gatti, Giulio Tagliaferro
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

classdef File_Rinex < Exportable

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
        
        trck_availability                            % boolean to get tracking availabilty the tracking are in Core_sky.group_delays_flags

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
        function this = File_Rinex(file_name, verbosity_lev, flag_header_only)
            % Creator of file_rinex simple parser            
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
            %   File_Rinex (file_name , verbosity_level, flag_header_only)
            %
            if nargin < 3 || isempty(flag_header_only)
                flag_header_only = false;
            end
            
            if nargin == 0
                % Empty File_Rinex;
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
                copy(r) = File_Rinex();
                copy(r).copyFrom(this(r));
            end
        end
        
        function copyFrom(this, file_rinex, id)
            % Copy from an object of the same type
            %
            % SYNTAX
            %   this.copyFrom(time)
            if nargin < 3 || isempty(id)
                id = 1 : numel(file_rinex.file_name_list);
            end
            % Get id within coo and time limits objects
            valid_lid = false(size(file_rinex.is_valid_list));
            valid_lid(id) = true;
            valid_lid = valid_lid & file_rinex.is_valid_list;
            valid_id = cumsum(file_rinex.is_valid_list);
            valid_id = valid_id(valid_lid);
            % Extract fields to copy
            this.is_valid       = all(file_rinex.is_valid_list(id));
            this.base_dir       = file_rinex.base_dir(id);
            this.file_name_list = file_rinex.file_name_list(id);
            this.ext            = file_rinex.ext(id);
            this.is_valid_list  = file_rinex.is_valid_list(id);
            this.is_composed    = file_rinex.is_composed;
            this.first_epoch    = file_rinex.first_epoch.getEpoch(valid_id);
            this.last_epoch     = file_rinex.last_epoch.getEpoch(valid_id);
            this.verbosity_lev  = file_rinex.verbosity_lev;
            this.eoh            = file_rinex.eoh(id);
            this.coo            = file_rinex.coo.getCopy();            
        end
    end

    methods
        function checkValidity(this, flag_header_only)
            % Update the status of validity of the files here pointed
            % Buffered read
            %
            % SYNTAX
            %   this.checkValidity();

            % pre allocate
            this.is_valid_list = false(1, numel(this.file_name_list));
            this.eoh = zeros(1, numel(this.file_name_list));
            this.first_epoch = GPS_Time();
            this.last_epoch = GPS_Time();
            this.coo = Coordinates;           % receiver coordinates
            log = Core.getLogger;
            % for each file present in the list
            for f = 1 : numel(this.file_name_list)
                marker_name = '';
                full_path = fullfile(this.base_dir{f}, [this.file_name_list{f} this.ext{f}]);

                % try to find the first and the last epoch stored in the file
                try
                    %%
                    this.is_valid_list(f) = true;
                    fid = fopen(full_path, 'rt');
                    if fid < 0
                        log.addError(['"' this.file_name_list{f} this.ext{f} '" appears to be missing'], this.verbosity_lev);
                        this.is_valid_list(f) = false;
                    else
                        buf = fread(fid, 1e4, '*char')';
                        if length(buf) > 65 && strcmp(buf(61:64), 'CRIN')
                            log.addError(sprintf('Check the following file, it seems to be hatanaka compressed\nDecompress "%s"', full_path));
                            this.is_valid_list(f) = false;
                            fclose(fid);
                        else
                            l = 1;
                            % read by buffer 10k char at a time
                            % fseek(fid, 0, 'bof'); buf = fread(fid, 1e4, '*char')';
                            
                            % detect windows carriage return
                            if ~isempty(find(buf(1:min(1000,numel(buf))) == 13, 1, 'first'))
                                has_cr = true;  % The file has carriage return - I hate you Bill!
                            else
                                has_cr = false;  % The file is UNIX standard
                            end
                            
                            % get new line separators
                            nl = regexp(buf, '\n')';
                            if nl(end) <  (numel(buf) - double(has_cr))
                                nl = [nl; numel(buf)];
                            end
                            lim = [[1; nl(1 : end - 1) + 1] (nl - 1 - double(has_cr))];
                            lim = [lim lim(:,2) - lim(:,1)];
                            while lim(end,3) < 3
                                lim(end,:) = [];
                            end
                            
                            % Check for meteo rinex
                            flag_met = buf(21) == 'M';
                            
                            line = buf(lim(l,1) : lim(l,2));
                            date_start = '';
                            trck_name = Core_Sky.group_delays_flags;
                            trck_availability = false(size(trck_name,1),1);
                            cur_trck_sys = '';
                            date_stop = '';
                            coo = '';
                            eof = false;
                            par_to_find = 4;
                            %while (par_to_find > 0) && isempty(strfind(line,'END OF HEADER')) && not(eof)
                            while ~((length(line) > 61) && (line(61) == 'E')) && not(eof)
                                l = l + 1;
                                if l >= size(lim, 1)
                                    % out of buffer
                                    % read block:
                                    buf = [buf char(fread(fid, 1e4))'];
                                    
                                    % get new line separators
                                    nl = regexp(buf, '\n')';
                                    if nl(end) <  (numel(buf) - double(has_cr))
                                        nl = [nl; numel(buf)];
                                    end
                                    lim = [[1; nl(1 : end - 1) + 1] (nl - 1 - double(has_cr))];
                                    lim = [lim lim(:,2) - lim(:,1)];
                                    while lim(end,3) < 3
                                        lim(end,:) = [];
                                    end
                                end
                                if l > size(lim, 1)
                                    eof = true;
                                else
                                    line = buf(lim(l,1) : lim(l,2));
                                    if par_to_find > 0
                                        if numel(line) >= 70 && isempty(marker_name) && line(64) == 'K' % read marker_name
                                            marker_name = strtrim(regexp(line, '.*(?=MARKER NAME)', 'match', 'once'));
                                            par_to_find = par_to_find - 1;
                                        end
                                        if numel(line) >= 76
                                            if line(61) == 'T'
                                                %tmp = regexp(line, '.*(?=GPS         TIME OF FIRST OBS)', 'match', 'once');
                                                if line(69) == 'F'
                                                    date_start = strtrim(line(1:48));
                                                    par_to_find = par_to_find - 1;
                                                elseif line(69) == 'L'
                                                    %tmp = regexp(line, '.*(?=GPS         TIME OF LAST OBS)', 'match', 'once');
                                                    date_stop = strtrim(line(1:48));
                                                    par_to_find = par_to_find - 1;
                                                end
                                            elseif line(76) == 'Y' || (numel(line) >= 79 && line(79) == 'V') %  SYS / # / OBS TYPES || # / TYPES OF OBSERV
                                                if line(1) ~= ' '
                                                    cur_trck_sys = line(1);
                                                else
                                                    cur_trck_sys = ' ';
                                                end
                                                trcks = strsplit(line(8:60),' ');
                                                for t = trcks
                                                    if ~isempty(t{1})
                                                        t = t{1};
                                                        t = [t repmat(' ',1,3-length(t))];
                                                        idx = trck_name(:,1) == cur_trck_sys & trck_name(:,2) == t(1) & trck_name(:,3) == t(2) & trck_name(:,4) == t(3);
                                                        trck_availability(idx) = true;
                                                        if strcmp(t,'P1 ')% rinex 2 codes does not specify contellation, assuming GPS
                                                            trck_availability(strLineMatch(trck_name,'GC1W')) = true;
                                                        elseif strcmp(t,'P2 ')
                                                            trck_availability(strLineMatch(trck_name,'GC2W')) = true;
                                                        elseif strcmp(t,'C1 ')
                                                            trck_availability(strLineMatch(trck_name,'GC1C')) = true;
                                                        elseif strcmp(t,'C2 ')
                                                            trck_availability(strLineMatch(trck_name,'GC2C')) = true;
                                                        end
                                                    end
                                                end
                                            elseif line(66) == 'Y' % "# / TYPES OF OBSERV"
                                                % temporary solution for RINEX 2
                                                % to be improved
                                                trcks = strsplit(line(8:60),' ');
                                                for cur_trck_sys = Core.getConstellationCollector.getActiveSysChar
                                                    for t = trcks
                                                        if ~isempty(t{1})
                                                            t = t{1};
                                                            if numel(t) == 2 % correct missing character for RINEX 2
                                                                % It should be done multi constellation
                                                                % GPS C1 -> C1C
                                                                % GPS C2 -> C2C
                                                                t = [t 'C'];
                                                            end
                                                            idx = trck_name(:,1) == cur_trck_sys & trck_name(:,2) == t(1) & trck_name(:,3) == t(2) & trck_name(:,4) == t(3);
                                                            trck_availability(idx) = true;
                                                        end
                                                    end
                                                end
                                            else
                                                %tmp = regexp(line, '.*(?=APPROX POSITION XYZ)', 'match', 'once');
                                                % character to recognize approximate position for met file: 'E' => sensor pos
                                                % character to recognize approximate position for met file: 'P' => approx position
                                                pos_char = iif(flag_met, 'E', 'P');
                                                if line(62) == pos_char % APPROXIMATE POSITION
                                                    coo = strtrim(line(1:43));
                                                    par_to_find = par_to_find - 1;
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                            
                            eoh = l; %end of header
                            if ~isempty(trck_availability)
                                this.trck_availability = trck_availability;
                            end
                            
                            % If I did not found date of start and date of end in the file header
                            % Try to read them
                            if ~isempty(date_start)
                                this.first_epoch.addEpoch(date_start, [], true);
                            else
                                l = l + 1;
                                if l > size(lim, 1)
                                    % out of buffer
                                    % read block:
                                    buf = [buf char(fread(fid, 1e4))'];
                                    
                                    % get new line separators
                                    nl = regexp(buf, '\n')';
                                    if nl(end) <  (numel(buf) - double(has_cr))
                                        nl = [nl; numel(buf)];
                                    end
                                    lim = [[1; nl(1 : end - 1) + 1] (nl - 1 - double(has_cr))];
                                    lim = [lim lim(:,2) - lim(:,1)];
                                    while lim(end,3) < 3
                                        lim(end,:) = [];
                                    end
                                end
                                if l <= size(lim,1)
                                    epoch_line = buf(lim(l,1) : lim(l,2));
                                    % try to guess the time format
                                    [id_start, id_stop] = regexp(epoch_line, '[.0-9]*');
                                    this.id_date = id_start(1) : id_stop(6); % save first and last char limits of the date in the line -> suppose it composed by 6 fields
                                    this.first_epoch.addEpoch(epoch_line(this.id_date), [], true);
                                else
                                    error(sprintf('"%s" seems corrupted', full_path));
                                end
                            end
                            log.addStatusOk(['"' this.file_name_list{f} this.ext{f} '" appears to be a valid RINEX'], this.verbosity_lev);
                            log.addMessage(sprintf('        first epoch found at: %s', this.first_epoch.last.toString()), this.verbosity_lev);
                            
                            if ~isempty(date_stop)
                                this.last_epoch.addEpoch(date_stop, [], true);
                                log.addMessage(sprintf('        last  epoch found at: %s', this.last_epoch.last.toString()), this.verbosity_lev);
                            else
                                if flag_header_only
                                    log.addWarning('Last epoch not found in header, search in file is not enabled\nThe last epoch has not been saved within the object', this.verbosity_lev);
                                    if ~isempty(this.first_epoch.last.getMatlabTime)
                                        this.last_epoch.addEpoch(this.first_epoch.last.getMatlabTime, [], true);
                                    end
                                else
                                    % go to the end of the file to search for the last epoch
                                    % to be sure to find at least one line containing a valid epoch, go to the end of the file minus 10000 characters
                                    fseek(fid, 0, 'eof');
                                    f_size = ftell(fid);
                                    fseek(fid, - min(1e4, f_size - 1), 'eof');
                                    buf = fread(fid, min(1e4, f_size - 1),'*char')';
                                    
                                    % Start searching for a valid epoch
                                    time = [];
                                    loop_n = 1;
                                    file_ko = false;
                                    while isempty(time) && ~isempty(buf) && ~file_ko
                                        if ~isempty(find(buf(1:min(1000,numel(buf))) == 13, 1, 'first'))
                                            has_cr = true;  % The file has carriage return - I hate you Bill!
                                        else
                                            has_cr = false;  % The file is UNIX standard
                                        end
                                        
                                        % get new line separators
                                        nl = regexp(buf, '\n')';
                                        if isempty(nl)
                                            file_ko = true;
                                        else                                            
                                            if nl(end) <  (numel(buf) - double(has_cr))
                                                nl = [nl; numel(buf)];
                                            end
                                            lim = [[1; nl(1 : end - 1) + 1] (nl - 1 - double(has_cr))];
                                            lim = [lim lim(:,2) - lim(:,1)];
                                            lim(lim(:,3) < this.id_date(end), :) = [];
                                            
                                            l = size(lim, 1);
                                            while (l >= 1) && isempty(time)
                                                line = buf(lim(l,1) : lim(l,2));
                                                if ~isempty(regexp(line(1:15),'( [0-9]{1,4} [ 0-9]{1}[0-9]{1} [ 0-9]{1}[0-9]{1})', 'once'))
                                                    % Check that the read epoch is within 30 days from the first epoch
                                                    % (if it's further in time it's probably a misleading false epoch line)
                                                    time = GPS_Time(line(this.id_date), [], true);
                                                end
                                                l = l - 1;
                                            end
                                            if ~file_ko
                                                if isempty(time) || time.isEmpty()
                                                    loop_n = loop_n + 1;
                                                    fseek(fid, loop_n * -10000, 'eof'); % If no valid time have been found try to go back more...
                                                    buf = fread(fid, 10000, '*char')';
                                                end
                                            end
                                        end
                                    end
                                    if file_ko
                                        log.addError(sprintf('Check the following file, it seems to be corrupted\nPath: "%s"', full_path));
                                        this.first_epoch = this.first_epoch.getEpoch(1 : this.last_epoch.length);
                                        this.is_valid_list(f) = false;                                        
                                    elseif isempty(time)
                                        log.addError(sprintf('Check the following file, it seems to be empty\nPath: "%s"', full_path));
                                        this.first_epoch = this.first_epoch.getEpoch(1 : this.last_epoch.length);
                                        this.is_valid_list(f) = false;
                                    else
                                        this.last_epoch.append(time);
                                        log.addMessage(sprintf('        last  epoch found at: %s', this.last_epoch.last.toString()), this.verbosity_lev);
                                    end
                                end
                            end
                            
                            if this.is_valid_list(f)
                                if ~isempty(coo)
                                    this.coo.append(Coordinates.fromStringXYZ(coo));
                                else
                                    this.coo.append(Coordinates.fromXYZ([0 0 0]));
                                end
                                this.eoh(f) = eoh;
                            end
                            
                            fclose(fid);
                        end
                    end
                catch ex
                    log.addWarning(['"' this.file_name_list{f} this.ext{f} '" appears to be a corrupted RINEX file or missing'], this.verbosity_lev);
                    %Core_Utils.printEx(ex);
                    this.is_valid_list(f) = false;
                    fclose(fid);
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
                log.addWarning('Some or all the RINEX files are corrupted or missing!!!', this.verbosity_lev);
            end
        end

        function checkValidityUnbuffered(this, flag_header_only)
            % Update the status of validity of the files here pointed
            % Unbuffered read => read line by line (old approach)
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
                
                % try to find the first and the last epoch stored in the file
                try
                    %%
                    fid = fopen(fullfile(this.base_dir{f}, [this.file_name_list{f} this.ext{f}]));
                    if fid < 0
                        this.log.addError(['"' this.file_name_list{f} this.ext{f} '" appears to be missing'], this.verbosity_lev);
                        this.is_valid_list(f) = false;
                    else
                        l = 1;
                        line = fgetl(fid);
                        date_start = '';
                        date_stop = '';
                        coo = '';
                        par_to_find = 4;
                        while par_to_find > 0 && isempty(strfind(line,'END OF HEADER')) && ischar(line) %#ok<*STREMP>
                            l = l + 1;
                            line = fgetl(fid);
                            if numel(line) > 70 && isempty(marker_name) && line(64) == 'K' % read marker_name
                                marker_name = strtrim(regexp(line, '.*(?=MARKER NAME)', 'match', 'once'));
                                par_to_find = par_to_find - 1;
                            end
                            if numel(line) > 76
                                if line(61) == 'T'
                                    %tmp = regexp(line, '.*(?=GPS         TIME OF FIRST OBS)', 'match', 'once');
                                    if line(69) == 'F'
                                        date_start = strtrim(line(1:48));
                                        par_to_find = par_to_find - 1;
                                    elseif line(69) == 'L'
                                        %tmp = regexp(line, '.*(?=GPS         TIME OF LAST OBS)', 'match', 'once');
                                        date_stop = strtrim(line(1:48));
                                        par_to_find = par_to_find - 1;
                                    end
                                else
                                    %tmp = regexp(line, '.*(?=APPROX POSITION XYZ)', 'match', 'once');
                                    if line(62) == 'P' % APPROXIMATE POSITION
                                        coo = strtrim(line(1:60));
                                        par_to_find = par_to_find - 1;
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
                                txt = fread(fid, 10000,'*char')';
                                
                                % Start searching for a valid epoch
                                time = [];
                                loop_n = 1;
                                while isempty(time) && ~isempty(txt)
                                    if ~isempty(find(txt(1:min(1000,numel(txt))) == 13, 1, 'first'))
                                        has_cr = true;  % The file has carriage return - I hate you Bill!
                                    else
                                        has_cr = false;  % The file is UNIX standard
                                    end
                                    
                                    % get new line separators
                                    nl = regexp(txt, '\n')';
                                    if nl(end) <  (numel(txt) - double(has_cr))
                                        nl = [nl; numel(txt)];
                                    end
                                    lim = [[1; nl(1 : end - 1) + 1] (nl - 1 - double(has_cr))];
                                    lim = [lim lim(:,2) - lim(:,1)];
                                    lim(lim(:,3) < 20, :) = [];
                                    
                                    l = size(lim, 1);
                                    while (l >= 1) && isempty(time)
                                        line = txt(lim(l,1) : lim(l,2));
                                        if ~isempty(regexp(line(1:15),'( [0-9]{1,4} [ 0-9]{1}[0-9]{1} [ 0-9]{1}[0-9]{1})', 'once'))
                                            % Check that the read epoch is within 30 days from the first epoch
                                            % (if it's further in time it's probably a misleading false epoch line)
                                            time = GPS_Time(line(this.id_date), [], true);
                                        end
                                        l = l - 1;
                                    end
                                    if isempty(time) || time.isEmpty()
                                        loop_n = loop_n + 1;
                                        fseek(fid, loop_n * -10000, 'eof'); % If no valid time have been found try to go back more...
                                        txt = fread(fid, 10000,'*char')';
                                    end
                                end
                                this.last_epoch.append(time);
                                this.log.addMessage(sprintf('        last  epoch found at: %s', this.last_epoch.last.toString()), this.verbosity_lev);
                            end
                        end
                        fclose(fid);
                        this.is_valid_list(f) = true;
                    end
                catch ex
                    this.log.addWarning(['"' this.file_name_list{f} this.ext{f} '" appears to be a corrupted RINEX file'], this.verbosity_lev);
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
        
        function checkCoordinates(rin_list, force_delete)
            flag_show = false;
            dist_lim = 1; % [m]
            for r = 1 : numel(rin_list)
                xyz = rin_list(r).coo.xyz;
                id_ko = (all(xyz == 0,2));
                xyz(id_ko,:) = [];
                if not(isempty(xyz))
                    xyz = bsxfun(@minus, xyz, median(xyz, 1, 'omitnan'));
                    dist = 1e-3 * (sqrt(xyz(:,1).^2 + xyz(:,2).^2 + xyz(:,3).^2)); % [km]
                    if any(dist > dist_lim)
                        % This set of file RINEX are from a receiver that have been moved
                        % for more than 100 meters
                        
                        marker_name = rin_list(r).marker_name{find(rin_list(r).is_valid_list, 1, 'first')};
                        if flag_show
                            rin_list(r).coo.showPositionENU
                            title(marker_name);
                        end
                        
                        log = Core.getLogger;
                        log.addWarning(sprintf('The RINEX of the station %s are anomalous. The antenna could have been moved or the files are from different omonimous receivers\n > MAX detected distance: %.3g [km]\n\nThe following files will be ignored:', marker_name, max(dist)));
                        id_coo = find(not(id_ko));
                        id_valid = find(rin_list(r).is_valid_list);
                        id_ko = find(dist > dist_lim)';
                        for i = id_ko
                            file_name = fullfile(rin_list(r).base_dir{id_valid(id_coo(i))}, [rin_list(r).file_name_list{id_valid(id_coo(i))}, rin_list(r).ext{id_valid(id_coo(i))}]);
                            if nargin == 2 && force_delete
                                delete(file_name);
                            end
                            log.addMonoMessage(log.indent(sprintf('- %s', file_name)));
                        end
                        rin_list(r).first_epoch.remEpoch(id_coo(id_ko));
                        rin_list(r).last_epoch.remEpoch(id_coo(id_ko));
                        rin_list(r).is_valid_list(id_valid(id_coo(id_ko))) = false;
                        rin_list(r).coo.rem(id_coo(id_ko));
                    end
                end
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
                    if file_number <= numel(this(r).is_valid_list)
                        validity = validity || this(r).is_valid_list(file_number);
                    end
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
       
        function printMissingFiles(rin_list, flag_download)
            % Show all the files that seems to be missing for the
            % processing
            if nargin == 1
                flag_download = false;
            end
            clc
            fprintf('These files are missing or currupted\n')
            for r = 1: numel(rin_list)
                id_ko = find(not(rin_list(r).is_valid_list));
                for i = id_ko
                    file_name = [rin_list(r).file_name_list{i}, rin_list(r).ext{i}];
                    if not(flag_download)
                        fprintf('%s\n', file_name);
                    end
                    if numel(file_name) == 12 % it is a RINEX2 file name
                        % Extracting year and doy
                        year = 2000 + str2double(file_name(end-2:end-1));
                        doy = str2double(file_name(5:7));
                        if flag_download
                            str_cmd = sprintf('/usr/bin/python3 /home/gred/Repositories/goget/goGet.py download -m %s --year-doy %04d %03d -o %s --no-wait', rin_list(r).marker_name{i}, year, doy, [rin_list(r).base_dir{i} filesep]);
                            fprintf('%s\n', str_cmd);
                        end
                    end
                end
                
                
            end
            
        end
        
    end
end
