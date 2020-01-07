% Ini_Manager(fileName, raw_data)
%
% =========================================================================
%   OBJECT Ini_Manager
% =========================================================================
%
% DESCRIPTION:
%   Object to create / read / modify an ini file
%
% EXTENDS: Ini_Reader
%
% INI EXAMPLE:
%   [Section 1]
%       array = [1 2 3]
%       string = "I'm a string"
%       stringCellArray = ["I" "am" "a string cell" "array"]
%   ; I'm a comment
%   # I'm also a comment
%   [Section 2]
%       number = 10
%   [Section n]     # comment
%       double = 10.4
%
% REQUIRES:
%   log:     Logger Class
%   cprintf:    http://www.mathworks.com/matlabcentral/fileexchange/24093-cprintf-display-formatted-colored-text-in-the-command-window
%

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
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

classdef Ini_Manager < handle

    properties (Constant, GetAccess = 'private')
        STD_COMMENT  = '#';       % Character Identifying the start of a comments
    end

    properties (GetAccess = 'private', SetAccess = 'private')
        log = Logger.getInstance();        % Handler to the log object
        raw_data = {};                 % Cell array containing the file, each cell is a line
    end

    properties (GetAccess = 'public', SetAccess = 'protected')
        c_comment  = [';', Ini_Manager.STD_COMMENT];  % Character Identifying the start of a comments
        file_name = 'config.ini';         % Name (and full path) of the ini
        fid = 0;                          % Handle of the ini file
        rw = 'r';                         % File access mode (r/w)
        read_status = false;              % Flag of reading status

        % Structure containing the data: section -> name / ( key -> (name / data )
        section = {};
        %section{1}.name = 'Section Name';
        %section{1}.key = {};
        %section{1}.key{1}.name = 'Key name';
        %section{1}.key{1}.data = 100;
    end

    % =========================================================================
    %    STANDARD METHODS
    % =========================================================================
    methods

        % Creator ---------------------------------------------------------
        function this = Ini_Manager(file_name, raw_data)
        % Creator of the class
        % SYNTAX:  Ini_Manager(<file_name>, <raw_data>) creator
            if (nargin == 0)
                file_name = '';
            end
            if isempty(file_name)
                file_name = '';
            end

            if (nargin == 1) && iscell(file_name) % ---------------- Parse RAW but do not save it
                this.setFileName('');
                this.raw_data = file_name;
                this.setRW('w');
                this.cleanRaw();
                this.parseData();
                this.setReadStatus(true)
                this.raw_data = {};
            elseif (nargin == 2) % --------------------------------- Parse and save RAW
                this.setFileName(file_name);
                this.raw_data = raw_data;
                this.setRW('w');
                this.writeFile();
                this.cleanRaw();
                this.parseData();
                this.setReadStatus(true)
                this.raw_data = {};
            else
                this.setFileName(file_name);
            end
        end

        % Destructor ------------------------------------------------------
        function delete(this) %#ok<INUSD>
            % Empty Destructor
        end

        % Cleaner ---------------------------------------------------------
        function clearAll(this)
            % Reset the object
            this.fid = 0;
            this.rw = 'r';
            this.file_name = 'config.ini';
            this.raw_data = {};
            this.section = {};
        end

    end

    % =========================================================================
    %  FILE
    % =========================================================================
    methods

        % Set file name ---------------------------------------------------
        function setFileName(this, file_name)
            % Set the complete path of the ini file
            this.file_name = file_name;
        end

        % Get file name ---------------------------------------------------
        function file_name = getFileName(this)
            % Get the complete path of the ini file
            file_name = this.file_name;
        end

    end

    % =========================================================================
    %  INI SYNTAX
    % =========================================================================
    methods

        % Set Comment Character -------------------------------------------
        function setCommentChar(this, character)
            % define the characters that define comment lines
            this.c_comment = unique([character this.STD_COMMENT]);
        end

        % Get Comment Character -------------------------------------------
        function characters = getCommentChar(this)
            % Get the characters that define comment lines
            characters = this.c_comment;
        end

    end

    % =========================================================================
    %  R/W MODES OF THE FILE
    % =========================================================================
    methods

        % Set RW mode -----------------------------------------------------
        function setRW(this,rwMode)
            % Set mode for open the ini file
            if ( ~this.readMode() && ~writeMode())
                this.rw = 'r';
            else
                this.rw = rwMode;
            end
        end

        % Get RW mode -----------------------------------------------------
        function rwMode = getRW(this)
            % Get the current mode of the file
            rwMode = this.rw;
        end

        % Read mode? ------------------------------------------------------
        function bool = readMode(this)
            % Is it in read mode?
            if ((this.rw(1)=='r'))
                bool = true;
            else
                bool = false;
            end
        end

        % Write mode? -----------------------------------------------------
        function bool = writeMode(this)
            % Is it in write mode?
            if ((this.rw(1)=='w'))
                bool = true;
            else
                bool = false;
            end
        end

        % Get read status mode --------------------------------------------
        function bool = getReadStatus(this)
            % Return if the file has been already parsed
            bool = this.read_status;
        end

        % Read File -------------------------------------------------------
        function status = readFile(this)
            % Force reading of the File
            status = true;
            % If the object already contains data - clean it
            if (this.getReadStatus())
                this.raw_data = {};
            end
            if ~exist(this.getFileName(),'file')
                this.log.addError(sprintf('INI reading of file "%s" failed', this.getFileName()))
                status = false;
            else
                this.fid = fopen(this.getFileName(), this.getRW());

                if (this.fid ~= -1)   % If reading is ok
                    this.raw_data = textscan(this.fid, '%s', 'delimiter', '\n', 'endOfLine', '\r\n');
                    fclose(this.fid);

                    this.read_status = true;
                    this.raw_data = this.raw_data{1};

                    this.log.addStatusOk('The INI file has been parsed correctly.', 15);

                    this.cleanRaw();      % Strip comments and spaces
                    this.parseData();     % Parse file
                    this.raw_data = {};   % clean RAW temp data
                else
                    this.log.addError(['INI file read failed:' ferror(this.fid)]);
                    this.setReadStatus(false);
                    this.raw_data = {};
                end
            end
        end

        % Write File ------------------------------------------------------
        function errStatus = writeFile(this)
            % Force reading of the File
            errStatus = false;
            if isempty(this.getFileName())
                this.log.addError(sprintf('INI writing  of file "%s" failed', this.getFileName()));
                errStatus = true;
            else
                if (this.fid ~= -1)   % If file access is ok
                    try
                        this.fid = fopen(this.getFileName(), this.getRW());
                        % Convert raw data to string
                        tmp_str = '';
                        for i = 1 : numel(this.raw_data)
                            tmp_str = [tmp_str this.raw_data{i} 10]; %#ok<AGROW>
                        end
                        fwrite(this.fid, tmp_str);
                        fclose(this.fid);
                    catch ex
                        errStatus = true;
                        this.log.addError(['INI file cannot be written (' this.file_name '): ' ex.message]);
                    end

                    this.log.addStatusOk('The INI file has been writted correctly', 10);
                else
                    this.log.addError(['INI file write failed:' ferror(this.fid)]);
                end
            end
        end

        % Update File (when needed) ---------------------------------------
        function reloaded = update(this, file_name, force_read)
            % Update the object when needed:
            %  - file_name changed
            %  - force flag == 1
            %  - INI not yet read
            if nargin == 2
                force_read = 0;
            end
            reloaded = 0;
            if (~strcmp(file_name,this.getFileName) || (force_read == 1) || ~this.getReadStatus())
                this.setFileName(file_name);
                this.readFile();
                reloaded = 1;
            end
        end

    end

    % =========================================================================
    %  SEARCH FUNCTIONS
    % =========================================================================
    methods

        % Search a section in the object Ini_Manager ----------------------
        function isPresent = containsSection(this, section)
            % Search a section in the object Ini_Manager
            s = 1;
            while ((s<=length(this.section)) && (s ~= 0))
                if (strcmp(this.section{s}.name,section))
                    s = 0;
                else
                    s = s + 1;    % go on with the search of the section
                end
            end
            if (s == 0)
                isPresent = true;
            else
                isPresent = false;
            end
        end

        % Search a key in the object Ini_Manager --------------------------
        function isPresent = containsKey(this, section, key)
            % Search a key in the object Ini_Manager
            s = 1;
            while ((s<=length(this.section)) && (s ~= 0))
                if (strcmp(this.section{s}.name,section))
                    while ((k<=length(this.section{s}.key)) && (k ~= 0))
                        if (strcmp(this.section{s}.key{k}.name,key))
                            k = 0;
                        else
                            k = k + 1;
                        end
                    end
                    s = 0;
                else
                    s = s + 1;    % go on with the search of the section
                end
            end
            if (k == 0)
                isPresent = true;
            else
                isPresent = false;
            end
        end

    end

    % =========================================================================
    %  GETTER OF THE PARAMETERS
    % =========================================================================
    methods

        % Get the list of sections present in the file --------------------
        function section_list = getSections(this)
            % Get the list of available sections
            if (~this.getReadStatus())
                %this.printWarning('File not yet read!\n');
                this.readFile();
            end

            section_list = cell(length(this.section),1);
            for s=1:length(this.section)
                section_list{s} = this.section{s}.name;
            end
        end

        % Return the presence of a section --------------------------------
        function isS = isSection(this, section)
            % Get the presence of a section
            s = 1;
            while ((s<=length(this.section)) && (s ~= 0))
                if (strcmp(this.section{s}.name,section))
                    s = 0;      % Stop searching
                else
                    s = s+1;    % go on with the search of the section
                end
            end
           isS = s == 0;
        end

        % Return the presence of a key ------------------------------------
        function isK = isKey(this, section, key)
            % Get the presence of a key
            s = 1;
            k = 1;
            while ((s<=length(this.section)) && (s ~= 0))
                if (strcmp(this.section{s}.name,section))
                    k = 1;
                    while ((k<=length(this.section{s}.key)) && (k ~= 0))
                        if (strcmp(this.section{s}.key{k}.name,key))
                            k = 0;
                        else
                            k = k + 1;
                        end
                    end
                    s = 0;
                else
                    s = s + 1;    % go on with the search of the section
                end
            end
            isK = k == 0;
        end

        % Get the list of keys present in the file ------------------------
        function key_list = getKeys(this, section)
            % Get the list of the keys available
            if (nargin == 1)
                section = 0;
            end

            if (~this.getReadStatus())
                %this.printWarning('File not yet read!\n');
                this.readFile();
            end

            key_list = {};
            l = 1;
            if (section == 0)
                for s = 1:length(this.section)
                    for p = 1 : length(this.section{s}.key)
                        key_list{l} = this.section{s}.key{p}.name; %#ok<AGROW>
                        l = l+1;
                    end
                end
            else
                s = 1;
                while ((s<=length(this.section)) && (s ~= 0))
                    if (strcmp(this.section{s}.name,section))
                        for p = 1:length(this.section{s}.key)
                            key_list{l} = this.section{s}.key{p}.name; %#ok<AGROW>
                            l = l+1;
                        end
                        s = 0;      % Stop searching
                    else
                        s = s+1;    % go on with the search of the section
                    end
                end
            end
        end

        % Get data --------------------------------------------------------
        function data = getData(this, section, key)
            % Get the value of a specified key
            if (~this.getReadStatus() && isempty(this.section))
               %this.printWarning('File not yet read!\n');
               this.readFile();
            end
            found = false;
            data = [];
            if (nargin == 2)
                key = section;
                % Search the key among all the sections
                s = 1;
                while ((s<=length(this.section)) && (s > 0))
                    p = 1;
                    while ((p <= length(this.section{s}.key)) && (p ~= 0))
                        if (strcmp(this.section{s}.key{p}.name,key))
                            data = this.section{s}.key{p}.data;
                            found = true;
                            p = 0;      % Stop searching key
                        else
                            p = p+1;    % go on with the search of the key
                        end
                    end
                    if (p == 0)
                        s = 0;     % Stop searching section
                    else
                        s = s+1;   % go on with the search of the section
                    end
                end
                if ~found
                    this.log.addWarning(['Key "' key '" not found while reading: "' this.file_name '"'], 10);
                    data = [];
                end
            else
                % Search the key of a specific section
                s = 1;
                while ((s<=length(this.section)) && (s > 0))
                    if (strcmp(this.section{s}.name,section))
                        p = 1;
                        while ((p <= length(this.section{s}.key)) && (p ~= 0))
                            if (strcmp(this.section{s}.key{p}.name,key))
                                data = this.section{s}.key{p}.data;
                                found = true;
                                p = 0; % Stop searching key
                            else
                                p = p+1;    % go on with the search of the key
                            end
                        end
                        s = 0;      % Stop searching section
                    else
                        s = s+1;    % go on with the search of the section
                    end
                end
                if ~found
                    this.log.addWarning(['Key "' key '" not found in section "' section '" while reading: "' this.file_name '"'], 100);
                    data = [];
                end
            end

        end

    end

    % =========================================================================
    %  MODIFIER FUNCTIONS
    % =========================================================================
    methods

        % Add a new section to the object ---------------------------------
        function addSection(this, new_section)
            % Add a new section to the object
            if ~this.isSection(new_section)
                newId = length(this.section)+1;
                this.section{newId}.name = new_section;
                this.section{newId}.key = {};
            end
        end

        % Add a new keys to the object ------------------------------------
        function addKey(this, section, key, data)
            % Add a new key to the object if the section exist
            if this.isKey(section, key)
                this.editKey(section, key, data);
            else
                s = 1;
                while ((s<=length(this.section)) && (s ~= 0))
                    if (strcmp(this.section{s}.name,section))
                        newId = length(this.section{s}.key)+1;
                        this.section{s}.key{newId}.name = key;
                        this.section{s}.key{newId}.data = data;
                        s = 0;
                    else
                        s = s+1;    % go on with the search of the section
                    end
                end
                if (s ~= 0)
                    this.printError(['Section "' section '" not found!\n']);
                end
            end
        end

        % Add a data to the object giving key and section -----------------
        function addData(this, section, key, data)
            % Add a data to the object an eventually create the section and key to store it
            this.addSection(section)
            this.addKey(section, key, data)
        end

        % Remove a section from the object Ini_Manager --------------------
        function rmSection(this, section)
            % Remove a section from the object Ini_Manager
            s = 1;
            while ((s<=length(this.section)) && (s ~= 0))
                if (strcmp(this.section{s}.name,section))
                    this.section(s) = [];
                    s = 0;
                else
                    s = s+1;    % go on with the search of the section
                end
            end
            if (s ~= 0)
                this.printError(['Section "' section '" not found!\n']);
            end
        end

        % Remove a key from the object Ini_Manager ------------------------
        function rmKey(this, section, key)
            % Remove a key from the object Ini_Manager
            s = 1;
            k = 1;
            while ((s<=length(this.section)) && (s ~= 0))
                if (strcmp(this.section{s}.name,section))
                    k = 1;
                    while ((k<=length(this.section{s}.key)) && (k ~= 0))
                        if (strcmp(this.section{s}.key{k}.name,key))
                            this.section{s}.key(k) = [];
                            k = 0;
                        else
                            k = k + 1;
                        end
                    end
                    s = 0;
                else
                    s = s + 1;    % go on with the search of the section
                end
            end
            if (k ~= 0)
                this.printError(['Key "' key '" not found!\n']);
            end
        end

        % Edit a key in the object Ini_Manager ----------------------------
        function editKey(this, section, key, data)
            this.setData(section, key, data)
        end

        function setData(this, section, key, data)
            % Edit a key in the object Ini_Manager
            s = 1;
            k = 1;
            if (nargin == 3)
                data  = key;
                key = section;
                % Search the key among all the sections
                while ((s<=length(this.section)) && (s > 0))
                    p = 1;
                    while ((p <= length(this.section{s}.key)) && (p ~= 0))
                        if (strcmp(this.section{s}.key{p}.name,key))
                            this.section{s}.key{p}.data = data;
                            p = 0;      % Stop searching key
                        else
                            p = p+1;    % go on with the search of the key
                        end
                    end
                    if (p == 0)
                        s = 0;     % Stop searching section
                    else
                        s = s+1;   % go on with the search of the section
                    end
                end
                if (isempty(data))
                    this.log.addWarning(['Key "' key '" not found while reading: "' this.file_name '"'], 10);
                end
            else
                while ((s<=length(this.section)) && (s ~= 0))
                    if (strcmp(this.section{s}.name,section))
                        k = 1;
                        while ((k<=length(this.section{s}.key)) && (k ~= 0))
                            if (strcmp(this.section{s}.key{k}.name,key))
                                this.section{s}.key{k}.data = data;
                                k = 0;
                            else
                                k = k + 1;
                            end
                        end
                        s = 0;
                    else
                        s = s + 1;    % go on with the search of the section
                    end
                end
                if (k ~= 0)
                    this.log.addWarning(['Key "' key '" not found while reading: "' this.file_name '"'], 10);
                end
            end
        end

    end

    % =========================================================================
    %  VISUALIZATION FUNCTIONS (not using LOGGER) - for DEBUG
    % =========================================================================
    methods

        % List Sections ---------------------------------------------------
        function listSections(this, color_mode)
            % List the list of available sections
            if (nargin == 1)
                color_mode = this.log.getColorMode();
            end

            if (~this.getReadStatus())
                %this.printWarning('File not yet read!\n');

                this.readFile();
            end

            fprintf('        List of sections:\n')
            for s=1:length(this.section)
                this.printSection(this.section{s}.name, color_mode);
            end
        end

        % List Key --------------------------------------------------------
        function listKeys(this, section, color_mode)
            % List the list of the keys available
            if (nargin == 1)
                section = 0;
                color_mode = this.log.getColorMode();
            end
            if (nargin == 2)
                if (ischar(section))
                    color_mode = this.log.getColorMode();
                else
                    color_mode = section;
                    section = 0;
                end
            end

            if (~this.getReadStatus() && isempty(this.section))
                %this.printWarning('File not yet read!\n');
                this.readFile();
            end

            if (section == 0)
                fprintf('        List of all the keys:\n')
                for s = 1:length(this.section)
                    this.printSection(this.section{s}.name, color_mode);
                    for p = 1:length(this.section{s}.key)
                        this.printData(this.section{s}.key{p}.name, this.section{s}.key{p}.data, color_mode)
                    end
                end
            else
                fprintf('        List of all the keys in section "%s":\n', section)
                s = 1;
                while ((s<=length(this.section)) && (s ~= 0))
                    if (strcmp(this.section{s}.name,section))
                        this.printSection(this.section{s}.name, color_mode);
                        for p = 1:length(this.section{s}.key)
                            this.printData(this.section{s}.key{p}.name, this.section{s}.key{p}.data, color_mode)
                        end
                        s = 0;      % Stop searching
                    else
                        s = s+1;    % go on with the search of the section
                    end
                end
            end
        end

        % Show data -------------------------------------------------------
        function showData(this, section, color_mode)
            % List the data contained in the ini
            if (nargin == 1)
                section = 0;
                color_mode = this.log.getColorMode();
            end
            if (nargin == 2)
                if (ischar(section))
                    color_mode = this.log.getColorMode();
                else
                    color_mode = section;
                    section = 0;
                end
            end

            if (~this.getReadStatus())
                %this.printWarning('File not yet read')
                this.readFile();
            end

            if (section == 0)
                fprintf('        List of all the key values:\n')
                for s = 1:length(this.section)
                    this.printSection(this.section{s}.name, color_mode);
                    for p = 1:length(this.section{s}.key)
                        this.printData(this.section{s}.key{p}.name, this.section{s}.key{p}.data, color_mode)
                    end
                end
            else
                fprintf('        List of all the key values in section "%s":\n', section)
                s = 1;
                while ((s <= length(this.section)) && (s ~= 0))
                    if (strcmp(this.section{s}.name,section))
                        this.printSection(this.section{s}.name, color_mode);
                        for p = 1:length(this.section{s}.key)
                            this.printData(this.section{s}.key{p}.name, this.section{s}.key{p}.data, color_mode)
                        end
                        s = 0;      % Stop searching
                    else
                        s = s+1;    % go on with the search of the section
                    end
                end
            end
        end

    end


    % =========================================================================
    %   CUSTOM DISPLAY UTILITIES (not using LOGGER) - for DEBUG
    % =========================================================================

    methods (Access = 'protected')

        % Display a Section -----------------------------------------------
        function printSection(this, section, color_mode)
            % Display a sections on a tree
            if (nargin == 2)
                color_mode = this.color_mode;
            end
            if (color_mode)
                cprintf('blue','         - ');
                cprintf('text', '"%s"\n', section);
            else
                fprintf('         - "%s"\n', section);
            end
        end

        % Display a Key ---------------------------------------------------
        function printKey(this, key, color_mode)
            % Display a key on a tree
            if (nargin == 2)
                color_mode = this.color_mode;
            end
            if (color_mode)
                cprintf('blue','             |- ');
                cprintf('text', '"%s"\n', key);
            else
                fprintf('             |- "%s"\n', key);
            end
        end

        % Display a Key ---------------------------------------------------
        function printData(this, key, data, color_mode)
            % Display the data of a key
            if (nargin == 3)
                color_mode = this.color_mode;
            end
            if (color_mode)
                cprintf('blue','             |- ');
                tmpData = data;
                if (iscell(tmpData))    % The only dataset that can be a cell array is an array of strings
                    cprintf('text','"%s" = [', key);
                    for l = 1:size(tmpData,2)
                        cprintf('text', ' "%s"', tmpData{l});
                    end
                    cprintf('text', ' ]\n');
                else
                    % if it is not an array of string...
                    if (ischar(tmpData))
                        cprintf('text','"%s" = "%s"\n',key, tmpData);
                    else
                        cprintf('text','"%s" = %s\n',key, num2str(tmpData));
                    end
                end
            else
                tmpData = data;
                if (iscell(tmpData))    % The only dataset that can be a cell array is an array of strings
                    fprintf('             |- "%s" = [', key);
                    for l = 1:size(tmpData,2)
                        fprintf(' "%s"', tmpData{l});
                    end
                    fprintf(' ]\n');
                else
                    % if it iks not an array of string...
                    if (ischar(tmpData))
                        fprintf('             |- "%s" = "%s"\n',key, num2str(tmpData));
                    else
                        fprintf('             |- "%s" = %s\n',key, num2str(tmpData));
                    end
                end
            end
        end

    end

    % =========================================================================
    %    PRIVATE UTILITIES
    % =========================================================================
    methods (Access = 'private')

        % Strip empty/comment lines ---------------------------------------
        function cleanRaw(this)
            % Strip empty lines
            this.raw_data((cellfun(@length, this.raw_data)==0)) = [];
            % Strip commented lines
            r = 1;
            while (r <= length(this.raw_data))
                if ismember(this.raw_data{r}(1), this.c_comment)
                    this.raw_data(r) = [];
                else
                    r = r + 1;
                end
            end

            for r=1:length(this.raw_data)
                % Strip inline comments

                start = regexp(this.raw_data{r}, ['([' this.c_comment ']).*'], 'start');
                if (~isempty(start))
                    this.raw_data{r}(start:end) = [];
                end

                start = regexp(this.raw_data{r}, ' *$', 'start');
                if (~isempty(start))
                    this.raw_data{r}(start:end) = [];
                end
            end
        end

        % Parse the file (cleaned) ----------------------------------------
        function parseData(this)
            p = 0;
            s = 0;
            this.section = {};
            for r=1:length(this.raw_data)
                sectionName = regexp(this.raw_data{r}, '(?<=^\[).*(?=\])', 'match');
                if (~isempty(sectionName))
                    % we have a new SECTION!
                    s = s+1;
                    p = 0;
                    this.section{s}.name = sectionName{1};
                    this.section{s}.key = [];
                else
                    % maybe this line contains a key, let's check
                    parName = regexp(this.raw_data{r}, '^(.*(?=\s\=))', 'match');
                    if (isempty(parName))
                        parName = regexp(this.raw_data{r}, '^(.*(?=\=))', 'match');
                    end
                    if (~isempty(parName))
                        % we have a PARAMETER!
                        p = p+1;
                        if (s == 0)
                            s = s + 1;
                            this.section{s}.name = 'Section';
                            this.section{s}.key = [];
                        end
                        this.section{s}.key{p}.name = parName{1};
                        % Get the DATA!
                        strData = regexp(this.raw_data{r}, '(?<=\=).*', 'match');
                        if (isempty(strData))
                            this.section{s}.key{p}.data = [];
                        else
                            tmpData = str2num(strData{1}); %#ok<ST2NM>
                            if (~isempty(tmpData) && isnumeric(tmpData))   % It's a number!
                                this.section{s}.key{p}.data = tmpData;
                            else
                                % It's a string!
                                % If it is a string properly formatted "string data"
                                tmpData = regexp(strData{1}, '(?<=\")[^"]*(?=\")', 'match');
                                if (~isempty(tmpData))
                                    if (size(tmpData,2) == 1)
                                        this.section{s}.key{p}.data = tmpData{1};
                                    else
                                        % Strip empty cells
                                        tmpData((cellfun(@length, tmpData)==0)) = [];
                                        % Strip space only cells
                                        for c=size(tmpData,2):-1:1
                                            if isempty(regexp(tmpData{c}, '[^ ]*', 'match'))
                                                tmpData(c) = [];
                                            end
                                        end
                                        this.section{s}.key{p}.data = tmpData;
                                    end
                                else
                                    % If it appears to be an empty string
                                    this.section{s}.key{p}.data = '';
                                end
                            end
                        end
                    end

                end

            end
        end

        % Set read status -------------------------------------------------
        function setReadStatus(this, bool)
            this.read_status = bool;
        end
    end

    % =========================================================================
    %  STATIC STRING GENERATION FROM INI TO STR
    % =========================================================================

    methods (Static)

        % cellStringtoString ----------------------------------------------
        function str = strCell2Str(value)
            % Converta cell of string to string
            % SYNTAX:
            %   str = strCell2Str(value)
            if ~isempty(value)
                if ~iscell(value)
                    value = {value};
                end
                str = strcat('"', value{1}, '"');
                for i = 2 : numel(value)
                    str = strcat(str, ' "', value{i}, '"');
                end
            else
                str = '[]';
            end
        end

        % toIniString -----------------------------------------------------
        function cell_str = toIniString(variable_name, value, format, cell_str)
            % Convert any variable to ini string format
            % SYNTAX:
            %   cell_str = toIniString(variable_name, value)
            %   cell_str = toIniString(variable_name, value, format)
            %   cell_str = toIniString(variable_name, value, cell_str)
            %   cell_str = toIniString(variable_name, value, format, cell_str)
            switch nargin
                case 1
                    error('Error in Ini_Manager.toIniString, too few parameters');
                case 2

                    format = '';
                    cell_str = {};
                case 3
                    if iscellstr(format)
                        cell_str = format;
                        format = '%g';
                    else
                        cell_str = {};
                    end
                case 4
                    % ok
            end

            if ischar(value) % is string
                cell_str{numel(cell_str) + 1} = [variable_name ' = "' value(:)' '"'];
            elseif isnumeric(value)
                if isempty(format)
                    format = '%g';
                else
                    format = strtrim(format);
                end
                if numel(value) > 1 % it's an array of values
                    format = [format ' '];
                    tmp = sprintf(format,value);
                    cell_str{numel(cell_str) + 1} = [variable_name ' = [' tmp(1:end-1) ']'];
                    clear tmp;
                else
                    cell_str{numel(cell_str) + 1} = [variable_name ' = ' sprintf(format,value)];
                end
            else % generic converter (may not work properly)
                toString = @(var) strtrim(regexprep(evalc('disp(var)'), '\n', ''));
                if iscell(value)
                    if ~isempty(value) && ischar(value{1})
                        cell_str{numel(cell_str) + 1} = [variable_name ' = [' Ini_Manager.strCell2Str(value) ']'];
                    else
                        cell_str{numel(cell_str) + 1} = [variable_name ' = [' toString(value) ']'];
                    end
                else
                    cell_str{numel(cell_str) + 1} = [variable_name ' = ' toString(value)];
                end
            end

            % I want a column array
            if size(cell_str,1) < size(cell_str,2)
                cell_str = cell_str';
            end
        end

        % toIniString Section ---------------------------------------------
        function cell_str = toIniStringSection(section_name, cell_str)
            % Add a section in ini string format
            % SYNTAX:
            %   cell_str = toIniStringSection(section_name, cell_str)
            if (nargin == 1)
                cell_str = {};
            end

            cell_str{numel(cell_str) + 1} = [Ini_Manager.STD_COMMENT '----------------------------------------------------------------------------------'];
            cell_str{numel(cell_str) + 1} = ['[' section_name ']'];
            cell_str{numel(cell_str) + 1} = [Ini_Manager.STD_COMMENT '----------------------------------------------------------------------------------'];

            % I want a column array
            if size(cell_str,1) < size(cell_str,2)
                cell_str = cell_str';
            end
        end

        % toIniString Comment ---------------------------------------------
        function cell_str = toIniStringComment(comment, cell_str)
            % Add a comment in ini string format
            % SYNTAX:
            %   cell_str = toIniStringComment(comment, cell_str)
            if (nargin == 1)
                cell_str = {};
            end
            line_lim = strfind(comment, 10);
            line_lim = [[1; line_lim'] [line_lim'; numel(comment)]];
            for i = 1 : size(line_lim, 1)            
                cell_str{numel(cell_str) + 1} = [ Ini_Manager.STD_COMMENT ' ' strrep(comment(line_lim(i, 1) : line_lim(i, 2)), char(10), '')];
            end

            % I want a column array
            if size(cell_str,1) < size(cell_str,2)
                cell_str = cell_str';
            end
        end

        % toIniString New Line --------------------------------------------
        function cell_str = toIniStringNewLine(cell_str)
            % Add a new line in ini string format
            % SYNTAX:
            %   cell_str = toIniStringSection(section_name, cell_str)
            if (nargin == 0)
                cell_str = {};
            end
            cell_str{numel(cell_str) + 1} = '';

            % I want a column array
            if size(cell_str,1) < size(cell_str,2)
                cell_str = cell_str';
            end
        end
    end
end
