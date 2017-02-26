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
%   logger:     Logger Class
%   cprintf:    http://www.mathworks.com/matlabcentral/fileexchange/24093-cprintf-display-formatted-colored-text-in-the-command-window
%

%--------------------------------------------------------------------------
%               ___ ___ ___ 
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.9.1
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

classdef Ini_Manager < handle
      
    properties (Constant, GetAccess = 'private')
        STD_COMMENT  = '#';       % Character Identifying the start of a comments
    end
      
    properties (GetAccess = 'private', SetAccess = 'private')
        logger = Logger.getInstance(); % Handler to the logger object
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
        function obj = Ini_Manager(file_name, raw_data)
        % Ini_Manager(fileName, verbosity) creator
            if (nargin == 0)
                file_name = '';
            end
            if isempty(file_name)
                file_name = '';
            end
            obj.setFileName(file_name);
            if nargin == 2
                obj.raw_data = raw_data;
                obj.setRW('w');
                obj.writeFile();
                obj.cleanRaw();
                obj.parseData();
                obj.setReadStatus(true)
                obj.raw_data = {};
            end                
        end
        
        % Distructor ------------------------------------------------------
        function delete(obj)
            % Empty Destructor
        end
        
        % Cleaner ---------------------------------------------------------
        function clearAll(obj)
            % Reset the object
            obj.fid = 0;
            obj.rw = 'r';
            obj.file_name = 'config.ini';
            obj.raw_data = {};
            obj.section = {};
        end
        
    end
    
    % =========================================================================
    %  FILE
    % =========================================================================
    methods
        
        % Set file name ---------------------------------------------------
        function setFileName(obj, file_name)
            % Set the complete path of the ini file
            obj.file_name = file_name;
        end
        
        % Get file name ---------------------------------------------------
        function file_name = getFileName(obj)
            % Get the complete path of the ini file
            file_name = obj.file_name;
        end
        
    end

    % =========================================================================
    %  INI SYNTAX
    % =========================================================================
    methods
        
        % Set Comment Character -------------------------------------------
        function setCommentChar(obj, character)
            % define the characters that define comment lines
            obj.c_comment = unique([character obj.STD_COMMENT]);
        end
        
        % Get Comment Character -------------------------------------------
        function characters = getCommentChar(obj)
            % Get the characters that define comment lines
            characters = obj.c_comment;
        end

    end

    % =========================================================================
    %  R/W MODES OF THE FILE
    % =========================================================================
    methods
                
        % Set RW mode -----------------------------------------------------
        function setRW(obj,rwMode)
            % Set mode for open the ini file
            if ( ~obj.readMode() && ~writeMode())
                obj.rw = 'r';
            else
                obj.rw = rwMode;
            end
        end
        
        % Get RW mode -----------------------------------------------------
        function rwMode = getRW(obj)
            % Get the current mode of the file
            rwMode = obj.rw;
        end
        
        % Read mode? ------------------------------------------------------
        function bool = readMode(obj)
            % Is it in read mode?
            if ((obj.rw(1)=='r'))
                bool = true;
            else
                bool = false;
            end
        end
        
        % Write mode? -----------------------------------------------------
        function bool = writeMode(obj)
            % Is it in write mode?
            if ((obj.rw(1)=='w'))
                bool = true;
            else
                bool = false;
            end
        end
        
        % Get read status mode --------------------------------------------        
        function bool = getReadStatus(obj)
            % Return if the file has been already parsed
            bool = obj.read_status;
        end
        
        % Read File -------------------------------------------------------
        function errStatus = readFile(obj)
            % Force reading of the File
            errStatus = false;            
            % If the object already contains data - clean it
            if (obj.getReadStatus())
                obj.raw_data = {};
            end
            if ~exist(obj.getFileName(),'file')
                obj.logger.addError('INI file name not set');                
                errStatus = true;
            else
                obj.fid = fopen(obj.getFileName(), obj.getRW());
                
                if (obj.fid ~= -1)   % If reading is ok
                    obj.raw_data = textscan(obj.fid, '%s', 'delimiter', '\n', 'endOfLine', '\r\n');                    
                    fclose(obj.fid);
                    
                    obj.read_status = true;
                    obj.raw_data = obj.raw_data{1};
                    
                    obj.logger.addStatusOk('The INI file has been parsed correctly.', 15);
                    
                    obj.cleanRaw();      % Strip comments and spaces
                    obj.parseData();     % Parse file
                    obj.raw_data = {};   % clean RAW temp data
                else
                    obj.logger.addError(['INI file read failed:' ferror(obj.fid)]);
                    obj.setReadStatus(false);
                    obj.raw_data = {};
                end
            end
        end
       
        % Write File -------------------------------------------------------
        function errStatus = writeFile(obj)
            % Force reading of the File
            errStatus = false;            
            if isempty(obj.getFileName())
                obj.logger.addError('INI file name not set');
                errStatus = true;
            else                
                if (obj.fid ~= -1)   % If file access is ok                    
                    try
                        obj.fid = fopen(obj.getFileName(), obj.getRW());
                        %fprintf(obj.fid, char(obj.raw_data));
                        fwrite(obj.fid, sprintf('%s\n', string(obj.raw_data)));
                        fclose(obj.fid);
                    catch ex
                        errStatus = true;
                        obj.logger.addError(['INI file cannot be written (' obj.file_name '): ' ex.message]);
                    end
                    
                    obj.logger.addStatusOk('The INI file has been writted correctly', 10);
                else
                    obj.logger.addError(['INI file write failed:' ferror(obj.fid)]);
                end
            end
        end
        
        % Update File (when needed) ---------------------------------------
        function reloaded = update(obj, file_name, force_read)
            % Update the object when needed:
            %  - file_name changed 
            %  - force flag == 1
            %  - INI not yet read
            if nargin == 2
                force_read = 0;
            end
            reloaded = 0;
            if (~strcmp(file_name,obj.getFileName) || (force_read == 1) || ~obj.getReadStatus())
                obj.setFileName(file_name);
                obj.readFile();
                reloaded = 1;
            end            
        end
        
    end
    
    % =========================================================================
    %  SEARCH FUNCTIONS
    % =========================================================================
    methods
        
        % Search a section in the object Ini_Manager ----------------------
        function isPresent = containsSection(obj, section)
            % Search a section in the object Ini_Manager
            s = 1;
            while ((s<=length(obj.section)) && (s ~= 0))
                if (strcmp(obj.section{s}.name,section))
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
        function isPresent = containsKey(obj, section, key)
            % Search a key in the object Ini_Manager
            s = 1;
            while ((s<=length(obj.section)) && (s ~= 0))
                if (strcmp(obj.section{s}.name,section))
                    while ((k<=length(obj.section{s}.key)) && (k ~= 0))
                        if (strcmp(obj.section{s}.key{k}.name,key))
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
        function section_list = getSections(obj)
            % Get the list of available sections
            if (~obj.getReadStatus())
                %obj.printWarning('File not yet read!\n');
                obj.readFile();
            end
            
            section_list = cell(length(obj.section),1);
            for s=1:length(obj.section)
                section_list{s} = obj.section{s}.name;
            end
        end 
        
        % Return the presence of a section --------------------------------
        function isS = isSection(obj, section)
            % Get the presence of a section
            s = 1;
            while ((s<=length(obj.section)) && (s ~= 0))
                if (strcmp(obj.section{s}.name,section))
                    s = 0;      % Stop searching
                else
                    s = s+1;    % go on with the search of the section
                end
            end
           isS = s == 0;
        end
        
        % Return the presence of a key ------------------------------------
        function isK = isKey(obj, section, key)
            % Get the presence of a key 
            s = 1;
            k = 1;
            while ((s<=length(obj.section)) && (s ~= 0))
                if (strcmp(obj.section{s}.name,section))
                    k = 1;
                    while ((k<=length(obj.section{s}.key)) && (k ~= 0))
                        if (strcmp(obj.section{s}.key{k}.name,key))
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
        function keyList = getKeys(obj, section)
            % Get the list of the keys available
            if (nargin == 1)
                section = 0;
            end
            
            if (~obj.getReadStatus())                
                %obj.printWarning('File not yet read!\n');
                obj.readFile();
            end
            
            keyList = {};
            l = 1;
            if (section == 0)
                for s = 1:length(obj.section)
                    for p = 1:length(obj.section{s}.key)
                        keyList{l} = obj.section{s}.key{p}.name;
                        l = l+1;
                    end
                end
            else
                s = 1;
                while ((s<=length(obj.section)) && (s ~= 0))
                    if (strcmp(obj.section{s}.name,section))
                        for p = 1:length(obj.section{s}.key)
                            keyList{l} = obj.section{s}.key{p}.name;
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
        function data = getData(obj, section, key)            
            % Get the value of a specified key
            if (~obj.getReadStatus() && isempty(obj.section))
               %obj.printWarning('File not yet read!\n');
               obj.readFile();
            end
            
            data = [];
            if (nargin == 2)
                key = section;
                % Search the key among all the sections
                s = 1;
                while ((s<=length(obj.section)) && (s > 0))
                    p = 1;
                    while ((p <= length(obj.section{s}.key)) && (p ~= 0))
                        if (strcmp(obj.section{s}.key{p}.name,key))
                            data = obj.section{s}.key{p}.data;
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
                    obj.logger.addWarning(['Key "' key '" not found while reading: "' obj.file_name '"'], 10);
                    data = [];
                end
            else                
                % Search the key of a specific section
                s = 1;
                while ((s<=length(obj.section)) && (s > 0))
                    if (strcmp(obj.section{s}.name,section))
                        p = 1;
                        while ((p <= length(obj.section{s}.key)) && (p ~= 0))
                            if (strcmp(obj.section{s}.key{p}.name,key))
                                data = obj.section{s}.key{p}.data;
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
                if (isempty(data))
                    obj.logger.addWarning(['Key "' key '" not found in section "' section '" while reading: "' obj.file_name '"'], 10);
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
        function addSection(obj, new_section)
            % Add a new section to the object
            if ~obj.isSection(new_section)
                newId = length(obj.section)+1;
                obj.section{newId}.name = new_section;
                obj.section{newId}.key = {};
            end
        end
        
        % Add a new keys to the object ------------------------------------
        function addKey(obj, section, key, data)
            % Add a new key to the object if the section exist
            if obj.isKey(section, key)
                obj.editKey(section, key, data);
            else
                s = 1;
                while ((s<=length(obj.section)) && (s ~= 0))
                    if (strcmp(obj.section{s}.name,section))
                        newId = length(obj.section{s}.key)+1;
                        obj.section{s}.key{newId}.name = key;
                        obj.section{s}.key{newId}.data = data;
                        s = 0;
                    else
                        s = s+1;    % go on with the search of the section
                    end
                end
                if (s ~= 0)
                    obj.printError(['Section "' section '" not found!\n']);
                end
            end
        end
        
        % Add a data to the object giving key and section -----------------
        function addData(obj, section, key, data)
            % Add a data to the object an eventually create the section and key to store it
            obj.addSection(section)
            obj.addKey(section, key, data)
        end

        % Remove a section from the object Ini_Manager --------------------
        function rmSection(obj, section)
            % Remove a section from the object Ini_Manager
            s = 1;
            while ((s<=length(obj.section)) && (s ~= 0))
                if (strcmp(obj.section{s}.name,section))
                    obj.section(s) = [];
                    s = 0;
                else
                    s = s+1;    % go on with the search of the section
                end
            end
            if (s ~= 0)
                obj.printError(['Section "' section '" not found!\n']);
            end
        end
        
        % Remove a key from the object Ini_Manager ------------------------
        function rmKey(obj, section, key)
            % Remove a key from the object Ini_Manager
            s = 1;
            k = 1;
            while ((s<=length(obj.section)) && (s ~= 0))
                if (strcmp(obj.section{s}.name,section))
                    k = 1;
                    while ((k<=length(obj.section{s}.key)) && (k ~= 0))
                        if (strcmp(obj.section{s}.key{k}.name,key))
                            obj.section{s}.key(k) = [];
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
                obj.printError(['Key "' key '" not found!\n']);
            end
        end
        
        % Edit a key in the object Ini_Manager ----------------------------
        function editKey(obj, section, key, data)
            % Edit a key in the object Ini_Manager
            s = 1;
            k = 1;
            while ((s<=length(obj.section)) && (s ~= 0))
                if (strcmp(obj.section{s}.name,section))
                    k = 1;
                    while ((k<=length(obj.section{s}.key)) && (k ~= 0))
                        if (strcmp(obj.section{s}.key{k}.name,key))
                            obj.section{s}.key{k}.data = data;
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
                obj.printError(['Key "' key '" not found!\n']);
            end
        end
        
    end
    
    % =========================================================================
    %  VISUALIZATION FUNCTIONS (not using LOGGER) - for DEBUG
    % =========================================================================
    methods
        
        % List Sections ---------------------------------------------------
        function listSections(obj, color_mode)
            % List the list of available sections
            if (nargin == 1)
                color_mode = obj.logger.getColorMode();
            end
            
            if (~obj.getReadStatus())
                %obj.printWarning('File not yet read!\n');

                obj.readFile();
            end
            
            fprintf('        List of sections:\n')
            for s=1:length(obj.section)
                obj.printSection(obj.section{s}.name, color_mode);
            end
        end
        
        % List Key --------------------------------------------------------
        function listKeys(obj, section, color_mode)
            % List the list of the keys available
            if (nargin == 1)
                section = 0;
                color_mode = obj.logger.getColorMode();
            end
            if (nargin == 2)
                if (ischar(section))
                    color_mode = obj.logger.getColorMode();
                else
                    color_mode = section;
                    section = 0;
                end                
            end
            
            if (~obj.getReadStatus() && isempty(obj.section))
                %obj.printWarning('File not yet read!\n');
                obj.readFile();
            end
            
            if (section == 0)
                fprintf('        List of all the keys:\n')
                for s = 1:length(obj.section)
                    obj.printSection(obj.section{s}.name, color_mode);
                    for p = 1:length(obj.section{s}.key)
                        obj.printData(obj.section{s}.key{p}.name, obj.section{s}.key{p}.data, color_mode)
                    end
                end
            else
                fprintf('        List of all the keys in section "%s":\n', section)
                s = 1;
                while ((s<=length(obj.section)) && (s ~= 0))
                    if (strcmp(obj.section{s}.name,section))
                        obj.printSection(obj.section{s}.name, color_mode);
                        for p = 1:length(obj.section{s}.key)
                            obj.printData(obj.section{s}.key{p}.name, obj.section{s}.key{p}.data, color_mode)
                        end
                        s = 0;      % Stop searching
                    else
                        s = s+1;    % go on with the search of the section
                    end
                end
            end
        end

        % Show data -------------------------------------------------------
        function showData(obj, section, color_mode)
            % List the data contained in the ini
            if (nargin == 1)
                section = 0;
                color_mode = obj.logger.getColorMode();
            end
            if (nargin == 2)
                if (ischar(section))
                    color_mode = obj.logger.getColorMode();
                else
                    color_mode = section;
                    section = 0;
                end                
            end            
            
            if (~obj.getReadStatus())
                %obj.printWarning('File not yet read')
                obj.readFile();
            end
            
            if (section == 0)
                fprintf('        List of all the key values:\n')
                for s = 1:length(obj.section)
                    obj.printSection(obj.section{s}.name, color_mode);
                    for p = 1:length(obj.section{s}.key)
                        obj.printData(obj.section{s}.key{p}.name, obj.section{s}.key{p}.data, color_mode)
                    end
                end
            else
                fprintf('        List of all the key values in section "%s":\n', section)
                s = 1;
                while ((s <= length(obj.section)) && (s ~= 0))
                    if (strcmp(obj.section{s}.name,section))
                        obj.printSection(obj.section{s}.name, color_mode);
                        for p = 1:length(obj.section{s}.key)
                            obj.printData(obj.section{s}.key{p}.name, obj.section{s}.key{p}.data, color_mode)                           
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
        function printSection(obj, section, color_mode)
            % Display a sections on a tree
            if (nargin == 2)
                color_mode = obj.color_mode;
            end
            if (color_mode)
                cprintf('blue','         - ');
                cprintf('text', '"%s"\n', section);
            else
                fprintf('         - "%s"\n', section);
            end
        end
        
        % Display a Key ---------------------------------------------------
        function printKey(obj, key, color_mode)
            % Display a key on a tree
            if (nargin == 2)
                color_mode = obj.color_mode;
            end
            if (color_mode)
                cprintf('blue','             |- ');
                cprintf('text', '"%s"\n', key);                
            else
                fprintf('             |- "%s"\n', key);
            end
        end
        
        % Display a Key ---------------------------------------------------
        function printData(obj, key, data, color_mode)
            % Display the data of a key
            if (nargin == 3)
                color_mode = obj.color_mode;
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
                        cprintf('text','"%s" = "%s"\n',key, num2str(tmpData));
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
        function cleanRaw(obj)
            % Strip empty lines
            obj.raw_data((cellfun(@length, obj.raw_data)==0)) = [];
            % Strip commented lines
            r = 1;
            while (r <= length(obj.raw_data))
                if ismember(obj.raw_data{r}(1), obj.c_comment)
                    obj.raw_data(r) = [];
                else
                    r = r + 1;
                end
            end
            
            for r=1:length(obj.raw_data)
                % Strip inline comments
                
                start = regexp(obj.raw_data{r}, ['([' obj.c_comment ']).*'], 'start');
                if (~isempty(start))
                    obj.raw_data{r}(start:end) = [];
                end
                
                start = regexp(obj.raw_data{r}, ' *$', 'start');
                if (~isempty(start))
                    obj.raw_data{r}(start:end) = [];
                end
            end
        end
        
        % Parse the file (cleaned) ----------------------------------------
        function parseData(obj)
            p = 0;
            s = 0;
            obj.section = {};
            for r=1:length(obj.raw_data)
                sectionName = regexp(obj.raw_data{r}, '(?<=^\[).*(?=\])', 'match');
                if (~isempty(sectionName))
                    % we have a new SECTION!
                    s = s+1;
                    p = 0;
                    obj.section{s}.name = sectionName{1};
                    obj.section{s}.key = [];
                else
                    % maybe this line contains a key, let's check
                    parName = regexp(obj.raw_data{r}, '^(.*(?=\s\=))', 'match');
                    if (isempty(parName))
                        parName = regexp(obj.raw_data{r}, '^(.*(?=\=))', 'match');
                    end
                    if (~isempty(parName))
                        % we have a PARAMETER!
                        p = p+1;
                        if (s == 0)
                            s = s + 1;
                            obj.section{s}.name = 'Section';
                            obj.section{s}.key = [];
                        end
                        obj.section{s}.key{p}.name = parName{1};
                        
                        % Get the DATA!
                        strData = regexp(obj.raw_data{r}, '(?<=\=).*', 'match');
                        if (isempty(strData))
                            obj.section{s}.key{p}.data = [];
                        else
                            tmpData = str2num(strData{1});
                            if (~isempty(tmpData))   % It's a number!
                                obj.section{s}.key{p}.data = tmpData;
                            else
                                % It's a string!
                                % If it is a string properly formatted "string data"
                                tmpData = regexp(strData{1}, '(?<=\")[^"]*(?=\")', 'match');
                                if (~isempty(tmpData))
                                    if (size(tmpData,2) == 1)
                                        obj.section{s}.key{p}.data = tmpData{1};
                                    else                                        
                                        % Strip empty cells
                                        tmpData((cellfun(@length, tmpData)==0)) = [];
                                        % Strip space only cells
                                        for c=size(tmpData,2):-1:1
                                            if isempty(regexp(tmpData{c}, '[^ ]*', 'match'))
                                                tmpData(c) = [];
                                            end
                                        end
                                        obj.section{s}.key{p}.data = tmpData;
                                    end
                                else
                                    % If it is not a string but apparently is
                                    % not a number
                                    obj.section{s}.key{p}.data = strData{1};
                                end
                            end
                        end
                    end
                    
                end
                
            end
        end 
        
        % Set read status -------------------------------------------------
        function setReadStatus(obj, bool)
            obj.read_status = bool;
        end        
    end
    
    % =========================================================================
    %  STATIC STRING GENERATION FROM INI TO STR
    % =========================================================================
    
    methods (Static)
        
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
            elseif isnumeric(value) % is string
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
                toString = @(var) strtrim(regexprep(evalc(['disp(var)']), '\n', ''));
                if iscell(value)
                    cell_str{numel(cell_str) + 1} = [variable_name ' = [' toString(value) ']'];
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

            cell_str{numel(cell_str) + 1} = [Ini_Manager.STD_COMMENT '--------------------------------------------------------------------------------'];
            cell_str{numel(cell_str) + 1} = ['[' section_name ']'];
            cell_str{numel(cell_str) + 1} = [Ini_Manager.STD_COMMENT '--------------------------------------------------------------------------------'];

            % I want a column array
            if size(cell_str,1) < size(cell_str,2)
                cell_str = cell_str';
            end
        end
        
        % toIniString Comment ---------------------------------------------
        function cell_str = toIniStringComment(comment, cell_str)
            % Add a comment in ini string format
            % SYNTAX: 
            %   cell_str = toIniStringSection(section_name, cell_str)
            if (nargin == 1)
                cell_str = {};
            end            
            cell_str{numel(cell_str) + 1} = [ Ini_Manager.STD_COMMENT ' ' comment];
            
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
