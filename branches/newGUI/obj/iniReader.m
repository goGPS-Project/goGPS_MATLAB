% =========================================================================
%   OBJECT IniReader
% =========================================================================
%
% DESCRIPTION:
%   Object to read an ini file
%
% INI EXAMPLE:
%   [Section 1]
%       array = [1 2 3]
%       string = "I'm a string"
%       stringCellArray = ["I" "am" "a string cell" "array"]
%   ; I'm a commen
%   # I'm also a comment
%   [Section 2]
%       number = 10
%   [Section n]     # comment
%       double = 10.4
%
% REQUIRES:
%   cprintf:    http://www.mathworks.com/matlabcentral/fileexchange/24093-cprintf-display-formatted-colored-text-in-the-command-window
%
% LIST of METHODS
%
% iniReader(fileName, verbosity)                creator
% delete(obj)                                   distructor (empty)
%
% FILE ----------------------------------------------------------------
%
%  setFileName(obj, fileName)                   set the complete path of the ini file
%  fileName = getFileName(obj)                  get the complete path of the ini file
%
% INI SYNTAX ----------------------------------------------------------
%
%  setCommentChar(obj, character)               define the characters that define comment lines
%
% R/W MODES OF THE FILE -----------------------------------------------
%
%  setRW(obj,rwMode)                            set mode for open the ini file
%  rwMode = getRW(obj)                          get the current mode of the file
%  bool = readMode(obj)                         is it in read mode?
%  bool = writeMode(obj)                        is it in write mode?
%
%  bool = getReadStatus(obj)                    return if the file has been already parsed
%
%  readFile(obj)                                force reading of the File
%  update(obj, filename, force)                 update the object when needed:
%                                                - filename changed 
%                                                - force flag == 1
%                                                - INI not yet read
%
% MANAGING LOGGING ----------------------------------------------------
%
%  setColorMode(obj, bool)                      set useage of colors in text output
%  bool = getColorMode(obj)                     get useage of colors in text output
%  setVerbosityLev(obj, verbosity)              set level of verbosity
%  verbosity = getVerbosityLev(obj)             get level of verbosity
%
% GETTER OF THE PARAMETERS --------------------------------------------
%  
%  isS = isSection(obj, section)                get the presence of a section 
%  isK = isKey(obj, section, key)               get the presence of a key 
%  sectionList = getSections(obj)               get the list of available sections
%  keyList = getKeys(obj, <section>)            get the list of the keys available
%  data = getData(obj, <section>, key)          get the value of a specified key
%
% MODIFIER FUNCTIONS --------------------------------------------------
%
%  addSection(obj, newSection)                  add a new section to the object
%  addKey(obj, section, key, data)              add a new keys to the object
%  rmSection(obj, section)                      remove a section from the object iniReader
%  rmKey(obj, section, key)                     remove a key from the object iniReader
%  editKey(obj, section, key, data)             edit a key in the object iniReader
%
% CONTAIN FUNCTION ----------------------------------------------------
%
% isPresent = containsSection(obj, section)     search a section in the object iniReader
% isPresent = containsKey(obj, section, key)    search a key in the object iniReader
%
% VISUALIZATION FUNCTIONS ---------------------------------------------
%
%  listSections(obj, <colorMode>)               list the list of available sections
%  listKeys(obj, <section>, <colorMode>)        list the list of the keys available
%  showData(obj, <section>, <colorMode>)        list the data contained in the ini
%
%
%----------------------------------------------------------------------------------------------
% Copyright (C) 2009-2013 Mirko Reguzzoni, Eugenio Realini
%----------------------------------------------------------------------------------------------
%
%    Code contributed by Andrea Gatti
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
%---------------------------------------------------------------------------------------------

classdef iniReader < handle
      
    properties (GetAccess = 'private', SetAccess = 'private')
        chrComment  = [';', '#']; % style of comments
        strFileName = 'config.ini';
        fid;
        rawData = {};
        rw = 'r';
        readStatus = false;
        verbosity = 1;
        colorMode = 1;
    
        section = {}
        %section{1}.name = 'Section Name';
        %section{1}.par = {};
        %section{1}.par{1}.name = 'Key name';
        %section{1}.par{1}.data = 100;
    end
    
    % =========================================================================
    %    PUBLIC METHODS
    % =========================================================================
    methods
        
        % Creator
        function obj = iniReader(fileName, verbosity)
            if (nargin < 1)
                fileName = '';
            end
            if isempty(fileName)
                fileName = '';
            end
            obj.setFileName(fileName);
            if (nargin == 2)
                obj.setVerbosityLev(verbosity);
            end            
        end
        
        % Distructor
        function delete(obj)
            
        end
        
        % Cleaner
        function clearAll(obj)
            obj.fid = 0;
            obj.rw = 'r';
            obj.strFileName = 'config.ini';
            obj.rawData = {};
            obj.verbosity = 1;
            obj.colorMode = 1;
            obj.section = {};
        end
        
        % ======================================================================
        %    Getter / Setter
        % ======================================================================
        
        % Set file name ---------------------------------------------------
        function setFileName(obj, fileName)
            obj.strFileName = fileName;
        end
        
        function fileName = getFileName(obj)
            fileName = obj.strFileName;
        end
        
        % Set Comment Character -------------------------------------------
        function setCommentChar(obj, character)
            obj.chrComment = character;
        end
        
        function characters = getCommentChar(obj)
            characters = obj.chrComment;
        end
        
        % Set RW mode -----------------------------------------------------
        function setRW(obj,rwMode)
            if ( ~obj.readMode() && ~writeMode())
                obj.rw = 'r';
            else
                obj.rw = rwMode;
            end
        end
        
        function rwMode = getRW(obj)
            rwMode = obj.rw;
        end
        
        function bool = readMode(obj)
            if ((obj.rw(1)=='r'))
                bool = true;
            else
                bool = false;
            end
        end
        
        function bool = writeMode(obj)
            if ((obj.rw(1)=='w'))
                bool = true;
            else
                bool = false;
            end
        end
        
        % Get read status mode --------------------------------------------        
        function bool = getReadStatus(obj)
            bool = obj.readStatus;
        end
        
        % Read File -------------------------------------------------------
        function errStatus = readFile(obj)
            errStatus = false;            
            % If the object already contains data - clean it
            if (obj.getReadStatus())
                obj.rawData = {};
            end
            if ~exist(obj.getFileName(),'file')
                obj.printError('File not found');
                errStatus = true;
            else
                obj.fid = fopen(obj.getFileName(), obj.getRW());
                
                if (obj.fid ~= -1)   % If reading is ok
                    obj.rawData = textscan(obj.fid, '%s', 'delimiter', '\n', 'endOfLine', '\r\n');
                    
                    fclose(obj.fid);
                    
                    obj.readStatus = true;
                    obj.rawData = obj.rawData{1};
                    
                    if (obj.getVerbosityLev)
                        obj.opStatus(1);
                        cprintf('The INI file has been read correctly.\n');
                    end
                    
                    obj.cleanRaw();      % Stirp comments and spaces
                    obj.parseData();     % Parse file
                    obj.rawData = {};    % clean RAW temp data
                else
                    obj.printError(ferror(obj.fid));
                    obj.setReadStatus(false);
                    obj.rawData = {};
                end
            end
        end
        
        % Update File (when needed) ---------------------------------------
        function reloaded = update(obj, filename, force)
            if nargin == 2
                force = 0;
            end
            reloaded = 0;
            if (~strcmp(filename,obj.getFileName) || (force == 1) || ~obj.getReadStatus())
                obj.setFileName(filename);
                obj.readFile();
                reloaded = 1;
            end            
        end
                
        % Set read status mode --------------------------------------------
        function setColorMode(obj, bool)
            obj.colorMode = bool;
        end
        
        function bool = getColorMode(obj)
            bool = obj.colorMode;
        end
        
        % Set verbosity level ---------------------------------------------
        function setVerbosityLev(obj, verbosity)
            obj.verbosity = verbosity;
        end
        
        function verbosity = getVerbosityLev(obj)
            verbosity = obj.verbosity;
        end
        
        % Get the list of sections present in the file --------------------        
        function sectionList = getSections(obj)
            if (~obj.getReadStatus())
                obj.printWarning('File not yet read!\n');
                obj.readFile();
            end
            
            sectionList = {};
            for s=1:length(obj.section)
                sectionList{s} = obj.section{s}.name;
            end
        end 
        
        % Return the presence of a section --------------------------------
        function isS = isSection(obj, section)
            s = 1;
            while ((s<=length(obj.section)) && (s ~= 0))
                if (strcmp(obj.section{s}.name,section))
                    s = 0;      % stop searching
                else
                    s = s+1;    % go on with the search of the section
                end
            end
           isS = s == 0;
        end
        
        % Return the presence of a key ------------------------------------
        function isK = isKey(obj, section, key)
            s = 1;
            while ((s<=length(obj.section)) && (s ~= 0))
                if (strcmp(obj.section{s}.name,section))
                    k = 1;
                    while ((k<=length(obj.section{s}.par)) && (k ~= 0))
                        if (strcmp(obj.section{s}.par{k}.name,key))
                            obj.section{s}.key{k} = [];
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
            if (nargin == 1)
                section = 0;
            end
            
            if (~obj.getReadStatus())                
                obj.printWarning('File not yet read!\n');
                obj.readFile();
            end
            
            keyList = {};
            l = 1;
            if (section == 0)
                for s = 1:length(obj.section)
                    for p = 1:length(obj.section{s}.par)
                        keyList{l} = obj.section{s}.par{p}.name;
                        l = l+1;
                    end
                end
            else
                s = 1;
                while ((s<=length(obj.section)) && (s ~= 0))
                    if (strcmp(obj.section{s}.name,section))
                        for p = 1:length(obj.section{s}.par)
                            keyList{l} = obj.section{s}.par{p}.name;
                            l = l+1;
                        end
                        s = 0;      % stop searching
                    else
                        s = s+1;    % go on with the search of the section
                    end
                end
            end
        end
        
        % Get data --------------------------------------------------------        
        function data = getData(obj, section, key)
            if (~obj.getReadStatus())
                obj.printWarning('File not yet read!\n');
                obj.readFile();
            end
            
            data = [];
            if (nargin == 2)
                key = section;
                % search the key among all the sections
                s = 1;
                while ((s<=length(obj.section)) && (s > 0))
                    p = 1;
                    while ((p <= length(obj.section{s}.par)) && (p ~= 0))
                        if (strcmp(obj.section{s}.par{p}.name,key))
                            data = obj.section{s}.par{p}.data;
                            p = 0;      % stop searching key
                        else
                            p = p+1;    % go on with the search of the key
                        end
                    end
                    if (p == 0)
                        s = 0;     % stop searching section
                    else
                        s = s+1;   % go on with the search of the section
                    end
                end
            else
                % search the key of a specific section
                s = 1;
                while ((s<=length(obj.section)) && (s > 0))
                    if (strcmp(obj.section{s}.name,section))
                        p = 1;
                        while ((p <= length(obj.section{s}.par)) && (p ~= 0))
                            if (strcmp(obj.section{s}.par{p}.name,key))
                                data = obj.section{s}.par{p}.data;
                                p = 0; % stop searching key
                            else
                                p = p+1;    % go on with the search of the key
                            end
                        end
                        s = 0;      % stop searching section
                    else
                        s = s+1;    % go on with the search of the section
                    end
                end
            end
            
            if (isempty(data))
                obj.printError(['Key "' key '" not found!\n']);
                data = [];
            end
        end
           
        % ======================================================================
        %    Modifiers functions
        % ======================================================================
        
        % Add a new section to the object ---------------------------------        
        function addSection(obj, newSection)
            if ~obj.isSection(newSection)
                newId = length(obj.section)+1;
                obj.section{newId}.name = newSection;
                obj.section{newId}.par = {};
            end
        end
        
        % Add a new keys to the object ------------------------------------
        function addKey(obj, section, key, data)
            if obj.isKey(section, key)
                obj.editKey(section, key, data);
            else
                s = 1;
                while ((s<=length(obj.section)) && (s ~= 0))
                    if (strcmp(obj.section{s}.name,section))
                        newId = length(obj.section{s}.par)+1;
                        obj.section{s}.par{newId}.name = key;
                        obj.section{s}.par{newId}.data = data;
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
        
        % Remove a section from the object iniReader ----------------------
        function rmSection(obj, section)
            s = 1;
            while ((s<=length(obj.section)) && (s ~= 0))
                if (strcmp(obj.section{s}.name,section))
                    obj.section{s} = [];
                    s = 0;
                else
                    s = s+1;    % go on with the search of the section
                end
            end
            if (s ~= 0)
                obj.printError(['Section "' section '" not found!\n']);
            end
        end
        
        % Remove a key from the object iniReader ----------------------
        function rmKey(obj, section, key)
            s = 1;
            while ((s<=length(obj.section)) && (s ~= 0))
                if (strcmp(obj.section{s}.name,section))
                    while ((k<=length(obj.section{s}.par)) && (k ~= 0))
                        if (strcmp(obj.section{s}.par{k}.name,key))
                            obj.section{s}.key{k} = [];
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
        
        % Edit a key in the object iniReader ----------------------
        function editKey(obj, section, key, data)
            s = 1;
            while ((s<=length(obj.section)) && (s ~= 0))
                if (strcmp(obj.section{s}.name,section))
                    k = 1;
                    while ((k<=length(obj.section{s}.par)) && (k ~= 0))
                        if (strcmp(obj.section{s}.par{k}.name,key))
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
        
        % ======================================================================
        %    Search functions
        % ======================================================================
        
        % Search a section in the object iniReader ------------------------
        function isPresent = containsSection(obj, section)
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
        
        % Search a key in the object iniReader ----------------------------
        function isPresent = containsKey(obj, section, key)
            s = 1;
            while ((s<=length(obj.section)) && (s ~= 0))
                if (strcmp(obj.section{s}.name,section))
                    while ((k<=length(obj.section{s}.par)) && (k ~= 0))
                        if (strcmp(obj.section{s}.par{k}.name,key))
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
        
        
        % ======================================================================
        %    Visualization functions
        % ======================================================================
        
        % List Sections ---------------------------------------------------
        function listSections(obj, colorMode)
            if (nargin == 1)
                colorMode = obj.getColorMode();
            end
            
            if (~obj.getReadStatus())
                obj.printWarning('File not yet read!\n');

                obj.readFile();
            end
            
            fprintf('        List of sections:\n')
            for s=1:length(obj.section)
                obj.printSection(obj.section{s}.name, colorMode);
            end
        end
        
        % List Key --------------------------------------------------------
        function listKeys(obj, section, colorMode)
            if (nargin == 1)
                section = 0;
                colorMode = obj.getColorMode();
            end
            if (nargin == 2)
                if (ischar(section))
                    colorMode = obj.colorMode;
                else
                    colorMode = section;
                    section = 0;
                end                
            end
            
            if (~obj.getReadStatus())
                obj.printWarning('File not yet read!\n');
                obj.readFile();
            end
            
            if (section == 0)
                fprintf('        List of all the keys:\n')
                for s = 1:length(obj.section)
                    obj.printSection(obj.section{s}.name, colorMode);
                    for p = 1:length(obj.section{s}.par)
                        obj.printData(obj.section{s}.par{p}.name, obj.section{s}.par{p}.data, colorMode)
                    end
                end
            else
                fprintf('        List of all the keys in section "%s":\n', section)
                s = 1;
                while ((s<=length(obj.section)) && (s ~= 0))
                    if (strcmp(obj.section{s}.name,section))
                        obj.printSection(obj.section{s}.name, colorMode);
                        for p = 1:length(obj.section{s}.par)
                            obj.printData(obj.section{s}.par{p}.name, obj.section{s}.par{p}.data, colorMode)
                        end
                        s = 0;      % stop searching
                    else
                        s = s+1;    % go on with the search of the section
                    end
                end
            end
        end

        % Show data -------------------------------------------------------
        function showData(obj, section, colorMode)
            if (nargin == 1)
                section = 0;
                colorMode = obj.colorMode;
            end
            if (nargin == 2)
                if (ischar(section))
                    colorMode = obj.colorMode;
                else
                    colorMode = section;
                    section = 0;
                end                
            end            
            
            if (~obj.getReadStatus())
                obj.printWarning('File not yet read')
                obj.readFile();
            end
            
            if (section == 0)
                fprintf('        List of all the key values:\n')
                for s = 1:length(obj.section)
                    obj.printSection(obj.section{s}.name, colorMode);
                    for p = 1:length(obj.section{s}.par)
                        obj.printData(obj.section{s}.par{p}.name, obj.section{s}.par{p}.data, colorMode)
                    end
                end
            else
                fprintf('        List of all the key values in section "%s":\n', section)
                s = 1;
                while ((s <= length(obj.section)) && (s ~= 0))
                    if (strcmp(obj.section{s}.name,section))
                        obj.printSection(obj.section{s}.name, colorMode);
                        for p = 1:length(obj.section{s}.par)
                            obj.printData(obj.section{s}.par{p}.name, obj.section{s}.par{p}.data, colorMode)                           
                        end
                        s = 0;      % stop searching
                    else
                        s = s+1;    % go on with the search of the section
                    end
                end
            end
        end
        
    end
    
    % =========================================================================
    %    PRIVATE METHODS
    % =========================================================================
    
    methods (Access = 'private')
        
        % strip empty/comment lines
        function cleanRaw(obj)
            % strip empty lines
            obj.rawData((cellfun(@length, obj.rawData)==0)) = [];
            % strip commented lines
            r = 1;
            while (r <= length(obj.rawData))
                if ismember(obj.rawData{r}(1), obj.chrComment);
                    obj.rawData(r) = [];
                else
                    r = r + 1;
                end
            end
            
            for r=1:length(obj.rawData)
                % strip comments
                
                start = regexp(obj.rawData{r}, ['([' obj.chrComment ']).*'], 'start');
                if (~isempty(start))
                    obj.rawData{r}(start:end) = [];
                end
                
                start = regexp(obj.rawData{r}, ' *$', 'start');
                if (~isempty(start))
                    obj.rawData{r}(start:end) = [];
                end
            end
        end
        
        % Parse the file (cleaned)
        function parseData(obj)
            s = 0;
            for r=1:length(obj.rawData)
                sectionName = regexp(obj.rawData{r}, '(?<=^\[).*(?=\])', 'match');
                if (~isempty(sectionName))
                    % we have a new SECTION!
                    s = s+1;
                    p = 0;
                    obj.section{s}.name = sectionName{1};
                    obj.section{s}.par = [];
                else
                    % maybe this line contains a key, let's check
                    parName = regexp(obj.rawData{r}, '^(.*(?=\s\=))', 'match');
                    if (isempty(parName))
                        parName = regexp(obj.rawData{r}, '^(.*(?=\=))', 'match');
                    end
                    if (~isempty(parName))
                        % we have a PARAMETER!
                        p = p+1;
                        obj.section{s}.par{p}.name = parName{1};
                        
                        % GET the DATA!
                        strData = regexp(obj.rawData{r}, '(?<=\=).*', 'match');
                        if (isempty(strData))
                            obj.section{s}.par{p}.data = [];
                        else
                            tmpData = str2num(strData{1});
                            if (~isempty(tmpData))   % It's a number!
                                obj.section{s}.par{p}.data = tmpData;
                            else
                                % It's a string!
                                % If it is a string properly formatted "string data"
                                tmpData = regexp(strData{1}, '(?<=\")[^"]*(?=\")', 'match');
                                if (~isempty(tmpData))
                                    if (size(tmpData,2) == 1)
                                        obj.section{s}.par{p}.data = tmpData{1};
                                    else                                        
                                        % strip empty cells
                                        tmpData((cellfun(@length, tmpData)==0)) = [];
                                        % strip space only cells
                                        for c=size(tmpData,2):-1:1
                                            if isempty(regexp(tmpData{c}, '[^ ]*', 'match'))
                                                tmpData(c) = [];
                                            end
                                        end
                                        obj.section{s}.par{p}.data = tmpData;
                                    end
                                else
                                    % If it is not a string but apparently is
                                    % not a number
                                    obj.section{s}.par{p}.data = strData{1};
                                end
                            end
                        end
                    end
                    
                end
                
            end
        end 
        
        % Set read status
        function setReadStatus(obj, bool)
            obj.readStatus = bool;
        end

        
% =========================================================================
%    DISPLAY UTILITIES
% =========================================================================

        % Display Warning
        function printWarning(obj, text, colorMode)
            if (nargin == 2)
                colorMode = obj.colorMode;
            end
            if (obj.getVerbosityLev > 0)
                obj.opstatus(0);
                if (colorMode)
                    cprintf('err', 'warning: ');
                    cprintf('text', [text '\n']);
                else
                    fprintf('warning: %s\n', text);
                end
            end
        end
        
        % Display Error
        function printError(obj, text, colorMode)
            if (nargin == 2)
                colorMode = obj.colorMode;
            end
            if (obj.getVerbosityLev > 0)
                obj.opStatus(0);
                if (colorMode)
                    cprintf('err', 'Error: ');
                    cprintf('text', [text '\n']);
                else
                    fprintf('Error: %s\n', text);
                end
            end
        end
        
        % Display a Section
        function printSection(obj, section, colorMode)
            if (nargin == 2)
                colorMode = obj.colorMode;
            end
            if (colorMode)
                cprintf('blue','         - ');
                cprintf('text', '"%s"\n', section);
            else
                fprintf('         - "%s"\n', section);
            end
        end
        
        % Display a Key
        function printKey(obj, key, colorMode)
            if (nargin == 2)
                colorMode = obj.colorMode;
            end
            if (colorMode)
                cprintf('blue','             |- ');
                cprintf('text', '"%s"\n', key);                
            else
                fprintf('             |- "%s"\n', key);
            end
        end
        
        % Display a Key
        function printData(obj, key, data, colorMode)
            if (nargin == 3)
                colorMode = obj.colorMode;
            end
            if (colorMode)
                cprintf('blue','             |- ');
                tmpData = data;
                if (iscell(tmpData))    % the only dataset that can be a cell array is an array of strings
                    cprintf('text','"%s" = [', key);
                    for (l = 1:size(tmpData,2))
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
                if (iscell(tmpData))    % the only dataset that can be a cell array is an array of strings
                    fprintf('             |- "%s" = [', key);
                    for l = 1:size(tmpData,2)
                        fprintf(' "%s"', tmpData{l});
                    end
                    fprintf(' ]\n');
                else
                    % if it is not an array of string...
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
%    DISPLAY UTILITIES
% =========================================================================

    methods(Static, Access = 'private')
        
        % Display a flag of operation status
        function opStatus(statusOk)
            cprintf('blue',' [ ');
            if (statusOk)
                cprintf('green','ok');
            else
                cprintf('red','!!');
            end
            cprintf('blue',' ] ');
        end        
    end
end
