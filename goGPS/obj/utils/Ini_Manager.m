% Ini_Manager(fileName, verbosity)
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
%   cprintf:    http://www.mathworks.com/matlabcentral/fileexchange/24093-cprintf-display-formatted-colored-text-in-the-command-window
%
%----------------------------------------------------------------------------------------------
%                           goGPS v0.5.9
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
% Written by:       Gatti Andrea
% Contributors:     Gatti Andrea, ...
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
classdef Ini_Manager < Ini_Reader
      
    properties (GetAccess = 'private', SetAccess = 'private')
    end
    
    % =========================================================================
    %    PUBLIC METHODS
    % =========================================================================
    
    methods
        
        % Creator
        function obj = Ini_Manager(file_name, verbosity)
        % Ini_Manager(fileName, verbosity) creator
            if (nargin == 0)
                file_name = '';
            end
            if isempty(file_name)
                file_name = '';
            end
            obj.setFileName(file_name);
            if (nargin == 2)
                obj.setVerbosityLev(verbosity);
            end            
        end
        
        % Distructor
        function delete(obj)
            obj.delete@Ini_Reader();            
        end
        
        % Cleaner
        function clearAll(obj)
            obj.clearAll@Ini_Reader();

            % Reset the object
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

        % Remove a section from the object IniReader ----------------------
        function rmSection(obj, section)
            % Remove a section from the object IniReader
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
        
        % Remove a key from the object IniReader ----------------------
        function rmKey(obj, section, key)
            % Remove a key from the object IniReader
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
        
        % Edit a key in the object IniReader ----------------------
        function editKey(obj, section, key, data)
            % Edit a key in the object IniReader
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
    %  STATIC STRING GENERATION FROM INI TO STR
    % =========================================================================
    
    methods (Static)
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
                cell_str{numel(cell_str) + 1} = [variable_name ' = "' value '"'];
            elseif isnumeric(value) % is string
                if isempty(format)
                    format = '%g';
                else
                    format = strtrim(format);
                end
                if numel(value) > 1 % it's an array of values
                    format = [format ' '];
                    cell_str{numel(cell_str) + 1} = [variable_name ' = [' sprintf(format,value) ']'];
                else
                    cell_str{numel(cell_str) + 1} = [variable_name ' = ' sprintf(format,value)];
                end
            else % generic converter (may not work properly)
                toString = @(var) evalc('disp(var)');
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
    end
end
