%   CLASS Settings_Interface
% =========================================================================
%
% DESCRIPTION
%   Abstract class with basic properties and methods that a setting class
%   must have
%
% COMMENTS
% Settings have been build with multiple inheritance
% A standard abstract interface have been created: Setting_Interface
% it add to each subclass the object log
% force the subclasses to implement three basic methods:
%  - import:       when a setting object of identical class, or that inherits from the same class is passed to this function, the relative parameters are copied in the calling object
%                  when an ini file is the input of the function, the object is updated with the settings contained into the file
%  - toString:     display the content of the object, in a human readable way, a goGPS user can "ask" for the value of the settings on screen
%  - export:       create a cell array of strings containing the settings in plain text ini format. The variable it's the raw data format of Ini_Manager
%

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.2 beta 1
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

classdef Settings_Interface < handle
    properties (SetAccess = protected, GetAccess = protected)
        log; % Handler to the log object
    end

    properties (Abstract)
    end

    methods  (Abstract)
        import(obj, settings);
        % This function import the settings of the current class from another setting object having the same properties, or an ini file

        str = toString(obj, str)
        % Display the content of the class in a human readable format

        str_cell = export(obj, str_cell)
        % Conversion to string ini format of the minimal information needed to reconstruct the obj
    end

    methods (Access = 'public')
        function initLogger(this)
            % Init the log object
            % SYNTAX: this.initLogger();
            this.log = Logger.getInstance();
        end

        function ini = save(this, file_path)
            % Save to a file (in INI fomat) the content of the Settings object
            % SYNTAX: <ini> = this.save(file_path);
            % return optionally the ini manager object used by the save function
            ini = Ini_Manager(file_path, this.export());
        end

        function importIniFile(this, file_path)
            % Import from an INI file the content of the Settings object
            % SYNTAX: this.importIniFile(file_path);
            ini = Ini_Manager(file_path);
            fnp = File_Name_Processor;
            prj_home = fnp.checkPath(ini.getData('prj_home'));
            if ~isempty(prj_home) && ~exist(prj_home, 'dir')
                [path_str, ~, ~] = fileparts(fnp.getFullDirPath(file_path));
                prj_home = fnp.getFullDirPath(strcat(path_str, [filesep '..' filesep]), pwd);
                if exist(prj_home, 'dir')
                    ini.setData('prj_home', prj_home);
                end
            end
            this.import(ini);
        end

    end

    % =========================================================================
    %  TEST PARAMETERS VALIDITY
    % =========================================================================
    methods (Access = 'protected')

        function checkLogicalField(this, field_name)
            % Check if a logical field of the object is a valid logical number
            % To make the function works it is needed to have defined the default
            % value of the field as a constant with same name but upper case
            % This superclass function must be copied in each child (to
            % have read / write permission on the parameters
            % SYNTAX: this.checkLogicalField(string_field_name);
            this.(field_name) = this.checkLogical(field_name, this.(field_name), this.(upper(field_name)));
        end

        function checkStringField(this, field_name, empty_is_valid, check_existence)
            % Check if a string field of the object is a valid string
            % To make the function works it is needed to have defined the default
            % value of the field as a constant with same name but upper case
            % This superclass function must be copied in each child (to
            % have read / write permission on the parameters
            % SYNTAX: this.checkStringField(string_field_name, <empty_is_valid == false>, <check_existence == false>);

            switch nargin
                case 2, this.(field_name) = this.checkString(field_name, this.(field_name), this.(upper(field_name)));
                case 3, this.(field_name) = this.checkString(field_name, this.(field_name), this.(upper(field_name)), empty_is_valid);
                case 4, this.(field_name) = this.checkString(field_name, this.(field_name), this.(upper(field_name)), empty_is_valid, check_existence);
                otherwise, error('Settings checkStringField called with the wrong number of parameters');
            end
        end

        function checkNumericField(this, field_name, limits, valid_val)
            % Check if a numeric field of the object is valid
            % To make the function works it is needed to have defined the default
            % value of the field as a constant with same name but upper case
            % This superclass function must be copied in each child (to
            % have read / write permission on the parameters
            % SYNTAX: this.checkNumericField(string_field_name, <limits>, <valid_values>);
            switch nargin
                case 2, this.(field_name) = this.checkNumber(field_name, this.(field_name), this.(upper(field_name)));
                case 3, this.(field_name) = this.checkNumber(field_name, this.(field_name), this.(upper(field_name)), limits);
                case 4, this.(field_name) = this.checkNumber(field_name, this.(field_name), this.(upper(field_name)), limits, valid_val);
                otherwise, error('Settings checkNumericField called with the wrong number of parameters');
            end
        end

        function checked_val = checkLogical(this, field_name, field_val, default_val)
            % Check if a logical is a valid logical
            % This superclass function must be called in each child
            % SYNTAX: checked_val = this.checkLogicalField(string_field_name);
            checked_val = default_val;
            if (~isempty(field_val)) && (~isnan(field_val))
                checked_val = logical(field_val);
            else
                this.log.addWarning(sprintf('The settings field %s is not valid => using default %d', field_name, checked_val));
            end
        end

        function checked_val = checkCellString(this, field_name, field_val, default_val, empty_is_valid, check_existence)
            % Check if a string is a valid string
            % SYNTAX: checked_val = this.checkString(string_field_name, <empty_is_valid == false>, <check_existence == false>);
            if (nargin < 5)
                empty_is_valid = false;
            end
            if (nargin < 6)
                check_existence = false;
            end

            checked_val = this.checkString(field_name, field_val, default_val, empty_is_valid, check_existence);
            if ~iscell(checked_val)
                checked_val = {checked_val};
            end
        end

        function checked_val = checkString(this, field_name, field_val, default_val, empty_is_valid, check_existence)
            % Check if a string is a valid string
            % SYNTAX: checked_val = this.checkString(string_field_name, <empty_is_valid == false>, <check_existence == false>);
            if (nargin < 5)
                empty_is_valid = false;
            end
            if (nargin < 6)
                check_existence = false;
            end

            checked_val = default_val;
            default_val = field_val;
            if ~isempty(field_val) && iscell(field_val) % A cell of strings must contain at least one string
                field_val = field_val{1};
            end
            if (ischar(field_val) || (empty_is_valid && isempty(field_val))) && ((~isempty(field_val)) || empty_is_valid) && ((exist(field_val,'file') || exist(field_val,'dir')) || ~check_existence)
                checked_val = default_val;
            else
                if iscell(checked_val)
                    this.log.addWarning(sprintf('The value "%s" of the settings field %s is not valid => using default "%s"', field_val, field_name, Ini_Manager.strCell2Str(checked_val)));
                else
                    this.log.addWarning(sprintf('The value "%s" of the settings field %s is not valid => using default "%s"', field_val, field_name, checked_val));
                end
            end
        end

        function checked_val = checkNumber(this, field_name, field_val, default_val, limits, valid_val)
            % Check if a number is valid
            % SYNTAX: checked_val = this.checkNumericField(string_variable_name, value, default_value, <limits>, <valid_values>);
            checked_val = default_val;
            if isnumeric(field_val) && (~isempty(field_val)) && (any(~isnan(field_val)))
                checked_val = field_val;
                % if I have limits to check => check for out of bound
                for i = 1 : numel(checked_val)
                    if (nargin >= 5) && (numel(limits) == 2) && ...
                            ((checked_val(i) > limits(2)) || (checked_val(i) < limits(1)))
                        checked_val(i) = max(limits(1), min(limits(2), checked_val(i)));
                        this.log.addWarning(sprintf('The value %g of the settings field %s is not within the valid limits (%g .. %g) => updating it to %g', field_val(i), field_name, limits(1),limits(2), checked_val(i)));
                    end
                    % if I have a set of values => check for set intersection
                    if (nargin >= 6) && (~isempty(valid_val)) && (~ismember(checked_val(i), valid_val))
                        checked_val(i) = default_val;
                        this.log.addWarning(sprintf('The value %g for the settings field %s is not valid => using default %g. It should be one of: %s', field_val(i), field_name, checked_val(i), sprintf('%g ', valid_val)));
                    end
                end
            else
                this.log.addWarning(sprintf('The settings field %s is not valid => using default %g', field_name, checked_val));
            end
        end
    end

    % =========================================================================
    %  TEST
    % =========================================================================

    methods (Access = 'protected')
        function testInterfaceRoutines(this)
            % test the class (Interface Routines)
            % SINTAX: this.testInterfaceRoutines();

            %try
                vl = this.log.getVerbosityLev();
                this.log.setVerbosityLev(1e3);
                test = this;
                raw_data = test.export();

                ini = Ini_Manager('test__.ini', raw_data);
                ini.showData();

                test.import(ini);
                test_copy = repmat(test,2,1); test_copy = test_copy(2);
                test.import(test_copy);
                clear test_copy;
                fprintf('\n');
                disp(test.toString());
                %delete('test__.ini');
            %catch ex
            %    this.log.addError(['Test failed: ' ex.message]);
            %end
            this.log.setVerbosityLev(vl);
        end
    end
end
