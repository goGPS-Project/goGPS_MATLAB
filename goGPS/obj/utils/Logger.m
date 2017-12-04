%   SINGLETON CLASS Logger
% =========================================================================
%
% DESCRIPTION
%   Manage messages output depending on the current interface
%   (up to now, only Command Windows is implemented)
%
% EXAMPLE
%   log = Logger.getInstance();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Logger

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

classdef Logger < handle

    properties (Constant, GetAccess = public)
        DEFAULT_VERBOSITY_LEV   = 9;  % Default verbosity level
        WARNING_VERBOSITY_LEV   = 3;  % Warning verbosity level
        ERROR_VERBOSITY_LEV     = 1;  % Error verbosity level
    end

    properties (GetAccess = 'private', SetAccess = 'protected')
        color_mode = true;            % Flag for coloured output messages (if true requires cprintf)
        verbosity = Logger.DEFAULT_VERBOSITY_LEV; % Verbosity level
        std_out = 0;                  % Define the standard output of the
    end

    methods (Access = private)
        % Guard the constructor against external invocation.  We only want
        % to allow a single instance of this class.  See description in
        % Singleton superclass.
        function this = Logger()
            % Initialisation of the variables
            if isdeployed
                this.color_mode = false;
            end
        end
    end

    methods (Static)
        function this = getInstance()
            % Concrete implementation.  See Singleton superclass.
            persistent unique_instance_logger__
            if isempty(unique_instance_logger__)
                this = Logger();
                unique_instance_logger__ = this;
            else
                this = unique_instance_logger__;
            end
        end
    end

    % =========================================================================
    %  MANAGING LOGGING
    % =========================================================================
    methods
        % Set read status mode --------------------------------------------
        function setColorMode(this, bool)
            % Set useage of colors in text output
            this.color_mode = bool;
        end

        function bool = getColorMode(this)
            % Get useage of colors in text output
            bool = this.color_mode;
        end

        % Set verbosity level ---------------------------------------------
        function setVerbosityLev(this, verbosity)
            % Set level of verbosity
            this.verbosity = verbosity;
        end

        function verbosity = getVerbosityLev(this)
            % Get level of verbosity
            verbosity = this.verbosity;
        end
    end

    % =========================================================================
    %  OUTPUT UTILITIES (respect verbosity)
    % =========================================================================
    methods
        function newLine(this, verbosity_level)
            % Send a message through the standard interface
            if (nargin < 2)
                verbosity_level = this.DEFAULT_VERBOSITY_LEV;
            end
            if (verbosity_level <= this.verbosity)
                fprintf('\n');
            end
        end

        function addMarkedMessage(this, text, verbosity_level)
            % Send a message through the standard interface
            if (nargin < 3)
                verbosity_level = this.DEFAULT_VERBOSITY_LEV;
            end
            text = strrep(text, char(10), char([10, 32 * ones(1,7)]));
            text = strrep(text, '\n', char([10, 32 * ones(1,7)]));
            if (verbosity_level <= this.verbosity)
                if this.color_mode
                    text = strrep(text, '\', '\\');
                    cprintf('Green','   **  ');
                    cprintf('text', strcat(text, '\n'));
                else
                    fprintf('   **  %s\n', text);
                end
            end
        end

        function addMessage(this, text, verbosity_level)
            % Send a message through the standard interface
            if (nargin < 3)
                verbosity_level = this.DEFAULT_VERBOSITY_LEV;
            end
            if (verbosity_level <= this.verbosity)
                text = strrep(text, char(10), char([10, 32]));
                text = strrep(text, '\n', char([10, 32]));
                fprintf(' %s\n', text);
            end
        end

        function addStatusOk(this, text, verbosity_level)
            % Send a message through the standard interface
            if (nargin < 3)
                verbosity_level = this.DEFAULT_VERBOSITY_LEV;
            end
            if (nargin < 2)
                text = '';
            end
            if (verbosity_level <= this.verbosity)
                this.printStatusOk(text);
            end
        end

        function addWarning(this, text, verbosity_level)
            % Send a warning through the standard interface
            if (nargin < 3)
                verbosity_level = this.WARNING_VERBOSITY_LEV;
            end
            if (verbosity_level <= this.verbosity)
                this.printWarning(text);
            end
        end

        function addError(this, text)
            % Send an error through the standard interface
            if (this.ERROR_VERBOSITY_LEV <= this.verbosity)
                this.printError(text);
            end
        end

        function status(this, status)
            % Display a flag of operation status ( -1 Err, 0 Ok, 1 Warning) on the previous line
            fprintf(char(08));
            this.opStatus(status, this.color_mode);
        end
    end


    % =========================================================================
    %    PRIVATE DISPLAY UTILITIES
    % =========================================================================

    methods (Access = 'protected')

        function printStatusOk(this, text, color_mode)
            % Display Warnings
            if (nargin == 2)
                color_mode = this.color_mode;
            end

            if isempty(text)
                fprintf('\b');
            else
                text = strrep(text, char(10), char([10, 32]));
                text = strrep(text, '\n', char([10, 32]));
                text = strrep(text, '\', '\\');
            end

            this.opStatus(0, color_mode);

            if (color_mode)
                cprintf('text', [text '\n']);
            else
                fprintf([text '\n']);
            end
        end

        function printWarning(this, text, color_mode)
            % Display Warnings
            if (nargin == 2)
                color_mode = this.color_mode;
            end
            this.opStatus(1, color_mode);
            text = strrep(text, char(10), char([10, 32 * ones(1,17)]));
            text = strrep(text, '\n', char([10, 32 * ones(1,17)]));
            text = strrep(text, '\', '\\');
            if (color_mode)
                cprintf([1 0.65 0], 'Warning: ');
                cprintf('text', [text '\n']);
            else
                fprintf(['WARNING: ' text '\n']);
            end
        end

        function printError(this, text, color_mode)
            % Display Errors
            if (nargin == 2)
                color_mode = this.color_mode;
            end
            this.opStatus(-1, color_mode);
            text = strrep(text, char(10), char([10, 32 * ones(1,15)]));
            text = strrep(text, '\n', char([10, 32 * ones(1,15)]));
            text = strrep(text, '\', '\\');
            if (color_mode)
                cprintf('err', 'Error: ');
                cprintf('text', [text '\n']);
            else
                fprintf(['ERROR: ' text '\n']);
            end
        end

    end

    % =========================================================================
    %    DISPLAY UTILITIES
    % =========================================================================

    methods(Static, Access = 'public')
        function opStatus(status, color_mode)
            % Display a flag of operation status ( -1 Err, 0 Ok, 1 Warning)
            if (nargin == 1)
                color_mode = true;
            end
            if color_mode
                cprintf('blue',' [ ');
                switch (status)
                    case 0, cprintf('Green','ok');
                    case 1, cprintf([1 0.65 0],'WW');
                    otherwise, cprintf('Red','!!');
                end
                cprintf('blue',' ] ');
            else
                fprintf(' [ ');
                switch (status)
                    case 0, fprintf('ok');
                    case 1, fprintf('WW');
                    otherwise, fprintf('!!');
                end
                fprintf(' ] ');
            end
        end
        
        function str = indent(str, n_spaces)
            str = strrep([char(ones(1,n_spaces) * 32) str], 10, char([10 ones(1, n_spaces) * 32]));
        end
    end

end
