%   SINGLETON CLASS Logger
% =========================================================================
%
% DESCRIPTION
%   Manage messages output depending on the current interface 
%   (up to now, only Command Windows is implemented)
%
% EXAMPLE
%   logger = Logger.getInstance();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Logger

%----------------------------------------------------------------------------------------------
%                           goGPS v0.9.1
% Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
% Written by:       Gatti Andrea
% Contributors:     Gatti Andrea, ...
%----------------------------------------------------------------------------------------------
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
%----------------------------------------------------------------------------------------------
classdef Logger < handle
    
    properties (Constant, GetAccess = private)
        DEFAULT_VERBOSITY_LEV   = 9;  % Default verbosity level
        WARNING_VERBOSITY_LEV   = 3;  % Warning verbosity level
        ERROR_VERBOSITY_LEV     = 1;  % Error verbosity level
    end

    properties (GetAccess = 'private', SetAccess = 'protected')
        color_mode = false;            % Flag for coloured output messages (if true requires cprintf)        
        verbosity = 50;                % Verbosity level 
    end
        
    methods (Access = private)
        % Guard the constructor against external invocation.  We only want
        % to allow a single instance of this class.  See description in
        % Singleton superclass.
        function obj = Logger()
            % Initialisation of the variables
        end
    end
    
    methods (Static)        
        function obj = getInstance()
            % Concrete implementation.  See Singleton superclass.
            persistent unique_instance_logger__
            if isempty(unique_instance_logger__)
                obj = Logger();
                unique_instance_logger__ = obj;
            else
                obj = unique_instance_logger__;
            end
        end
    end
    
    % =========================================================================
    %  MANAGING LOGGING
    % =========================================================================
    methods    
        % Set read status mode --------------------------------------------
        function setColorMode(obj, bool)
            % Set useage of colors in text output
            obj.color_mode = bool;
        end
        
        function bool = getColorMode(obj)
            % Get useage of colors in text output
            bool = obj.color_mode;
        end
        
        % Set verbosity level ---------------------------------------------
        function setVerbosityLev(obj, verbosity)
            % Set level of verbosity
            obj.verbosity = verbosity;
        end
        
        function verbosity = getVerbosityLev(obj)
            % Get level of verbosity
            verbosity = obj.verbosity;
        end
    end

    % =========================================================================
    %  OUTPUT UTILITIES (respect verbosity)
    % =========================================================================
    methods       
        function addMessage(obj, text, verbosity_level)
            % Send a message through the standard interface
            if (nargin < 3)
                verbosity_level = obj.DEFAULT_VERBOSITY_LEV;
            end
            if (verbosity_level <= obj.verbosity)
                fprintf('%s\n', text);
            end
        end
        
        function addStatusOk(obj, text, verbosity_level)
            % Send a message through the standard interface
            if (nargin < 3)
                verbosity_level = obj.DEFAULT_VERBOSITY_LEV;
            end
            if (verbosity_level <= obj.verbosity)
                obj.printStatusOk(text);
            end
        end
        
        function addWarning(obj, text)
            % Send a warning through the standard interface
            if (obj.WARNING_VERBOSITY_LEV <= obj.verbosity)
                obj.printWarning(text);
            end
        end
        
        function addError(obj, text)
            % Send a warning through the standard interface
            if (obj.ERROR_VERBOSITY_LEV <= obj.verbosity)
                obj.printError(text);
            end
        end
    end
    
    % =========================================================================
    %    PRIVATE DISPLAY UTILITIES
    % =========================================================================

    methods (Access = 'protected')

        function printStatusOk(obj, text, color_mode)
            % Display Warnings
            if (nargin == 2)
                color_mode = obj.color_mode;
            end
            obj.opStatus(0, color_mode);
            if (color_mode)
                cprintf('text', [text '\n']);
            else
                fprintf('%s\n', text);
            end
        end
        
        function printWarning(obj, text, color_mode)
            % Display Warnings
            if (nargin == 2)
                color_mode = obj.color_mode;
            end
            obj.opStatus(1, color_mode);
            if (color_mode)
                cprintf('SystemCommands', 'Warning: ');
                cprintf('text', [text '\n']);
            else
                fprintf('Warning: %s\n', text);
            end
        end
        
        function printError(obj, text, color_mode)
            % Display Errors
            if (nargin == 2)
                color_mode = obj.color_mode;
            end
            obj.opStatus(-1, color_mode);
            if (color_mode)
                cprintf('err', 'Error: ');
                cprintf('text', [text '\n']);
            else
                fprintf('Error: %s\n', text);
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
                    case 1, cprintf('SystemCommands','WW');
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
    end
    
end
