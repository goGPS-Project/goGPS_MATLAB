%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.6.0 alpha 2 - nightly
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Originally written by:       Mirko Reguzzoni, Eugenio Realini
%  Contributors:                Gatti Andrea, Giulio Tagliaferro, ...
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

% clear all variables
% NOTE: using only 'clearvars' does not clear global variables, while using
% 'clear all' removes breakpoints
clearvars -global -except ini_settings_file use_gui; %  exceptions for goGPSgo
clearvars -except ini_settings_file use_gui; % exceptions for goGPSgo

% if the plotting gets slower than usual, there might be problems with the
% Java garbage collector. In case, you can try to use the following
% command:
%
% java.lang.System.gc() %clear the Java garbage collector
%
% or:
%
% clear java

% add all the subdirectories to the search path
if (~isdeployed)
    addpath(genpath(pwd));
end

%----------------------------------------------------------------------------------------------
% INIT
%----------------------------------------------------------------------------------------------

% Init Core
core = Core.getInstance();

if exist('ini_settings_file', 'var')
    core.importIniFile(ini_settings_file);
end

%----------------------------------------------------------------------------------------------
% INTERFACE TYPE DEFINITION
%----------------------------------------------------------------------------------------------

if exist('use_gui', 'var')
    mode_user = use_gui;
else
    mode_user =   1; % user interface type
    % mode_user = 0 --> use text interface
    % mode_user = 1 --> use GUI
end

%----------------------------------------------------------------------------------------------
% INTERFACE STARTUP
%----------------------------------------------------------------------------------------------

if (mode_user == 1)

    % Now there's a unique interface for goGPS
    % to be compatible among various OSs the property "unit" of all the
    % elements must be set to "pixels"
    % (default unit is "character", but the size of a character is OS dependent)
    ui = Core_UI.getInstance();
    ui.openGUI();

    if ~ui.isGo()
        return
    end
end

%% GO goGPS - here the computations start

core.prepareProcessing();
core.go();
