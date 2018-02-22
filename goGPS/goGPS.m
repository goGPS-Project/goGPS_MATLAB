%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.6.0 alpha 1 - nightly
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

% Init output interfaces (singletons)


%if mode_user == 1
%    w_bar.setOutputType(1); % 0 means text, 1 means GUI, 5 both
%else
%w_bar.setOutputType(0); % 0 means text, 1 means GUI, 5 both
%end

%----------------------------------------------------------------------------------------------
% INTERFACE STARTUP
%----------------------------------------------------------------------------------------------

% Set global variable for goGPS obj mode
clearvars -global goObj;

if (mode_user == 1)
    
    % Now there's a unique interface for goGPS
    % to be compatible among various OSs the property "unit" of all the
    % elements must be set to "pixels"
    % (default unit is "character", but the size of a character is OS dependent)
    
    ok_go = core.openGUI();
    if (~ok_go)
        return
    end        
end

%% GO goGPS - here the computations start

core.prepareProcessing();

[state, log, w_bar] = core.getUtilities();

% start evaluating computation time
tic;

state.showTextMode();

sky = Core_Sky.getInstance();
for s = 1 : state.getSessionCount()
    %-------------------------------------------------------------------------------------------
    % SESSION START
    %-------------------------------------------------------------------------------------------
    
    fprintf('\n--------------------------------------------------------------------------\n');
    log.addMessage(sprintf('Starting session %d of %d', s, state.getSessionCount()));
    fprintf('--------------------------------------------------------------------------\n');
    
    % Init sky
    fr = File_Rinex(state.getTrgPath(1, s), 100);
    cur_date_start = fr.first_epoch.last();
    cur_date_stop = fr.last_epoch.first();
    sky.initSession(cur_date_start, cur_date_stop);
        
    clear rec;  % handle to all the receivers
    clear mst;
    r = 0;
    for i = 1 : state.getMstCount()
        log.newLine();
        log.addMessage(sprintf('Reading master %d of %d', i, state.getMstCount()));
        fprintf('--------------------------------------------------------------------------\n\n');
        
        r = r + 1;
        mst(i) = Receiver(state.getConstellationCollector(), state.getMstPath(i, s)); %#ok<SAGROW>
        mst(i).preProcessing();
        rec(r) = mst(i); %#ok<SAGROW>
    end
    
    clear ref;
    for i = 1 : state.getRefCount()
        log.newLine();
        log.addMessage(sprintf('Reading reference %d of %d', i, state.getRefCount));
        fprintf('--------------------------------------------------------------------------\n\n');
        
        r = r + 1;
        ref(i) = Receiver(state.getConstellationCollector(), state.getRefPath(i, s)); %#ok<SAGROW>
        ref(i).preProcessing();
        rec(r) = ref(i);        
    end
    
    clear trg;
    for i = 1 : state.getTrgCount()
        log.newLine();
        log.addMessage(sprintf('Reading target %d of %d', i, state.getTrgCount));
        fprintf('--------------------------------------------------------------------------\n\n');
        
        r = r + 1;
        trg(i) = Receiver(state.getConstellationCollector(), state.getTrgPath(i, s)); %#ok<SAGROW>        
        trg(i).preProcessing();
        rec(r) = trg(i);        
    end
    
    fprintf('--------------------------------------------------------------------------\n');
    log.newLine();
    log.addMarkedMessage('Syncing times, computing reference time');
    [p_time, id_sync] = Receiver.getSyncTime(rec, state.obs_type, state.getProcessingRate());
    
    for i = 1 : state.getTrgCount()
        trg(i).staticPPP([], id_sync{i}(:,1));
        % Align data and recompute solution
        %dt_i0 = trg(i).dt; 
        %trg(i).applyDtRec(dt_i0);         
        %trg(i).staticPPP([], id_sync{i}(:,1));
        %trg(i).dt = trg(i).dt + dt_i0;
    end
    
    trg_list(:,s) = trg; %#ok<SAGROW>
end
    
