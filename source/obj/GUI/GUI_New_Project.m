%   CLASS GUI_New_Project
% =========================================================================
%
% DESCRIPTION
%   class to manage the new project window of goGPS
%
% EXAMPLE
%   ui = Core_UI.getInstance();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Core_UI


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
%  Contributors:      Andrea Gatti
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

classdef GUI_New_Project < GUI_Unique_Win
    properties (Constant)
        WIN_NAME = 'goGPS_New_Prj_Win';
    end
    
    %% PROPERTIES SINGLETON POINTERS
    % ==================================================================================================================================================
    properties % Utility Pointers to Singletons
        log
        state
    end
    
    %% PROPERTIES GUI
    % ==================================================================================================================================================
    properties
        w_edt      % Handle of the main window
        win         % Handle to this window
        
        dir_base    % handle to dir EditBox
        dir_rin     % source of observation RINEX files
        rin_op      % Type of operation on the Observation Folder
        prj_name    % handle to prj EditBox
        prj_type    % handle to prj type List
        
        check_boxes % List of chgoGPS
        pop_ups     % List of drop down menu
        edit_texts  % List of editable text
        edit_texts_array % list of editable text array
        ceckboxes
    end
    
    %% METHOD CREATOR
    % ==================================================================================================================================================
    methods (Static)
        function this = GUI_New_Project(w_edt)
            % GUI_MAIN object creator
            this.init();
            this.openGUI();
            this.w_edt = w_edt;
        end
    end
    %% METHODS INIT
    % ==================================================================================================================================================
    methods
        function init(this)
            this.log = Core.getLogger();
            this.state = Core.getState();
        end
        
        function openGUI(this)
            % Main Window ----------------------------------------------------------------------------------------------
            
            % If there is still an old logging wondow still open, close it
            old_win = this.getUniqueWinHandle();
            if ~isempty(old_win)
                delete(old_win);
            end
            
            win = figure( 'Name', 'Create an new project', ...
                'Visible', 'off', ...
                'DockControls', 'off', ...
                'MenuBar', 'none', ...
                'ToolBar', 'none', ...
                'NumberTitle', 'off', ...
                'Position', [0 0 750 250]);
            
            win.UserData.name = this.WIN_NAME;
            this.win = win;
            
            % empty check boxes
            this.check_boxes = {};
            this.pop_ups = {};
            this.edit_texts = {};
            this.edit_texts_array = {};
            
            % Main Window ----------------------------------------------------------------------------------------------
            
            if isunix && not(ismac())
                win.Position(1) = round((win.Parent.ScreenSize(3) - win.Position(3)) / 2);
                win.Position(2) = round((win.Parent.ScreenSize(4) - win.Position(4)) / 2);
            else
                win.OuterPosition(1) = round((win.Parent.ScreenSize(3) - win.OuterPosition(3)) / 2);
                win.OuterPosition(2) = round((win.Parent.ScreenSize(4) - win.OuterPosition(4)) / 2);
            end
            
            try
                main_bv = uix.VBox('Parent', win, ...
                    'Padding', 5, ...
                    'BackgroundColor', Core_UI.DARK_GREY_BG);
            catch
                log = Core.getLogger;
                log.setOutMode(1,[],0); % to plot a Warning I need to disable GUI and enable
                log.addError('Please install GUI Layout Toolbox (https://it.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox)');
                open('GUI Layout Toolbox 2.3.4.mltbx');
                log.newLine();
                log.addWarning('After installation re-run goGPS');
                close(win);
                return;
            end
            top_bh = uix.HBox( 'Parent', main_bv);
            
            left_bv = uix.VBox('Parent', top_bh, ...
                'Padding', 5, ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
            
            % Logo/title box -------------------------------------------------------------------------------------------
            
            logo_bg_color = Core_UI.DARK_GREY_BG;
            logo_g = uix.Grid('Parent', left_bv, ...
                'Padding', 5, ...
                'BackgroundColor', logo_bg_color);
            
            logo_ax = axes( 'Parent', logo_g);
            logo_g.Widths = 64;
            logo_g.Heights = 64;
            [logo, transparency] = Core_UI.getLogo();
            logo(repmat(sum(logo,3) == 0,1,1,3)) = 0;
            logo = logo - 20;
            image(logo_ax, ones(size(logo)), 'AlphaData', transparency);
            logo_ax.XTickLabel = [];
            logo_ax.YTickLabel = [];
            axis off;
            
            Core_UI.insertEmpty(left_bv, Core_UI.DARK_GREY_BG);
            
            % Main Panel -----------------------------------------------------------------------------------------------
            
            panel_g_border = uix.VBox('Parent', top_bh, ...
                'Padding', 5, ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
            %panel = uix.BoxPanel('Parent', panel_border, 'Title', 'Settings' );
            
            new_field = uix.HBox('Parent', panel_g_border, ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
            Core_UI.insertEmpty(new_field, Core_UI.DARK_GREY_BG);
            new_field.Widths = 25;
            txt = Core_UI.insertText(new_field, {'Create a new project folder containing an initial configuration file'}, 9, [], [], 'left');
            txt.FontWeight = 'bold';
            
            fnp = File_Name_Processor();
            % ProjectType
            new_field = uix.HBox('Parent', panel_g_border, ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
            Core_UI.insertEmpty(new_field, Core_UI.DARK_GREY_BG);
            new_field.Widths = 25;
            [~, this.prj_type] = Core_UI.insertPopUpDark(new_field, 'Project type', {'PPP - Precise Point Positioning', 'NET (short baselines) no iono - no tropo', 'NET (medium baselines) no iono', 'NET (long baselines) iono - free network'}, 'prj_type', @this.none, [143 300]);
            
            % Folder
            [~, dir_base] = Core_UI.insertDirBoxDark(panel_g_border, 'Where to create', 'prj_home', @this.none, [25 150 -1 25]);
            dir_base.String = fnp.getFullDirPath((fullfile(this.state.getHomeDir, '..')));
            this.dir_base = dir_base;
            
            % Project Name
            new_field = uix.HBox('Parent', panel_g_border, ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
            Core_UI.insertEmpty(new_field, Core_UI.DARK_GREY_BG);
            new_field.Widths = 25;
            [~, prj_name] = Core_UI.insertEditBoxDark(new_field, 'Project name', 'prj_name', '', @this.none, [150 285]);
            prj_name.HorizontalAlignment = 'left';
            prj_name.FontWeight = 'bold';
            prj_name.String = 'New_Project';
            this.prj_name = prj_name;
            
            Core_UI.insertEmpty(panel_g_border, Core_UI.DARK_GREY_BG);
            
            % Source of RINEX
            [~, dir_rin] = Core_UI.insertDirBoxDark(panel_g_border, 'Observations folder', 'prj_home', @this.none, [25 150 -1 25]);
            dir_path = fnp.getFullDirPath((fullfile(Core.getInstallDir, '../data/project/default_PPP/RINEX')));
            if exist(dir_path, 'dir') == 7
                dir_rin.String = dir_path;
            else
                dir_rin.String = fnp.getFullDirPath((fullfile(this.state.getHomeDir, '..')));
            end
            this.dir_rin = dir_rin;
            
            % ProjectType
            new_field = uix.HBox('Parent', panel_g_border, ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
            Core_UI.insertEmpty(new_field, Core_UI.DARK_GREY_BG);
            new_field.Widths = 25;
            [~, this.rin_op] = Core_UI.insertPopUpDark(new_field, '', {'Copy the observations folder into the new project', 'Move the observations folder into the new project', 'Keep observations in the current folder', 'Do not add any receiver now'}, 'rin_op', @this.onUpdateObsType, [143 300]);
            panel_g_border.Heights = [35 30 30 30 15 30 30];
            
            % Botton Panel ---------------------------------------------------------------------------------------------
            bottom_bh = uix.HBox( 'Parent', main_bv, ...
                'Padding', 5, ...
                'Spacing', 5, ...
                'BackgroundColor', 0.14 * [1 1 1]);
            
            bottom_bhl = uix.HButtonBox( 'Parent', bottom_bh, ...
                'Spacing', 5, ...
                'HorizontalAlignment', 'left', ...
                'BackgroundColor', 0.14 * [1 1 1]);
            
            bottom_bhr = uix.HButtonBox( 'Parent', bottom_bh, ...
                'ButtonSize', [120 28] , ...
                'Spacing', 5, ...
                'HorizontalAlignment', 'right', ...
                'BackgroundColor', 0.14 * [1 1 1]);
            
            exit_but = uicontrol( 'Parent', bottom_bhl, ...
                'String', 'Cancel', ...
                'Callback', @this.close); %#ok<NASGU>
            
            ok_but = uicontrol( 'Parent', bottom_bhr, ...
                'String', 'Create New Project', ...
                'Callback', @this.okCreate); %#ok<NASGU>
            
            % Manage dimension -------------------------------------------------------------------------------------------
            
            main_bv.Heights = [-1 30];
            %session_height = sum(left_bv.Children(2).Children(1).Heights);
            left_bv.Heights = [82 -1];
            top_bh.Widths = [82 -1];
            
            this.win.Visible = 'on';
            %uiwait(this.win);
        end
    end
    %% METHODS INSERT
    % ==================================================================================================================================================
    methods
        
    end
    %% METHODS getters
    % ==================================================================================================================================================
    methods
        function ok_go = isGo(this)
            ok_go = this.ok_go;
        end
    end
    
    %% METHODS EVENTS
    % ==================================================================================================================================================
    methods (Access = public)
        function none(this, caller, event)
        end
        
        function close(this, caller, event)
            close(this.win);
        end
        
        function onUpdateObsType(this, caller, event)
            if this.rin_op.Value == 4
                this.dir_rin.Enable = 'off'; % if no RINEX must be added the input directory must be disabled
            else
                this.dir_rin.Enable = 'on';
            end
        end
        
        function okCreate(this, caller, event)
            try
                log = Core.getLogger;
                log.addMarkedMessage(sprintf('Creating a new project: "%s"\ninto: "%s"', this.prj_name.String, this.dir_base.String));
                Core_Utils.createEmptyProject(this.dir_base.String, this.prj_name.String, this.prj_type.Value);
                this.init();
                state = Core.getCurrentSettings;
                obs_path = this.dir_rin.String;
                if this.rin_op.Value < 4 && ~exist(obs_path, 'dir') == 7
                    log.addError('Observations folder does not exist!!!');
                else
                    switch this.rin_op.Value
                        case 1 % copy
                            log.addMessage(log.indent(sprintf('Start copying "%s" to "%s"', [obs_path filesep '*'], state.getObsDir)));
                            [flag_ok, msg, msg_id] = copyfile([obs_path filesep '*'], state.getObsDir, 'f');
                            if flag_ok
                                log.addStatusOk('Data copied successfully!');
                            else
                                log.addError(msg);
                            end
                        case 2 % move
                            log.addMessage(log.indent(sprintf('Start moving "%s" to "%s"', [obs_path filesep '*'], state.getObsDir)));
                            [flag_ok, msg, msg_id] = movefile([obs_path filesep '*'], state.getObsDir, 'f');
                            if flag_ok
                                log.addStatusOk('Data moved successfully!');
                            else
                                log.addError(msg);
                            end
                        case 3 % keep the files there
                            state.setObsDir(obs_path);
                        otherwise % do nothing
                    end
                    
                    if this.rin_op.Value < 4                        
                        state.setObsName(Core_Utils.getStationList(state.getObsDir, 'oO', true));
                    end
                end
                        
                close(this.win);
            catch ex
                error(ex.message);
            end
            
            % Update main goGPS settings interface
            try
                this.w_edt.init();
                this.w_edt.updateUI();
                rf = Core.getReferenceFrame;
                rf.init(state.getCrdFile);
                this.updateCooTable();
            catch
                % Windows have been closed?
            end
        end
    end
end
