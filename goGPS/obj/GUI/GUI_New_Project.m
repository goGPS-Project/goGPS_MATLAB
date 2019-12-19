%   CLASS GUI_New_Project
% =========================================================================
%
% DESCRIPTION
%   class to manages the new project window of goGPS
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
%    |___/                    v 1.0 beta 5 Merry Christmas
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

classdef GUI_New_Project < handle
    
    %% PROPERTIES SINGLETON POINTERS
    % ==================================================================================================================================================
    properties % Utility Pointers to Singletons
        log
        state
    end
    
    %% PROPERTIES GUI
    % ==================================================================================================================================================
    properties
        w_main      % Handle of the main window 
        win         % Handle to this window
        
        dir_base    % handle to dir EditBox
        prj_name    % handle to prj EditBox
        prj_type    % handle to prj type List
        
        check_boxes % List of chgoGPS
        pop_ups     % List of drop down menu
        edit_texts  % List of editable text
        edit_texts_array % list of editable text array
        ceckboxes
    end    
    
    %% PROPERTIES STATUS
    % ==================================================================================================================================================
    properties (GetAccess = private, SetAccess = private)
        ok_go = false;
    end
    
    %% METHOD CREATOR
    % ==================================================================================================================================================
    methods (Static)
        function this = GUI_New_Project(w_main)
            % GUI_MAIN object creator
            this.init();
            this.openGUI();
            this.w_main = w_main;
        end
    end    
    %% METHODS INIT
    % ==================================================================================================================================================
    methods
        function init(this)
            this.log = Core.getLogger();
            this.state = Core.getState();
        end
        
        function status_ok = openGUI(this)
            % WIN CONFIGURATION
            % L| N|    W
            %
            %
            % ----------
            % b      b b
            %
            if ~isempty(this.w_main) && isvalid(this.w_main)
                close(this.w_main);
            end
            
            status_ok = true;
            this.ok_go = false;
            
            % empty check boxes
            this.check_boxes = {};
            this.pop_ups = {};
            this.edit_texts = {};
            this.edit_texts_array = {};
            
            % Main Window ----------------------------------------------------------------------------------------------
            
            win = figure( 'Name', 'Create an empty project', ...
                'Visible', 'off', ...
                'MenuBar', 'none', ...
                'ToolBar', 'none', ...
                'NumberTitle', 'off', ...
                'Position', [0 0 750 200]);
            
            this.win = win;
            
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
                this.log.addError('Please install GUI Layout Toolbox (https://it.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox)');
                open('GUI Layout Toolbox 2.3.1.mltbx');
                this.log.newLine();
                this.log.addWarning('After installation re-run goGPS');
                close(win);
                status_ok = false
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
            
            fnp = File_Name_Processor();
            % ProjectType
            [~, this.prj_type] = Core_UI.insertPopUpDark(panel_g_border, 'Project type', {'PPP', 'NET (short baselines) no iono - no tropo', 'NET  (medium baselines) no iono', 'NET (long baselines) iono free'}, 'prj_type', @this.none, [110 300]);
            
            % Folder
            [~, dir_base] = Core_UI.insertDirBoxDark(panel_g_border, 'Where to create', 'prj_home', @this.none, [120 -1 25]);                      
            dir_base.String = [fnp.getFullDirPath((fullfile(this.state.getHomeDir, '..'))) filesep];
            this.dir_base = dir_base;
            
            % Project Name
            [~, prj_name] = Core_UI.insertEditBoxDark(panel_g_border, 'Project name', 'prj_name', '', @this.none, [120 -1]);
            prj_name.HorizontalAlignment = 'left';
            prj_name.FontWeight = 'bold';
            prj_name.String = 'New_Project';
            this.prj_name = prj_name;
            
            panel_g_border.Heights = [30 30 30];
            
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
                'Spacing', 5, ...
                'HorizontalAlignment', 'right', ...
                'BackgroundColor', 0.14 * [1 1 1]);
            
            exit_but = uicontrol( 'Parent', bottom_bhl, ...
                'String', 'Cancel', ...
                'Callback', @this.close); %#ok<NASGU> 
            
            ok_but = uicontrol( 'Parent', bottom_bhr, ...
                'String', 'ok', ...
                'FontAngle', 'italic', ...
                'Callback', @this.okCreate, ...
                'FontWeight', 'bold'); %#ok<NASGU>
            
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
        
        function okCreate(this, caller, event)
            try
                Core_Utils.createEmptyProject(this.dir_base.String, this.prj_name.String, this.prj_type.Value);
                this.init();
                close(this.win);
            catch ex
                error(ex.message);
            end
            
            % Update main goGPS settings interface
            try
                this.w_main.init();
                this.w_main.updateUI(); 
            catch
                % Windows have been closed?
            end            
        end
    end
end
