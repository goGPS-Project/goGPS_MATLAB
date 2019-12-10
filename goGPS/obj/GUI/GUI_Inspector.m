%   CLASS GUI_Inspector
% =========================================================================
%
% DESCRIPTION
%   class to manages the inspector window of goGPS
%
% EXAMPLE
%   ui = GUI_Inspector.getInstance();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Core_UI


%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 4 ION
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

classdef GUI_Inspector < handle
    
    properties (Constant, Access = 'protected')
        BG_COLOR = Core_UI.DARK_GREY_BG;
    end
    
    %% PROPERTIES GUI
    % ==================================================================================================================================================
    properties
        win         % Handle to this window
        rec_tbl     % Handle to the table with all the receivers
        j_cmd       % Handle to the j_cmd java component
        
        edit_texts = {} % Handle to all the edit boxes
        flag_list = {}  % Handle to all the FLAGS
        
        flag_auto_exec      % Handle to automatic exec flag
        flag_auto_export    % Handle to automatic exeport flag
    end    
    
    %% PROPERTIES STATUS
    % ==================================================================================================================================================
    properties (GetAccess = private, SetAccess = private)
    end
    
    %% METHOD CREATOR
    % ==================================================================================================================================================
    methods  (Static, Access = private)
        function this = GUI_Inspector()
            % GUI_MAIN object creator
            this.openGUI();
            this.init();
        end
    end         
    
    methods (Static, Access = public)
        function this = getInstance()
            % Get the persistent instance of the class
            persistent unique_instance_GUI_Inspector__
            
            if isempty(unique_instance_GUI_Inspector__) || ~ishandle(unique_instance_GUI_Inspector__.win)
                this = GUI_Inspector();
                unique_instance_GUI_Inspector__ = this;
            else
                this = unique_instance_GUI_Inspector__;
                this.init();
            end
        end
    end
    
    %% METHODS INIT
    % ==================================================================================================================================================
    methods
        function init(this)
            this.updateRecList();
            this.updateUI();
        end
        
        function openGUI(this)
            % Main Window ---------------------------------------------------------------------------------------------
            
            win = figure( 'Name', 'goGPS inspector', ...
                'Visible', 'off', ...
                'MenuBar', 'none', ...
                'ToolBar', 'none', ...
                'NumberTitle', 'off', ...
                'Position', [0 0 1040, 640], ...
                'Resize', 'on');
            
            this.win = win;
            
            if isunix && not(ismac())
                 % win.Position(1) = round((win.Parent.ScreenSize(3) - win.Position(3)));
                % win.Position(2) = round((win.Parent.ScreenSize(4) - win.Position(4)));
                % centered
                win.Position(1) = round((win.Parent.ScreenSize(3) - win.Position(3)) / 2);
                win.Position(2) = round((win.Parent.ScreenSize(4) - win.Position(4)) / 2);
            else
                % top right
                % win.OuterPosition(1) = round((win.Parent.ScreenSize(3) - win.OuterPosition(3)));
                % win.OuterPosition(2) = round((win.Parent.ScreenSize(4) - win.OuterPosition(4)));
                % centered
                win.OuterPosition(1) = round((win.Parent.ScreenSize(3) - win.OuterPosition(3)) / 2);
                win.OuterPosition(2) = round((win.Parent.ScreenSize(4) - win.OuterPosition(4)) / 2);
            end
            
            try
                main_vb = uix.VBox('Parent', win, ...
                    'Padding', 5, ...
                    'BackgroundColor', Core_UI.DARKER_GREY_BG);                
            catch
                log = Core.getLogger;
                log.addError('Please install GUI Layout Toolbox (https://it.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox)');
                open('GUI Layout Toolbox 2.3.1.mltbx');
                log.newLine();
                log.addWarning('After installation re-run goGPS');
                close(win);
                return;
            end
            top_bh = uix.HBox('Parent', main_vb);
            
            logo_GUI_Inspector.BG_COLOR = Core_UI.DARK_GREY_BG;
            left_tbv = uix.VBox('Parent', top_bh, ...
                'BackgroundColor', logo_GUI_Inspector.BG_COLOR, ...
                'Padding', 5);
            
            % Logo/title box ------------------------------------------------------------------------------------------
            
            logo_g = uix.Grid('Parent', left_tbv, ...
                'Padding', 5, ...
                'BackgroundColor', logo_GUI_Inspector.BG_COLOR);
            
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
                        
            Core_UI.insertEmpty(left_tbv, logo_GUI_Inspector.BG_COLOR);
            left_tbv.Heights = [82 -1];
            
            % Title Panel ---------------------------------------------------------------------------------------------
            right_tvb = uix.VBox('Parent', top_bh, ...
                'Padding', 5, ...
                'BackgroundColor', logo_GUI_Inspector.BG_COLOR);

            top_bh.Widths = [106 -1];
            
            title = uix.HBox('Parent', right_tvb, ...
                'BackgroundColor', logo_GUI_Inspector.BG_COLOR);
            
            txt = this.insertBoldText(title, 'goGPS', 12, Core_UI.LBLUE, 'left');
            txt.BackgroundColor = logo_GUI_Inspector.BG_COLOR;
            title_l = uix.VBox('Parent', title, 'BackgroundColor', GUI_Inspector.BG_COLOR);
            title.Widths = [60 -1];
            Core_UI.insertEmpty(title_l, logo_GUI_Inspector.BG_COLOR)
            txt = this.insertBoldText(title_l, ['- software V' Core.GO_GPS_VERSION], 9, [], 'left');
            txt.BackgroundColor = logo_GUI_Inspector.BG_COLOR;
            title_l.Heights = [2, -1];
            
            % Top Panel -----------------------------------------------------------------------------------------------
            
            Core_UI.insertEmpty(right_tvb, logo_GUI_Inspector.BG_COLOR)
            txt = this.insertText(right_tvb, {'A GNSS processing software powered by GReD'}, 9, [], 'left');
            txt.BackgroundColor = logo_GUI_Inspector.BG_COLOR;
            right_tvb.Heights = [25 3 -1];
            
            Core_UI.insertEmpty(main_vb, Core_UI.DARKER_GREY_BG)
            
            % Main Panel ----------------------------------------------------------------------------------------------
            
            main_hb = uix.HBox('Parent', main_vb, ...
                'Padding', 5, ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);

            % Bottom Panel --------------------------------------------------------------------------------------------
                        
            Core_UI.insertEmpty(main_vb, Core_UI.DARKER_GREY_BG)
            bottom = uix.HBox('Parent', main_vb, ...
                'Padding', 5, ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
            [grp, this.edit_texts{end + 1}, this.flag_list{end + 1}] = Core_UI.insertDirBox(bottom, 'Out directory', 'out_dir', @this.onEditChange, [25 100 -1 25]);
            grp.BackgroundColor = Core_UI.DARK_GREY_BG;
            grp.Children(3).BackgroundColor = Core_UI.DARK_GREY_BG;
            grp.Children(3).ForegroundColor = [1 1 1];

            list_but = uix.HButtonBox( 'Parent', bottom, ...
                'ButtonSize', [100 28] , ...
                'VerticalAlignment', 'top', ...
                'HorizontalAlignment', 'right', ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
            
            load_core_but = uicontrol( 'Parent', list_but, ...
                'String', 'Load Core', ...
                'Callback', @this.onLoadCore); %#ok<NASGU>
            
            bottom.Widths = [-1 105];
            main_vb.Heights = [84 5 -1 5 35];
            
            % Middle Tab Panel ----------------------------------------------------------------------------------------
            
            this.insertRecList(main_hb);
            Core_UI.insertEmpty(main_hb, Core_UI.DARK_GREY_BG)
            tab_pnl_left = uix.TabPanel('Parent', main_hb, ...
                'TabWidth', 90, ...
                'Padding', 5, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG, ...
                'SelectionChangedFcn', @this.onTabChange);            
            
            this.insertTabSimplePlots(tab_pnl_left);
            this.insertTabMaps(tab_pnl_left);
            
            tab_pnl_left.TabTitles = {'Simple Plots', 'Maps'};

            % Right Panel ---------------------------------------------------------------------------------------------
            
            Core_UI.insertEmpty(main_hb, Core_UI.DARK_GREY_BG);
            
            tab_pnl_right = uix.TabPanel('Parent', main_hb, ...
                'TabWidth', 90, ...
                'Padding', 5, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG, ...
                'SelectionChangedFcn', @this.onTabChange);
                        
            this.j_cmd = this.insertCommandEditor(tab_pnl_right);
            
            tab_pnl_right.TabTitles = {'Commands'};
            
            main_hb.Widths = [190 5 360 5 -1];
            
            % Manage dimension ----------------------------------------------------------------------------------------
            
                        
            this.win.Visible = 'on';    
        end
        
        function close(this)
            if ~isempty(this.win) && ishandle(this.win)
                close(this.win);
            end
        end
    end
    %% METHODS INSERT
    % ==================================================================================================================================================
    methods (Static)
        function txt = insertBoldText(parent, title, font_size, color, alignment)
            if nargin < 4 || isempty(color)
                color = Core_UI.WHITE;
            end
            if nargin < 5 || isempty(alignment)
                alignment = 'center';
            end
            txt = uicontrol('Parent', parent, ...
                'Style', 'Text', ...
                'String', title, ...
                'ForegroundColor', color, ...
                'HorizontalAlignment', alignment, ...
                'FontSize', Core_UI.getFontSize(font_size), ...
                'FontWeight', 'bold', ...
                'BackgroundColor', GUI_Inspector.BG_COLOR);
        end

        function txt = insertText(parent, title, font_size, color, alignment)
            if nargin < 4 || isempty(color)
                color = Core_UI.WHITE;
            end
            if nargin < 5 || isempty(alignment)
                alignment = 'center';
            end
            txt = uicontrol('Parent', parent, ...
                'Style', 'Text', ...
                'String', title, ...
                'ForegroundColor', color, ...
                'HorizontalAlignment', alignment, ...
                'FontSize', Core_UI.getFontSize(font_size), ...
                'BackgroundColor', GUI_Inspector.BG_COLOR);
        end
    end
    
    methods (Access = 'private')
        function insertRecList(this, container)
            lv_box = uix.VBox('Parent', container, ...
                'Padding', 0, ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);            
            tv_text = uix.VBox( 'Parent', lv_box, ...
                'Padding', 0, ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
            
            list_title = uicontrol('Parent', tv_text, ...
                'Style', 'Text', ...
                'String', 'Receiver List', ...
                'ForegroundColor', Core_UI.WHITE, ...
                'HorizontalAlignment', 'left', ...
                'FontSize', Core_UI.getFontSize(9), ...
                'FontWeight', 'bold', ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
            
            list_but = uix.HButtonBox( 'Parent', tv_text, ...
                'Spacing', 5, ...
                'HorizontalAlignment', 'left', ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
                        
            select_but = uicontrol( 'Parent', list_but, ...
                'String', 'Select All', ...
                'Callback', @this.onSelectAll); %#ok<NASGU>
            clear_but = uicontrol( 'Parent', list_but, ...
                'String', 'Clear All', ...
                'Callback', @this.onUnselectAll); %#ok<NASGU>           
                        
            tv_text.Heights = [20 25];                        
            
            rec_g = uix.Grid('Parent', lv_box, ...
                'Padding', 0, ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
            
            this.rec_tbl = uitable('Parent', rec_g);
            this.rec_tbl.RowName = {}; 
            this.rec_tbl.ColumnName = {''; 'N'; 'Name'; 'SYS'};
            colTypes = {'logical', 'char', 'char', 'char'};
            this.rec_tbl.ColumnFormat = colTypes;
            this.rec_tbl.ColumnEditable = [true false false false];
            this.rec_tbl.ColumnWidth = {20, 45, 60, 46};        
            lv_box.Heights = [50 -1]; 
        end

        function insertTabSimplePlots(this, container)
            cmd_bg = Core_UI.LIGHT_GREY_BG;
            tab = uix.HBox('Parent', container, ...
                'Padding', 0, ...
                'BackgroundColor', cmd_bg);
             
            v_left = uix.VBox('Parent', tab, ...
                'Padding', 0, ...
                'BackgroundColor', cmd_bg);
            
            % --------------------------------------------------------
            % COMMAND SUGGESTIOS
            % --------------------------------------------------------
            eg_box = uix.VBox('Parent', v_left);
                        
            % uicontrol('Parent', eg_box, ...
            %     'Style', 'Text', ...
            %     'String', 'Execution examples:', ...
            %     'ForegroundColor', Core_UI.BLACK, ...
            %     'HorizontalAlignment', 'left', ...
            %     'FontSize', Core_UI.getFontSize(9), ...
            %     'BackgroundColor', cmd_bg);            
                        
            Core_UI.insertEmpty(eg_box);
            
            % Work-Space
            % --------------------------------------------------------
            last_sss_pnl = Core_UI.insertPanelLight(eg_box, 'Plots of the last computed session (Work-Space)');
            workspace_box = uix.VButtonBox('Parent', last_sss_pnl, ...
                'ButtonSize', [340 28] , ...
                'Spacing', 0, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'BackgroundColor', cmd_bg);
            
            but_line1 = uix.HButtonBox('Parent', workspace_box, ...
                'ButtonSize', [165 28] , ...
                'Spacing', 5, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'BackgroundColor', cmd_bg);
            
            uicontrol( 'Parent', but_line1, ...
                'String', 'Data Availability', ...
                'UserData', {'SHOW T@ DA'}, ...
                'Callback', @this.onInsertCommand);
            
            uicontrol( 'Parent', but_line1, ...
                'String', 'Observation Stats', ...
                'UserData', {'SHOW T@ OBS_STAT'}, ...
                'Callback', @this.onInsertCommand);
            
            but_line2 = uix.HButtonBox('Parent', workspace_box, ...
                'ButtonSize', [165 28] , ...
                'Spacing', 5, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'BackgroundColor', cmd_bg);

            uicontrol( 'Parent', but_line2, ...
                'String', 'SNR (polar)', ...
                'UserData', {'SHOW T@ SNR'}, ...
                'Callback', @this.onInsertCommand);
            
            uicontrol( 'Parent', but_line2, ...
                'String', 'SNR (polar - interpolated)', ...
                'UserData', {'SHOW T@ SNRI'}, ...
                'Callback', @this.onInsertCommand);
            
            but_line3 = uix.HButtonBox('Parent', workspace_box, ...
                'ButtonSize', [165 28] , ...
                'Spacing', 5, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'BackgroundColor', cmd_bg);
            
            uicontrol( 'Parent', but_line3, ...
                'String', 'Clock Errors', ...
                'UserData', {'SHOW T@ CKW'}, ...
                'Callback', @this.onInsertCommand);
           
            but_line4 = uix.HButtonBox('Parent', workspace_box, ...
                'ButtonSize', [165 28] , ...
                'Spacing', 5, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'BackgroundColor', cmd_bg);
           
            uicontrol( 'Parent', but_line4, ...
                'String', 'Outliers and CS', ...
                'UserData', {'SHOW T@ OCS'}, ...
                'Callback', @this.onInsertCommand);
            
            uicontrol( 'Parent', but_line4, ...
                'String', 'Outliers and CS (polar)', ...
                'UserData', {'SHOW T@ OCSP'}, ...
                'Callback', @this.onInsertCommand);
                        
            % --------------------------------------------------------
            
            Core_UI.insertEmpty(eg_box);
            
            % Output
            % --------------------------------------------------------
            
            out_pnl = Core_UI.insertPanelLight(eg_box, 'Plots of the stored outputs');
            out_box = uix.VButtonBox('Parent', out_pnl, ...
                'ButtonSize', [340 28] , ...
                'Spacing', 0, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'BackgroundColor', cmd_bg);
            
            but_line = uix.HButtonBox('Parent', out_box, ...
                'ButtonSize', [165 28] , ...
                'Spacing', 5, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'BackgroundColor', cmd_bg);
            
            uicontrol( 'Parent', but_line, ...
                'String', 'Clock Errors', ...
                'UserData', {'SHOW T@ CK'}, ...
                'Callback', @this.onInsertCommand);
                        
            but_line = uix.HButtonBox('Parent', out_box, ...
                'ButtonSize', [165 28] , ...
                'Spacing', 5, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'BackgroundColor', cmd_bg);
            
            uicontrol( 'Parent', but_line, ...
                'String', 'ZTD', ...
                'UserData', {'SHOW T@ ZTD'}, ...
                'Callback', @this.onInsertCommand);
            
            uicontrol( 'Parent', but_line, ...
                'String', 'ZHD', ...
                'UserData', {'SHOW T@ ZHD'}, ...
                'Callback', @this.onInsertCommand);
                        
            uicontrol( 'Parent', but_line, ...
                'String', 'ZWD', ...
                'UserData', {'SHOW T@ ZWD'}, ...
                'Callback', @this.onInsertCommand);
            
            uicontrol( 'Parent', but_line, ...
                'String', 'PWV', ...
                'UserData', {'SHOW T@ PWV'}, ...
                'Callback', @this.onInsertCommand);

            but_line = uix.HButtonBox('Parent', out_box, ...
                'ButtonSize', [165 28] , ...
                'Spacing', 5, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'BackgroundColor', cmd_bg);
            
            uicontrol( 'Parent', but_line, ...
                'String', 'ZTD vs Height', ...
                'UserData', {'SHOW T@ ZTD_VSH'}, ...
                'Callback', @this.onInsertCommand);
            
            uicontrol( 'Parent', but_line, ...
                'String', 'ZWD vs Height', ...
                'UserData', {'SHOW T@ ZWD_VSH'}, ...
                'Callback', @this.onInsertCommand);

            but_line = uix.HButtonBox('Parent', out_box, ...
                'ButtonSize', [165 28] , ...
                'Spacing', 5, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'BackgroundColor', cmd_bg);
            
            uicontrol( 'Parent', but_line, ...
                'String', 'Pressure / Temp. / Humidity', ...
                'UserData', {'SHOW T@ PTH'}, ...
                'Callback', @this.onInsertCommand);
            
            % --------------------------------------------------------
            
            %scroller = uix.ScrollingPanel('Parent', eg_box);
            %container = uix.Grid('Parent', scroller, ...
            %    'BackgroundColor', Core_UI.LIGHT_GREY_BG);

            eg_box.Heights = [5, 140, 5, -1];
        end
        
        function insertTabMaps(this, container)
            cmd_bg = Core_UI.LIGHT_GREY_BG;
            tab = uix.HBox('Parent', container, ...
                'Padding', 0, ...
                'BackgroundColor', cmd_bg);
             
            v_left = uix.VBox('Parent', tab, ...
                'Padding', 0, ...
                'BackgroundColor', cmd_bg);
            
            % --------------------------------------------------------
            % COMMAND SUGGESTIOS
            % --------------------------------------------------------
            eg_box = uix.VBox('Parent', v_left);
                        
            % uicontrol('Parent', eg_box, ...
            %     'Style', 'Text', ...
            %     'String', 'Execution examples:', ...
            %     'ForegroundColor', Core_UI.BLACK, ...
            %     'HorizontalAlignment', 'left', ...
            %     'FontSize', Core_UI.getFontSize(9), ...
            %     'BackgroundColor', cmd_bg);                        
                        
            Core_UI.insertEmpty(eg_box);            
            
            % Output
            % --------------------------------------------------------
            
            out_pnl = Core_UI.insertPanelLight(eg_box, 'Static Maps');
            out_box = uix.VButtonBox('Parent', out_pnl, ...
                'ButtonSize', [340 28] , ...
                'Spacing', 0, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'BackgroundColor', cmd_bg);
                        
            uicontrol( 'Parent', out_box, ...
                'String', 'Map of GNSS stations (satellite image background)', ...
                'UserData', {'SHOW T@ MAP'}, ...
                'Callback', @this.onInsertCommand);
            
            uicontrol( 'Parent', out_box, ...
                'String', 'Browsable map of GNSS stations (satellite image background)', ...
                'UserData', {'SHOW T@ L_MAP'}, ...
                'Callback', @this.onInsertCommand);
            
            uicontrol( 'Parent', out_box, ...
                'String', 'Map of GNSS stations (DTM background)', ...
                'UserData', {'SHOW T@ DTM_MAP'}, ...
                'Callback', @this.onInsertCommand);                        
            
            % --------------------------------------------------------
            
            %scroller = uix.ScrollingPanel('Parent', eg_box);
            %container = uix.Grid('Parent', scroller, ...
            %    'BackgroundColor', Core_UI.LIGHT_GREY_BG);

            eg_box.Heights = [5, -1];
        end
        
        function j_cmd = insertCommandEditor(this, container)
            % COMMAND LIST EDITOR
            % --------------------------------------------------------
                        
            cmd_bg = Core_UI.LIGHT_GREY_BG;
            
            cmd_box = uix.VBox('Parent', container, ...
                'Padding', 0, ...
                'BackgroundColor', cmd_bg);
            
            this.flag_auto_exec = uicontrol('Parent', cmd_box,...
                'Style', 'checkbox',...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG, ...
                'FontSize', Core_UI.getFontSize(9), ...
                'String', 'Immediate execution');
            
            
            this.flag_auto_export = uicontrol('Parent', cmd_box,...
                'Style', 'checkbox',...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG, ...
                'FontSize', Core_UI.getFontSize(9), ...
                'String', 'Add automatic plot export');
            
            Core_UI.insertEmpty(cmd_box);
            
            uicontrol('Parent', cmd_box, ...
                'Style', 'Text', ...
                'String', 'Insert here the goGPS command list:', ...
                'FontWeight', 'normal', ...
                'ForegroundColor', Core_UI.BLACK, ...
                'HorizontalAlignment', 'left', ...
                'FontSize', Core_UI.getFontSize(9), ...
                'BackgroundColor', cmd_bg);

            j_cmd = com.mathworks.widgets.SyntaxTextPane;
            codeType = j_cmd.M_MIME_TYPE;  % j_settings.contentType='text/m-MATLAB'
            j_cmd.setContentType(codeType);
            str = '';
            j_cmd.setText(str);
            % Create the ScrollPanel containing the widget
            j_scroll_settings = com.mathworks.mwswing.MJScrollPane(j_cmd);
            % Inject edit box with the Java Scroll Pane into the main_window
            [panel_j, panel_h] = javacomponent(j_scroll_settings, [1 1 1 1], cmd_box);
                    
            % HELP
            Core_UI.insertEmpty(cmd_box, cmd_bg);
            
            list_but = uix.HButtonBox( 'Parent', cmd_box, ...
                'Spacing', 5, ...
                'HorizontalAlignment', 'right', ...
                'ButtonSize', [120 28] , ...
                'BackgroundColor', cmd_bg);
                        
            help = uicontrol( 'Parent', list_but, ...
                'String', 'HELP', ...
                'Callback', @this.onOpenCommandHelp);
            clr = uicontrol( 'Parent', list_but, ...
                'String', 'Clear', ...
                'Callback', @this.onClearCommands);
            check = uicontrol( 'Parent', list_but, ...
                'String', 'Check validity', ...
                'Callback', @this.onCheckValidity);
            exec = uicontrol( 'Parent', list_but, ...
                'String', 'EXEC', ...
                'Callback', @this.exec);
            
            cmd_box.Heights = [Core_UI.LINE_HEIGHT, Core_UI.LINE_HEIGHT, 5, Core_UI.LINE_HEIGHT, -1, 2, 30];
        
            % --------------------------------------------------------
        end
    end
    %% METHODS setters
    % ==================================================================================================================================================
    methods
        function addMessage(this, text, type)
            % Add a message to the logger
            % 
            % INPUT
            %   text    text in HTML format
            %   type    'm'     marked message
            %           'w'     warning message
            %           'e'     error message
            %           otherwise normal
            %
            % SYNTAX
            %   this.addHTML(text, type)
            
            if nargin < 3 || isempty(type)
                type = 'n';
            end
            Core_UI.guiAddMessage(this.jedt, text, type);
        end
        
        function addHTML(this, html_txt)
            % Add a message to the logger
            % 
            % INPUT
            %   html_txt    text in HTML format
            %
            % SYNTAX
            %   this.addHTML(html_txt)

            
            Core_UI.guiAddHTML(this.jedt, html_txt);
        end
        
        function updateRecList(this)
            % Update rec table
            core = Core.getCurrentCore;
            this.rec_tbl.Data = cell(numel(core.rec), 4);
            tmp = cell(numel(core.rec), 4);
            for r = 1 : numel(core.rec)
                name = core.rec(r).getMarkerName4Ch;                                
                
                tmp{r,1} = true;
                tmp{r,2} = sprintf('%s style="font-weight: bold; font-size: 9px; color: #6666FF; ">%d', '<html><tr><td width=9999 align=center ', r);
                tmp{r,3} = sprintf('%s style="font-weight: bold; font-size: 9px; color: #6666FF; ">%s', '<html><tr><td width=9999 align=center ', upper(name(1:4)));
                tmp{r,4} = sprintf('%s style="font-weight: bold; font-size: 9px; color: #444444; ">%s', '<html><tr><td width=9999 align=center ', core.rec(r).getActiveSys());
            end
            this.rec_tbl.Data = tmp;
        end
        
        function updateUI(this)
            state = Core.getState();
            for i = 1 : length(this.edit_texts)
                value = state.getProperty(this.edit_texts{i}.UserData);
                if ~isempty(value)
                    this.edit_texts{i}.String = value;
                end
            end
            
            this.updateFlagList
        end
        
        function updateFlagList(this)
            Core_UI.checkFlag(this.flag_list)
        end               

        function cleaned_cmd_list = checkCommands(this)
            % Check Commands
            cleaned_cmd_list = {};
            if isempty(char(this.j_cmd.getText))
                this.j_cmd.setText('% Write here the commands to be executed');
            else
                cmd_list = textscan(strrep(char(this.j_cmd.getText),'%','#'),'%s','Delimiter', '\n');
                cmd = Core.getCommandInterpreter();
                if ~isempty(cmd_list)
                    [cleaned_cmd_list, err_list, execution_block, sss_list, trg_list, key_lev, flag_push, flag_parallel] = cmd.fastCheck(cmd_list{1});
                    key_lev = key_lev - (diff([0 key_lev]) > 0);
                    
                    cid = 0; % index running on valid commands
                    for c = 1 : numel(cmd_list{1})
                        cid = cid + ~err_list(c);
                        
                        cur_cmd = cmd_list{1}{c};
                        if (length(cur_cmd) > 1) && (cur_cmd(1) ~= '#') && err_list(c)
                            cur_cmd = ['# ' cur_cmd ' - ERROR: CMD UNKNOWN']; %#ok<AGROW>
                            cmd_list{1}{c} = cur_cmd;
                        elseif ~err_list(c)
                            cmd_list{1}{c} = sprintf('%s%s', char(32 * ones(1,3 * key_lev(cid))), strtrim(cleaned_cmd_list{cid}));
                        end
                    end
                    str = strrep(strCell2Str(cmd_list{1}, 10),'#','%');
                    this.j_cmd.setText(str);
                end
            end
            Core.getLogger.addMarkedMessage('The command validity has been checked');
        end

        function trg_list = getTargetList(this)
            % Get the target list from the Receiver list checked boxes
            %
            % SYNTAX
            %   trg_list = this.getTargetList()
            
            trg = [this.rec_tbl.Data{:,1}];
            trg_list = '';
            if all(trg)
                trg_list = 'T*';
            elseif any(trg)
                trg = sprintf('%d,', find(trg));
                trg_list = ['T', trg(1:end-1)];
            end
        end
    end
    
    %% METHODS EVENTS
    % ==================================================================================================================================================
    methods (Access = public)         
        function onSelectAll(this, caller, event)
            % Select all the receivers
            for r = 1 : size(this.rec_tbl.Data, 1)
                this.rec_tbl.Data{r,1} = true;
            end
        end
        
        function onLoadCore(this, caller, event)
            % Load state settings
            
            if ~exist('core', 'var')
                addPathGoGPS;
            end
            core_dir = Core.getState.getOutDir();
            
            [file_name, path_name] = uigetfile({'*.mat;','goGPS core from previous session (*.mat)';}, 'Choose file with saved core', core_dir);
            
            if path_name ~= 0 % if the user pressed cancelled, then we exit this callback
                % get the extension (mat/ini):
                [~, ~, ext] = fileparts(file_name);
                
                % build the path name of the file to be loaded
                core_file = fullfile(path_name, file_name);
                if strcmp(ext, '.mat')
                    Core.getLogger.addMarkedMessage(sprintf('Start loading core from "%s"', core_file));
                    load(core_file);
                    if ~exist('core', 'var')
                        Core.getLogger.addError(sprintf('No core variable found in file', core_file));
                    else
                        core.setCurrentCore(core);
                        this.init();
                        Core.getLogger.addStatusOk(sprintf('Core successifully loaded', core_file));
                    end
                else
                    Core.getLogger.addError('Unrecognized input file format!');
                end
            end
        end
        
        function onUnselectAll(this, caller, event)
            % Select all the receivers
            for r = 1 : size(this.rec_tbl.Data, 1)
                this.rec_tbl.Data{r,1} = false;
            end
        end
        
        function onTabChange(this, caller, event)
        end
        
        function onOpenCommandHelp(this, caller, event)
            % Open Help Window
            GUI_Command_Help;
        end

        function onCheckValidity(this, caller, event)
            % Check Validity of the command window
            this.checkCommands();
        end
        
        function onClearCommands(this, caller, event)
            % Clear the command window
            this.j_cmd.setText('');
            this.checkCommands();
        end
        
        function exec(this, caller, event)
            cleaned_cmd_list = this.checkCommands();
            if ~isempty(cleaned_cmd_list)
                Core.getCurrentCore.exec(cleaned_cmd_list);
            end
        end
                
        function onInsertCommand(this, caller, event)
            cmd_list = {};
            txt = char(this.j_cmd.getText);
            if ~isempty(txt)
                cmd_list = textscan(strrep(txt, '%', '#'), '%s', 'Delimiter', '\n');
                cmd_list = cmd_list{1};
            end
            new_cmd = strrep(caller.UserData(:), 'T@', this.getTargetList());
            for l = 1 : numel(new_cmd)
                if strcmp(new_cmd{l}(1:4), 'SHOW')
                    if this.flag_auto_export.Value
                        new_cmd{l} = [new_cmd{l} ' -e=".png"'];
                    end
                end
            end
            cmd_list = [cmd_list; new_cmd];
            str = strrep(strCell2Str(cmd_list, 10),'#','%');
            this.j_cmd.setText(str);
            this.checkCommands();
            
            % If immediate execution is required            
            if this.flag_auto_exec.Value
                core = Core.getCurrentCore();
                cmd = core.getCommandInterpreter();
                core.exec(cmd.fastCheck(new_cmd));
            end            
        end
        
        function onEditChange(this, caller, event)
            % Manage edit box change events
            state = Core.getState;
            prop = state.getProperty(caller.UserData);
            if ~isnumeric(prop)
                state.setProperty(caller.UserData, caller.String);
            else
                state.setProperty(caller.UserData, str2num(caller.String));
            end
            
            state.check();
            caller.String = state.getProperty(caller.UserData);            
            this.updateFlagList();
        end
    end
end
