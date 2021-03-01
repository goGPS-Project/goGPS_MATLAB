%   CLASS GUI_Downloader
% =========================================================================
%
% DESCRIPTION
%   class to manage the user interface of goGPS
%
% EXAMPLE
%   ui = GUI_Downloader.getInstance();
%   ui.openGUI();
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
%  Copyright (C) 2021 Geomatics Research & Development srl (GRed)
%  Written by:        Andrea Gatti
%  Contributors:      Andrea Gatti, Giulio Tagliaferro
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

classdef GUI_Downloader < GUI_Unique_Win   
    properties (Constant)
        WIN_NAME = 'goGPS_Downloader';
    end
        
    %% PROPERTIES GUI
    % ==================================================================================================================================================
    properties
        win         % Handle of the main window
        menu        % Handle of the menu
        go_but      % Handle to goButton
        
        info_g      % Info group
        rec_tbl     % Receiver table
        session_panel % Panel of the session definition 
        session_info % Session info
        ui_sss_start
        ui_sss_stop
        
        ini_path    % ini path text box
        check_boxes % List of all the checkboxes
        check_boxes_dwn % List of all the download checkboxes
        pop_ups     % List of drop down menu
        rpop_up     % Remote resources pup-up
        ropref      % Remote Orbit Preferences
        ripref      % Remote Iono preferences
        rv2pref     % Remote VMF source preferences
        edit_texts  % List of editable text
        edit_texts_array % list of editable text array
        flag_list   % list of all the flags
        
        uip         % User Interface Pointers
    end    
    %% PROPERTIES STATUS
    % ==================================================================================================================================================
    properties (GetAccess = private, SetAccess = private)
        ok_go = false;
    end
    
    %% METHOD CREATOR
    % ==================================================================================================================================================
    methods (Static, Access = private)
        function this = GUI_Downloader(flag_wait)
            % GUI_MAIN object creator
            this.init();
            this.openGUI(flag_wait);
        end
    end    
    
    methods (Static, Access = public)
        function this = getInstance(flag_wait)
            if (nargin < 1) || isempty(flag_wait)
                flag_wait = false;
            end
                
            % Get the persistent instance of the class
            persistent unique_instance_gui_downloader__
            
            if isempty(unique_instance_gui_downloader__)
                this = GUI_Downloader(flag_wait);
                unique_instance_gui_downloader__ = this;
                if isvalid(this.win) && flag_wait
                    uiwait(this.win);
                end
            else
                this = unique_instance_gui_downloader__;
                this.init();
                this.openGUI(flag_wait);
                if isvalid(this.win) && flag_wait
                    uiwait(this.win);
                end
            end
        end
        
        function closeGUI()
            fh_list = get(groot, 'Children');
            fig_handle = [];
            
            % bad code writing style but fast
            for f = 1 : numel(fh_list)
                try
                    if isfield(fh_list(f).UserData, 'name') && strcmp(fh_list(f).UserData.name, GUI_Downloader.WIN_NAME)
                        fig_handle = fh_list(f);
                        break
                    end
                catch
                end
            end
            if ~isempty(fig_handle)
                delete(fig_handle);
            end
        end

    end
    %% METHODS INIT
    % ==================================================================================================================================================
    methods                
        function init(this)
        end
                        
        function openGUI(this, flag_wait)
            % WIN CONFIGURATION
            % L| N|    W
            %
            %
            % ----------
            % b      b b
            %
                                    
            t0 = tic();
            state = Core.getCurrentSettings;
            this.ok_go = false;

            log = Core.getLogger;

            % Get the old goGPS window
            if flag_wait
               % If goGPS is started rec have been reinitialized, close the inspector then
                fh_list = get(groot, 'Children');
                
                % bad code writing style but fast
                for f = 1 : numel(fh_list)
                    try
                        if isfield(fh_list(f).UserData, 'name') && strcmp(fh_list(f).UserData.name, GUI_Inspector.WIN_NAME)
                            % If there are lone Edit figures close them
                            close(fh_list(f));
                        end
                    catch ex
                        Core_Utils.printEx(ex);
                    end
                end
            end
            old_win = this.getUniqueWinHandle();

            if ~isempty(old_win)
                log.addMarkedMessage('Resetting the old Downloader Window');
                win = old_win;
                this.go_but.Enable = 'on';
                
                if strcmp(this.win.Visible, 'off')
                    try
                        msg_pos = GUI_Msg.getPosition;
                        win.Position(1) = sum(msg_pos([1 3]));
                        win.Position(2) = sum(msg_pos([2 4])) - win.Position(4);
                    catch
                        if isunix && not(ismac())
                            win.Position(1) = round((win.Parent.ScreenSize(3) - win.Position(3)) / 2);
                            win.Position(2) = round((win.Parent.ScreenSize(4) - win.Position(4)) / 2);
                        else
                            win.OuterPosition(1) = round((win.Parent.ScreenSize(3) - win.OuterPosition(3)) / 2);
                            win.OuterPosition(2) = round((win.Parent.ScreenSize(4) - win.OuterPosition(4)) / 2);
                        end
                    end
                end
            else
                log.addMarkedMessage('Opening a new Edit Settings Window');
                                                   
                % Main Window ----------------------------------------------------------------------------------------------

                win = figure( 'Name', sprintf('goDownloader using %s @ %s', state.getPrjName, state.getHomeDir), ...
                    'Visible', 'off', ...
                    'DockControls', 'off', ...
                    'MenuBar', 'none', ...
                    'ToolBar', 'none', ...
                    'NumberTitle', 'off', ...
                    'Renderer', 'opengl', ...
                    'Position', [0 0 1040, 640]);
                win.UserData.name = this.WIN_NAME;
                % Center the window
                try
                    msg_pos = GUI_Msg.getPosition;
                    win.Position(1) = sum(msg_pos([1 3]));
                    win.Position(2) = sum(msg_pos([2 4])) - win.Position(4);
                catch
                    if isunix && not(ismac())
                        win.Position(1) = round((win.Parent.ScreenSize(3) - win.Position(3)) / 2);
                        win.Position(2) = round((win.Parent.ScreenSize(4) - win.Position(4)) / 2);
                    else
                        win.OuterPosition(1) = round((win.Parent.ScreenSize(3) - win.OuterPosition(3)) / 2);
                        win.OuterPosition(2) = round((win.Parent.ScreenSize(4) - win.OuterPosition(4)) / 2);
                    end
                end
                win.Visible = 'on';
                
                this.win = win;            

                % empty cur_lists
                this.check_boxes = {}; % List of all the checkboxes
                this.check_boxes_dwn = {}; % List of all the download checkboxes
                this.pop_ups = {};     % List of drop down menu
                this.rpop_up = {};     % Remote resources pup-up
                this.ropref = {};      % Remote Orbit Preferences
                this.ripref = {};      % Remote Iono preferences
                this.rv2pref = {};     % Remote VMF source preferences
                this.edit_texts = {};  % List of editable text
                this.edit_texts_array = {}; % list of editable text array
                this.flag_list = {};   % list of all the flags                

                try
                    main_bv = uix.VBox('Parent', win, ...
                        'Padding', 5, ...
                        'BackgroundColor', Core_UI.DARK_GREY_BG);
                catch ex
                    log.setOutMode(1,[],0); % to plot a Warning I need to disable GUI and enable
                    log.addError('Please install GUI Layout Toolbox (https://it.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox)');
                    open('GUI Layout Toolbox 2.3.4.mltbx');
                    log.newLine();
                    log.addWarning('After installation re-run goGPS');
                    delete(win);
                    return;
                end
                top_bh = uix.HBox( 'Parent', main_bv, 'BackgroundColor', Core_UI.DARK_GREY_BG);

                bottom_bh = uix.HBox( 'Parent', main_bv, ...
                    'Padding', 5, ...
                    'Spacing', 5, ...
                    'BackgroundColor', Core_UI.DARKER_GREY_BG);
                main_bv.Heights = [-1 30];

                left_bv = uix.VBox('Parent', top_bh, ...
                    'Padding', 5, ...
                    'BackgroundColor', Core_UI.DARK_GREY_BG);

                panel_g_border = uix.VBox('Parent', top_bh, ...
                    'Padding', 5, ...
                    'BackgroundColor', Core_UI.DARK_GREY_BG);
                %panel = uix.BoxPanel('Parent', panel_border, 'Title', 'Settings' );
                top_bh.Widths = [210 -1];
                
                % Set-up menu ----------------------------------------------------------------------------------------------

                this.addGoMenu();
                
                % Logo/title box -------------------------------------------------------------------------------------------

                Core_UI.insertLogoGUI(left_bv);
                left_bv.Heights = 94;
                
                % Left Panel -----------------------------------------------------------------------------------------------

                Core_UI.insertEmpty(left_bv, Core_UI.LIGHT_GREY_BG);
                
                % Session selector -----------------------------------------------------------------------------------------
                this.session_panel = Core_UI.insertPanelLight(left_bv, 'Sessions');
                sss_box_v = uix.VBox('Parent', this.session_panel, ...
                    'BackgroundColor', Core_UI.LIGHT_GREY_BG);
                sss_box_h = uix.HBox('Parent', sss_box_v, ...
                    'BackgroundColor', Core_UI.LIGHT_GREY_BG);
                
                sss_box_l = uix.VBox('Parent', sss_box_h, ...
                    'BackgroundColor', Core_UI.LIGHT_GREY_BG);
                
                date_g = uix.Grid( 'Parent', sss_box_l, ...
                    'BackgroundColor', Core_UI.LIGHT_GREY_BG);
                uicontrol('Parent', date_g, ...
                    'Style', 'Text', ...
                    'String', 'Start', ...
                    'FontSize', Core_UI.getFontSize(8), ...
                    'BackgroundColor', Core_UI.LIGHT_GREY_BG, ...
                    'ForegroundColor', Core_UI.BLACK);
                uicontrol('Parent', date_g, ...
                    'Style', 'Text', ...
                    'String', 'Stop', ...
                    'FontSize', Core_UI.getFontSize(8), ...
                    'BackgroundColor', Core_UI.LIGHT_GREY_BG, ...
                    'ForegroundColor', Core_UI.BLACK);
                ts = state.getSessionsStart();
                te = state.getSessionsStop();
                if te.isempty() || ts.isempty()
                    ts = GPS_Time.now();
                    te = GPS_Time.now();
                end
                this.ui_sss_start = Core_UI.insertVDateSpinnerHour(date_g, ts, @this.onSessionChange);
                this.ui_sss_stop = Core_UI.insertVDateSpinnerHour(date_g, te, @this.onSessionChange);
                date_g.Heights = [1, 1] .* 2.2 * Core_UI.LINE_HEIGHT;
                date_g.Widths = [46, -1];
                Core_UI.insertEmpty(left_bv, Core_UI.LIGHT_GREY_BG);
                left_bv.Heights = [left_bv.Heights(1) 10 128 -1];
                
                % Download selector ----------------------------------------------------------------------------------------

                dwn_panel = Core_UI.insertPanelLight(left_bv, 'Resources to Download');

                Core_UI.insertEmpty(left_bv, Core_UI.DARK_GREY_BG);
                left_bv.Heights = [left_bv.Heights(1:3)' 10 200 -1];
                vb = uix.VBox('Parent', dwn_panel, ...
                    'BackgroundColor', Core_UI.LIGHT_GREY_BG);
                this.check_boxes_dwn{end+1} = Core_UI.insertCheckBoxLight(vb, 'igs14.atx', 'chkbox_igsatx', []); this.check_boxes_dwn{end}.Value = false;                
                this.check_boxes_dwn{end+1} = Core_UI.insertCheckBoxLight(vb, 'CRX', 'chkbox_crx', []); this.check_boxes_dwn{end}.Value = true;
                this.check_boxes_dwn{end+1} = Core_UI.insertCheckBoxLight(vb, 'Ephemerides, clocks, ERP', 'chkbox_eph', []); this.check_boxes_dwn{end}.Value = true;
                this.check_boxes_dwn{end+1} = Core_UI.insertCheckBoxLight(vb, 'Ionospheric maps', 'chkbox_iono', []); this.check_boxes_dwn{end}.Value = true;
                this.check_boxes_dwn{end+1} = Core_UI.insertCheckBoxLight(vb, 'Broadcast Ionosphere', 'chkbox_iono_brdc', []); this.check_boxes_dwn{end}.Value = true;
                this.check_boxes_dwn{end+1} = Core_UI.insertCheckBoxLight(vb, 'Biases', 'chkbox_bias', []); this.check_boxes_dwn{end}.Value = true;
                this.check_boxes_dwn{end+1} = Core_UI.insertCheckBoxLight(vb, 'VMF', 'chkbox_vmf', []); this.check_boxes_dwn{end}.Value = true;
                this.check_boxes_dwn{end+1} = Core_UI.insertCheckBoxLight(vb, 'Athmospeheric loading', 'chkbox_atm', []); this.check_boxes_dwn{end}.Value = true;
                
                % Main Panel -----------------------------------------------------------------------------------------------
                
                wait_box = uix.VBox( 'Parent', panel_g_border, ...
                    'Padding', 200, ...
                    'BackgroundColor', Core_UI.LIGHT_GREY_BG);
                this.insertWaitBox(wait_box, 'Building interface...');
                % On linux I have to repeat this operation or the wait
                % box will not be centered, let's always do it, it is safer!
                wait_box.Padding = 190;
                drawnow;
                delete(wait_box)
                tab_panel = uix.TabPanel('Parent', panel_g_border, ...
                    'TabWidth', 100, ...
                    'Padding', 5, ...
                    'BackgroundColor', Core_UI.LIGHT_GREY_BG);
                                
                % Tab creation  --------------------------------------------------------------------------------------------
                
                % remote resource ini
                this.insertTabRemoteResource(tab_panel)
                                
                % Tabs settings --------------------------------------------------------------------------------------------

                tab_panel.TabTitles = {'Resources'};
                tab_panel.Selection = 1;
                
                % Botton Panel ---------------------------------------------------------------------------------------------
                
                bottom_bhl = uix.HButtonBox( 'Parent', bottom_bh, ...
                    'ButtonSize', [165 28] , ...
                    'Spacing', 5, ...
                    'HorizontalAlignment', 'left', ...
                    'BackgroundColor', Core_UI.DARKER_GREY_BG);

                ini_name_box = uix.HBox( 'Parent', bottom_bh, ...
                    'Padding', 2, ...
                    'BackgroundColor', Core_UI.DARKER_GREY_BG);

                uicontrol('Parent', ini_name_box, ...
                    'Style', 'Text', ...
                    'String', ' Current INI path:', ...
                    'ForegroundColor', Core_UI.LIGHT_GREY_BG, ...
                    'HorizontalAlignment', 'left', ...
                    'FontSize', Core_UI.getFontSize(8), ...
                    'BackgroundColor', Core_UI.DARKER_GREY_BG);   

                this.ini_path = uicontrol('Parent', ini_name_box, ...
                    'Style', 'Text', ...
                    'String', 'last_settings.ini', ...
                    'ForegroundColor', Core_UI.LIGHT_GREY_BG, ...
                    'HorizontalAlignment', 'left', ...
                    'FontSize', Core_UI.getFontSize(8), ...
                    'BackgroundColor', Core_UI.DARKER_GREY_BG);            

                ini_name_box.Widths = [100 -1];

                bottom_bhr = uix.HButtonBox( 'Parent', bottom_bh, ...
                    'Spacing', 5, ...
                    'ButtonSize', [185 28] , ...
                    'HorizontalAlignment', 'right', ...
                    'BackgroundColor', Core_UI.DARKER_GREY_BG);

                exit_but = uicontrol( 'Parent', bottom_bhl, ...
                    'String', 'Exit', ...
                    'Callback', @this.close); %#ok<NASGU>

                load_but = uicontrol( 'Parent', bottom_bhr, ...
                    'String', 'Load', ...
                    'Callback', @this.loadState); %#ok<NASGU>

                % Show go button only if I'm executing the interface from goGPS script
                this.go_but = uicontrol( 'Parent', bottom_bhr, ...
                    'String', 'Download!', ...
                    'FontAngle', 'italic', ...
                    'Enable', 'on', ...
                    'Callback', @this.go, ...
                    'FontWeight', 'bold');
                
                %session_height = sum(left_bv.Children(2).Children(1).Heights);
                bottom_bh.Widths = [60 -1 260];
                            
                set(win, 'CloseRequestFcn', @this.close);                
            end
            
            this.updateUI();            
                        
            this.win.Visible = 'on';

            t_win = toc(t0);
            cm = log.getColorMode();
            log.setColorMode(false);
            log.addStatusOk(sprintf('Downloader GUI initialization completed in %.2f seconds\n', t_win));
            log.setColorMode(cm);
            this.bringOnTop();            
        end
        

    end
    %% METHODS INSERT
    % ==================================================================================================================================================
    methods
        function insertWaitBox(this, container, status)
            % Insert a wait spinner
            % Thanks to undocumented MATLAB blog
            iconsClassName = 'com.mathworks.widgets.BusyAffordance$AffordanceSize';
            iconsSizeEnums = javaMethod('values',iconsClassName);
            SIZE_32x32 = iconsSizeEnums(2);  % (1) = 16x16,  (2) = 32x32
            j_spinner = com.mathworks.widgets.BusyAffordance(SIZE_32x32, 'Loading...');  % icon, label
            
            j_spinner.setPaintsWhenStopped(true);  % default = false
            j_spinner.useWhiteDots(false);         % default = false (true is good for dark backgrounds)
            % DEPRECATE!!!
            warning off
            tmp = javacomponent(j_spinner.getComponent, [0,0,80,150], container);
            warning on
            tmp.setBackground(java.awt.Color(Core_UI.LIGHT_GREY_BG(1),Core_UI.LIGHT_GREY_BG(2),Core_UI.LIGHT_GREY_BG(3)));
            
            j_spinner.start;
            j_spinner.setBusyText(status);
            drawnow
        end        
        
        function insertResources(this, container)
            resources_BG = Core_UI.LIGHT_GREY_BG;
            tab = uix.VBox('Parent', container, ...
                'Padding', 5, ...
                'BackgroundColor', resources_BG);
            
            uicontrol('Parent', tab, ...
                'Style', 'Text', ...
                'String', 'Select Computational Center:', ...
                'ForegroundColor', Core_UI.BLACK, ...
                'HorizontalAlignment', 'left', ...
                'FontSize', Core_UI.getFontSize(9), ...
                'BackgroundColor', resources_BG);
            
            this.uip.tab_res = tab;
        end
        
        function insertTabRemoteResource(this, container)
            tab = uix.VBox('Parent', container);
            
            state = Core.getCurrentSettings;
            tab_bv = uix.VBox( 'Parent', tab, ...
                'Spacing', 5, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            
%             uicontrol('Parent', tab_bv, ...
%                 'Style', 'Text', ...
%                 'HorizontalAlignment', 'left', ...
%                 'String', 'Remote Resources ini file contains download locations - not editable from GUI', ...
%                 'BackgroundColor', Core_UI.LIGHT_GREY_BG, ...
%                 'ForegroundColor', Core_UI.BLACK, ...
%                 'FontSize', Core_UI.getFontSize(10), ...
%                 'FontWeight', 'bold');
            
            uicontrol('Parent', tab_bv, ...
                'Style', 'Text', ...
                'HorizontalAlignment', 'left', ...
                'String', ['Remote resources INI path: ' state.getRemoteSourceFile], ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG, ...
                'ForegroundColor', 0.3 * ones(3, 1), ...
                'FontSize', Core_UI.getFontSize(7.5));
            
            Core_UI.insertHBarLight(tab_bv);
            
            try
                r_man = Remote_Resource_Manager.getInstance(state.getRemoteSourceFile());
                [tmp, this.rpop_up{end+1}] = Core_UI.insertPopUpLight(tab_bv, 'Orbit Center', r_man.getCenterListExtended, 'selected_orbit_center', @this.onResourcesPopUpChange, [200 -1]);                
            catch
                str = sprintf('[!!] Resource file missing:\n"%s"\nnot found\n\ngoGPS may not work properly', state.getRemoteSourceFile);
            end
            
            box_opref = uix.HBox( 'Parent', tab_bv, ...
                'Spacing', 5, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            
            uicontrol('Parent', box_opref, ...
                'Style', 'Text', ...
                'HorizontalAlignment', 'left', ...
                'String', 'Center orbit type preference', ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG, ...
                'ForegroundColor', Core_UI.BLACK, ...
                'FontSize', Core_UI.getFontSize(9));
          
            this.ropref = {};
            this.ropref{1} = Core_UI.insertCheckBoxLight(box_opref, 'Final', 'orbit1', @this.onResourcesPrefChange);
            this.ropref{2} = Core_UI.insertCheckBoxLight(box_opref, 'Rapid', 'orbit2', @this.onResourcesPrefChange);
            this.ropref{3} = Core_UI.insertCheckBoxLight(box_opref, 'Ultra rapid', 'orbit3', @this.onResourcesPrefChange);
            this.ropref{4} = Core_UI.insertCheckBoxLight(box_opref, 'Broadcast', 'orbit4', @this.onResourcesPrefChange);
            this.ropref{5} = Core_UI.insertCheckBoxLight(box_opref, 'Real-time', 'orbit5', @this.onResourcesPrefChange);
            box_opref.Widths = [250 -1 -1 -1 -1 -1];
            
            try
                r_man = Remote_Resource_Manager.getInstance(state.getRemoteSourceFile());
                [tmp, this.rpop_up{end+1}] = Core_UI.insertPopUpLight(tab_bv, 'Iono Center', r_man.getCenterListExtended(1), 'selected_iono_center', @this.onResourcesPopUpChange, [200 -1]);                
            catch
                str = sprintf('[!!] Resource file missing:\n"%s"\nnot found\n\ngoGPS may not work properly', state.getRemoteSourceFile);
            end
            
            box_v1pref = uix.HBox( 'Parent', tab_bv, ...
                'Spacing', 5, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            
            uicontrol('Parent', box_v1pref, ...
                'Style', 'Text', ...
                'HorizontalAlignment', 'left', ...
                'String', 'Center iono type preference', ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG, ...
                'ForegroundColor', Core_UI.BLACK, ...
                'FontSize', Core_UI.getFontSize(9));
            
            this.ripref = {};
            this.ripref{1} = Core_UI.insertCheckBoxLight(box_v1pref, 'Final', 'iono1', @this.onResourcesPrefChange);
            this.ripref{1}.TooltipString = 'Final ionospheric map';
            this.ripref{2} = Core_UI.insertCheckBoxLight(box_v1pref, 'Rapid', 'iono2', @this.onResourcesPrefChange);
            this.ripref{2}.TooltipString = 'Rapid ionospheric map';
            this.ripref{3} = Core_UI.insertCheckBoxLight(box_v1pref, 'P1', 'iono3', @this.onResourcesPrefChange);
            this.ripref{3}.TooltipString = 'Ionospheric map predicted one day ahead ';
            this.ripref{4} = Core_UI.insertCheckBoxLight(box_v1pref, 'P2', 'iono4', @this.onResourcesPrefChange);
            this.ripref{4}.TooltipString = 'Ionospheric map predicted two day ahead ';
            this.ripref{5} = Core_UI.insertCheckBoxLight(box_v1pref, 'Broadcast', 'iono5', @this.onResourcesPrefChange);
            this.ripref{5}.TooltipString = 'Klobuchar ionospheric parameters';
            box_v1pref.Widths = [250 -1 -1 -1 -1 -1];
                 
            % vmf source
            box_v2pref = uix.HBox( 'Parent', tab_bv, ...
                'Spacing', 5, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            
            uicontrol('Parent', box_v2pref, ...
                'Style', 'Text', ...
                'HorizontalAlignment', 'left', ...
                'String', 'VMF source preference', ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG, ...
                'ForegroundColor', Core_UI.BLACK, ...
                'FontSize', Core_UI.getFontSize(9));
            
            this.rv2pref = {};
            this.rv2pref{1} = Core_UI.insertCheckBoxLight(box_v2pref, 'Operational', 'vmfs1', @this.onResourcesPrefChange);
            this.rv2pref{2} = Core_UI.insertCheckBoxLight(box_v2pref, 'ERA-Interim', 'vmfs2', @this.onResourcesPrefChange);
            this.rv2pref{3} = Core_UI.insertCheckBoxLight(box_v2pref, 'Forecast', 'vmfs3', @this.onResourcesPrefChange);
            box_v2pref.Widths = [250 -1 -1 -1 ];
            
            [~, this.pop_ups{end+1}] = Core_UI.insertPopUp(tab_bv, 'Mapping function', state.MF_LABEL, 'mapping_function', @this.onPopUpChange, [], Core_UI.LIGHT_GREY_BG);
            
            % Resource tree
            bottom_box = uix.VBox( 'Parent', tab_bv, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            
            rr_box = uix.VBox( 'Parent', bottom_box, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);

            but_line = uix.HButtonBox('Parent', rr_box, ...
                'ButtonSize', [165 28] , ...
                'Spacing', 5, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            
            uicontrol( 'Parent', but_line, ...
                'String', 'Open resource tree', ...
                'TooltipString', 'Open the inspector of the resources locations', ...
                'Callback', @this.openRRI);
            
            uicontrol( 'Parent', but_line, ...
                'String', 'Show orbit availability', ...
                'TooltipString', 'Show orbits and clock that are already on the disk', ...
                'Callback', @this.showOrbitsAvailability);
            
            uicontrol( 'Parent', but_line, ...
                'String', 'Download Orbits now', ...
                'TooltipString', 'Download orbits and clock if available', ...
                'Callback', @this.downloadOrbits);
            
            rr_box.Heights = [30];
            
            Core_UI.insertEmpty(bottom_box);
            
            dir_but_box = uix.VBox( 'Parent', bottom_box, ...
                'Spacing', 5, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            
            dir_reset = uicontrol( 'Parent', dir_but_box, ...
                'String', 'Reset all resource paths', ...
                'Callback', @this.resetResDir);

            dir_box = uix.VBox( 'Parent', bottom_box, ...
                'Spacing', 5, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            bottom_box.Heights = [-1 2 28 28*8];
                         
            [~, this.edit_texts{end + 1}, this.edit_texts{end + 2}, this.flag_list{end + 1}] = Core_UI.insertDirFileBox(dir_box, 'Antex (ATX) filename', 'atx_dir', 'atx_name', @this.onEditChange, [28 130 -3 5 -1 25]);
            [~, this.edit_texts{end + 1}, this.edit_texts{end + 2}, this.flag_list{end + 1}] = Core_UI.insertDirFileBox(dir_box, 'CRX path', 'crx_dir', 'crx_name', @this.onEditChange, [28 130 -3 5 -1 25]);
            [~, this.edit_texts{end + 1}, this.flag_list{end + 1}] = Core_UI.insertDirBox(dir_box, 'Eph local dir', 'eph_dir', @this.onEditChange, [28 130 -1 25]);
            [~, this.edit_texts{end + 1}, this.flag_list{end + 1}] = Core_UI.insertDirBox(dir_box, 'Clk local dir', 'clk_dir', @this.onEditChange, [28 130 -1 25]);
            [~, this.edit_texts{end + 1}, this.flag_list{end + 1}] = Core_UI.insertDirBox(dir_box, 'ERP local dir', 'erp_dir', @this.onEditChange, [28 130 -1 25]);
            [~, this.edit_texts{end + 1}, this.flag_list{end + 1}] = Core_UI.insertDirBox(dir_box, 'IONO local dir', 'iono_dir', @this.onEditChange, [28 130 -1 25]);
            [~, this.edit_texts{end + 1}, this.flag_list{end + 1}] = Core_UI.insertDirBox(dir_box, 'Biases local dir', 'bias_dir', @this.onEditChange, [28 130 -1 25]);
            [~, this.edit_texts{end + 1}, this.flag_list{end + 1}] = Core_UI.insertDirBox(dir_box, 'VMF local dir', 'vmf_dir', @this.onEditChange, [28 130 -1 25]);
            [~, this.edit_texts{end + 1}, this.flag_list{end + 1}] = Core_UI.insertDirBox(dir_box, 'ATM local dir', 'atm_load_dir', @this.onEditChange, [28 130 -1 25]);

            
            tab_bv.Heights = [15 2 22 18 22 18 18 22 -1];
            this.uip.tab_rr = tab;            
        end
                                        
        function updateCooTable(this)
            % Update the table of coordinates (CRD file interface)
            rf = Core.getReferenceFrame();
            if ~rf.isValid && (exist(Core.getState.getCrdFile, 'file') == 2)
                rf.init();
            end
            this.coo_tbl.Data = rf.getEntryCell();
            this.coo_tbl.RowName = {};
        end
        
        function insertTabProcessing(this, container)
            color_bg = Core_UI.LIGHT_GREY_BG_NOT_SO_LIGHT;
            
            tab_panel = uix.TabPanel('Parent', container, ...
                'TabWidth', 110, ...
                'Padding', 5, ...
                'BackgroundColor', color_bg, ...
                'SelectionChangedFcn', @this.onTabChange);
            
            % Insert tabs
            tab_sat = uix.VBox('Parent', tab_panel, ...
                'Padding', 2, ...
                'BackgroundColor', color_bg);
            tab_atm = uix.VBox('Parent', tab_panel, ...
                'Padding', 2, ...
                'BackgroundColor', color_bg);
            tab_prepro = uix.VBox('Parent', tab_panel, ...
                'Padding', 2, ...
                'BackgroundColor', color_bg);
            tab_generic = uix.VBox('Parent', tab_panel, ...
                'Padding', 2, ...
                'BackgroundColor', color_bg);
            tab_ppp = uix.VBox('Parent', tab_panel, ...
                'Padding', 2, ...
                'BackgroundColor', color_bg);
            tab_net = uix.VBox('Parent', tab_panel, ...
                'Padding', 2, ...
                'BackgroundColor', color_bg);
            tab_panel.TabTitles = {'Data selection', 'Atmosphere', 'Pre-Processing', 'Generic Options', 'PPP Parameters', 'NET Parameters'};
            
            this.insertDataSelection(tab_sat, color_bg);
            this.insertTabAtmosphere(tab_atm, color_bg)

            this.insertTabPrePro(tab_prepro, color_bg);            

            this.insertGenericOpt(tab_generic, color_bg);
            this.insertPPP(tab_ppp, color_bg);
            this.insertNET(tab_net, color_bg);
                        
            this.uip.tab_reg = tab_panel;
        end
        
        function insertDataSelection(this, tab_sat, color_bg)
            dopt_vbox = uix.VBox('Parent', tab_sat,...
                'BackgroundColor', color_bg);
            
            uicontrol('Parent', dopt_vbox, ...
                'Style', 'Text', ...
                'String', 'Data to keep during processing (if present in the receiver data)', ...
                'ForegroundColor', Core_UI.BLACK, ...
                'HorizontalAlignment', 'left', ...
                'FontSize', Core_UI.getFontSize(9), ...
                'BackgroundColor', color_bg);
            
            Core_UI.insertEmpty(dopt_vbox, color_bg);
            ss_panel  = this.insertSatSelector(dopt_vbox, color_bg); %#ok<NASGU>
            Core_UI.insertEmpty(dopt_vbox, color_bg);
            
            err_box_g = uix.VBox('Parent', dopt_vbox, ...
                'BackgroundColor', color_bg);

            dopt_vbox.Heights = [18 10 168 10 -1];

            field_dim = [280 40 5 50];
            [grd, this.edit_texts{end+1}] = Core_UI.insertEditBox(err_box_g, 'Min satellites per epoch', 'min_n_sat', 'n', @this.onEditChange, field_dim, color_bg);
            [grd, this.edit_texts{end+1}] = Core_UI.insertEditBox(err_box_g, 'Min percentage of required epochs [0-100]', 'min_p_epoch', '%', @this.onEditChange, field_dim, color_bg);
            ttip = 'This is not kept in case of snooping';
            if verLessThan('matlab','9.5')
                grd.Children(end).TooltipString = ttip;
            else
                grd.Children(end).Tooltip = ttip;
            end

            [~, this.edit_texts{end+1}] = Core_UI.insertEditBox(err_box_g, 'Data cut-off angle', 'cut_off', 'deg', @this.onEditChange, field_dim, color_bg);
            [~, this.edit_texts{end+1}] = Core_UI.insertEditBox(err_box_g, 'SNR absolute threshold', 'abs_snr_thr', 'dBHz', @this.onEditChange, field_dim, color_bg);
            [grd, this.edit_texts{end+1}] = Core_UI.insertEditBox(err_box_g, 'SNR scaled threshold', 'scaled_snr_thr', 'dBHz', @this.onEditChange, field_dim, color_bg);
            ttip = 'Different trackings have different scaling factor, rescale them w.r.t. the code error level of the first frequency/tracking';
            if verLessThan('matlab','9.5')
                grd.Children(end).TooltipString = ttip;
            else
                grd.Children(end).Tooltip = ttip;
            end
            [~, this.edit_texts{end+1}] = Core_UI.insertEditBox(err_box_g, 'Min arc length', 'min_arc', 's', @this.onEditChange, field_dim, color_bg);
            Core_UI.insertEmpty(err_box_g, color_bg);            
            err_box_g.Heights = [Core_UI.LINE_HEIGHT * ones(6,1); -1];            
        end        
    end
    %% METHODS getters
    % ==================================================================================================================================================
    methods
        function ok_go = isGo(this)
            ok_go = this.ok_go;
        end
    end
    %% METHODS UI getters
    % ==================================================================================================================================================
    methods
        function [sss_start, sss_stop, validity_check] = getSessionsLimits(this)
            % get Start session and stop from UI
            % check validity if sss_start > sss_stop then stop = start
            %
            % SYNTAX:
            %   [sss_start, sss_stop, validity_check] = getSessionLimits(this)
            %
            state = Core.getCurrentSettings();
            validity_check = true;
            
            try
                date = this.ui_sss_start.Children(2).JavaPeer.getDate;
            catch
                date = [];
            end
            if isempty(date)
                sss_start = state.getSessionsStartExt;
            else
                sss_start = GPS_Time([date.getYear+1900 (date.getMonth + 1) date.getDate 0 0 0]);
            end
            try
                hh_mm_ss = this.ui_sss_start.Children(1).Children(1).String;
            catch
                hh_mm_ss = [];
            end
            if ~isempty(hh_mm_ss)
                time_parts = regexp(hh_mm_ss,'(?<hour>\d+):(?<minute>\d+):(?<second>\d+)','names');
                if isempty(time_parts)
                    time_parts = struct('hour', '00', 'minute', '00', 'second', '00');
                end
                sss_start.addSeconds(str2num(time_parts.hour)*3600 + str2num(time_parts.minute)*60 + str2num(time_parts.second));
                this.ui_sss_start.Children(1).Children(1).String = sss_start.toString('HH:MM:SS');
                date = this.ui_sss_stop.Children(2).JavaPeer.getDate;
            end
            if isempty(date)
                sss_stop = state.getSessionsStopExt;
            else
                sss_stop = GPS_Time([date.getYear+1900 (date.getMonth + 1) date.getDate 00 00 00]);
            end
            try
                hh_mm_ss = this.ui_sss_stop.Children(1).Children(1).String;
            catch
                hh_mm_ss = [];
            end
            if ~isempty(hh_mm_ss)
                time_parts = regexp(hh_mm_ss,'(?<hour>\d+):(?<minute>\d+):(?<second>\d+)','names');
                if isempty(time_parts)
                    time_parts = struct('hour', '23', 'minute', '59', 'second', '59');
                end
                sss_stop.addSeconds(str2num(time_parts.hour)*3600 + str2num(time_parts.minute)*60 + str2num(time_parts.second));
                this.ui_sss_stop.Children(1).Children(1).String = sss_stop.toString('HH:MM:SS');
            end
            if sss_stop <= sss_start
                validity_check = false;
                sss_stop = GPS_Time(floor(sss_start.getMatlabTime) + 86399/86400);
                this.ui_sss_stop.Children(1).Children(1).String = sss_stop.toString('HH:MM:SS');
            end
        end
        
        function type = getDownloadableItems(this)
            % Get from GUI the resources to be downloaded (or checked)
            type = {'igsatx', 'crx', 'eph', 'iono', 'iono_brdc', 'bias', 'vmf', 'atm'};
            flag_type = false(numel(type), 1);
            for i = 1 : numel(type)
                flag_type(i) = this.check_boxes_dwn{i}.Enable(2) == 'n' && this.check_boxes_dwn{i}.Value;
            end
            type = type(flag_type);
        end
    end
    %% METHODS EVENTS
    % ==================================================================================================================================================
    methods (Access = public)
        function onSessionChange(this, caller, event)
            % Manage the event of session modification (UI)
            %
            % SYNTAX:
            %   this.onSessionChange()
            %
            [sss_start, sss_stop, validity_check] = getSessionsLimits(this);
            
            try
                if ~validity_check
                    this.ui_sss_stop.Children(2).JavaPeer.setDate(java.util.Date(sss_stop.toString('yyyy/mm/dd')));
                end
                
                state = Core.getCurrentSettings();
                status_change = false;
                if sss_start - state.getSessionsStart() ~= 0
                    status_change = true;
                    state.setSessionStart(sss_start);
                end
                if sss_stop - state.getSessionsStop() ~= 0
                    status_change = true;
                    state.setSessionStop(sss_stop);
                end
                if status_change
                    this.checkFlag();
                end
            catch
                % The object is still being built
            end
        end
        
        function onCheckBoxChange(this, caller, event)
            Core.getCurrentSettings.setProperty(caller.UserData, caller.Value);
            this.updateCheckBoxFromState(); % refresh duplicated checkboxes
            if strcmp(caller.UserData, 'flag_clock_align')
                clear Core_Sky
            end
                
        end
        
        function onPopUpChange(this, caller, event)
            if isprop(Core.getCurrentSettings,[upper(caller.UserData) '_UI2INI'])
                value = Core.getCurrentSettings.([upper(caller.UserData) '_UI2INI'])(caller.Value);
            elseif strcmpi(caller.UserData(1:6),'tparam') && isempty(strfind(caller.UserData,'ztd')) && isempty(strfind(caller.UserData,'grad'))
                value = Core.getCurrentSettings.(['TIME_PARAMETRIZATION_UI2INI'])(caller.Value);
            else
                value = caller.Value;
            end
            Core.getCurrentSettings.setProperty(caller.UserData, value);
            % Update VMF Resolution Preferences
            if strcmpi(caller.UserData,'mapping_function')
                if (caller.Value == Prj_Settings.MF_VMF1 || caller.Value == Prj_Settings.MF_VMF3_1 || caller.Value == Prj_Settings.MF_VMF3_5)
                    r_man = Remote_Resource_Manager.getInstance();
                    state = Core.getCurrentSettings();
                    
                    % Update VMF Source Preferences
                    available_orbit = r_man.getVMFSourceType();
                    flag_preferred_orbit = true(3,1);
                    for i = 1 : 3
                        this.rv2pref{i}.Enable = iif(available_orbit(i), 'on', 'off');
                        flag_preferred_orbit(i) = available_orbit(i) && logical(this.rv2pref{i}.Value);
                    end
                    state.setPreferredVMFSource(flag_preferred_orbit)
                else
                    for i = 1 : 3
                        this.rv2pref{i}.Enable = 'off';
                    end
                end
                
                % Enable VMF if used
                this.check_boxes_dwn{6}.Enable = iif(Core.getState.isVMF,'on', 'off');
            end
        end
        
        function onResourcesPopUpChange(this, caller, event)
            
            state = Core.getCurrentSettings;
            if strcmp(caller.UserData, 'selected_orbit_center')
                % Particular case selected_center is in GUI with full description of the center
                % Use caller.Value and r_man.getCenterList();
                r_man = Remote_Resource_Manager.getInstance();
                
                % read current center
                [center_list, center_ss] = r_man.getCenterList();
                state.setCurCenter(center_list{caller.Value});
                this.updateUI;
            elseif strcmp(caller.UserData, 'selected_iono_center')
                % Particular case selected_center is in GUI with full description of the center
                % Use caller.Value and r_man.getCenterList();
                r_man = Remote_Resource_Manager.getInstance();
                
                % read current center
                [center_list, center_ss] = r_man.getCenterList(1);
                state.setProperty(caller.UserData, center_list{caller.Value});
            else
                state.setProperty(caller.UserData, caller.String(caller.Value));
            end
            
            % Set resources preferences
            r_man = Remote_Resource_Manager.getInstance();
            
            % Update Iono Preferences
            available_iono = r_man.getIonoType();
            flag_preferred_iono = true(numel(this.ripref),1);
            for i = 1 : numel(this.ripref)
                this.ripref{i}.Enable = iif(available_iono(i), 'on', 'off');
                flag_preferred_iono(i) = available_iono(i) && logical(this.ripref{i}.Value);
            end
            state.setPreferredIono(flag_preferred_iono)
            
            % Update Orbit Preferences
            available_orbit = r_man.getOrbitType(state.getRemoteCenter());
            flag_preferred_orbit = true(5,1);
            for i = 1 : numel(this.ropref)
                this.ropref{i}.Enable = iif(available_orbit(i), 'on', 'off');
                flag_preferred_orbit(i) = available_orbit(i) && logical(this.ropref{i}.Value);
            end
            state.setPreferredOrbit(flag_preferred_orbit)
            
            this.updateResourcePopUpsState();
        end
        
        function onResourcesPrefChange(this, caller, event)
            % Set resources preferences          
            r_man = Remote_Resource_Manager.getInstance();
            state = Core.getCurrentSettings;
            if strcmp(caller.UserData(1:4), 'iono')
                % Update Iono Preferences
                available_iono = r_man.getIonoType();
                flag_preferred_iono = true(numel(this.ripref),1);
                for i = 1 : numel(this.ripref)
                    flag_preferred_iono(i) = available_iono(i) && logical(this.ripref{i}.Value);
                end
                state.setPreferredIono(flag_preferred_iono)
            elseif strcmp(caller.UserData(1:4), 'orbi')
                % Update Orbit Preferences
                available_orbit = r_man.getOrbitType(state.getRemoteCenter());
                flag_preferred_orbit = true(5,1);
                for i = 1 : numel(this.ropref)
                    flag_preferred_orbit(i) = available_orbit(i) && logical(this.ropref{i}.Value);
                end
                state.setPreferredOrbit(flag_preferred_orbit)
            elseif strcmp(caller.UserData(1:4), 'vmfs')
                % Update Orbit Preferences
                available_orbit = r_man.getVMFSourceType();
                flag_preferred_orbit = true(3,1);
                for i = 1 : numel(this.rv2pref)
                    flag_preferred_orbit(i) = available_orbit(i) && logical(this.rv2pref{i}.Value);
                end
                state.setPreferredVMFSource(flag_preferred_orbit)
            end
                                    
            this.updateResourcePopUpsState();
        end
        
        function onEditChange(this, caller, event)
            state = Core.getCurrentSettings;
            prop = state.getProperty(caller.UserData);
            if ~isnumeric(prop)
                state.setProperty(caller.UserData, caller.String);
            else
                state.setProperty(caller.UserData, str2num(caller.String));
            end
            
            state.check();
            caller.String = state.getProperty(caller.UserData);            
            this.checkFlag();
            
            if strcmp(caller.UserData, 'crd_name') || strcmp(caller.UserData, 'crd_dir')
                rf = Core.getReferenceFrame;
                rf.init(state.getCrdFile);
                this.updateCooTable();
            end
            if strcmp(caller.UserData, 'obs_name') || strcmp(caller.UserData, 'obs_dir')
                this.updateRecList()
            end
        end
        
        function onEditArrayChange(this, caller, event)
            state = Core.getCurrentSettings;
            prop = state.getProperty(caller.UserData);
            n_child = length(caller.Parent.Children);
            array = [];
            for i = n_child : -1 : 1
                child =  caller.Parent.Children(i);
                if strcmp(child.Style, 'edit')
                    val = str2num(child.String);
                    if isempty(val)
                        val =0;
                    end
                    array = [array val];
                end
            end
            state.setProperty(caller.UserData, array);
            state.check();
            this.updateEditArrayFromState(caller.Parent);
        end
        
        function openRRI(this, caller, event)
            GUI_Remote_Resources.getInstance(this.win);
        end
        
        function showOrbitsAvailability(this, caller, event)
            sky = Core.getCoreSky; 
            sky.showOrbitsAvailability;
        end
        
        function downloadOrbits(this, caller, event)
            fw = File_Wizard;
            fw.downloadResource('eph',Core.getState.getSessionsStartExt, Core.getState.getSessionsStopExt);
        end
        
        function resetResDir(this, caller, event)
            state = Core.getCurrentSettings;
            state.atx_dir = '';
            state.atx_name = '';
            state.geoid_dir = '';
            state.geoid_name = state.GEOID_NAME;
            state.crx_dir = '';
            state.crx_name = state.CRX_NAME;
            state.eph_dir = '';
            state.clk_dir = '';
            state.erp_dir = '';
            state.iono_dir = '';
            state.igrf_dir = '';
            state.bias_dir = '';
            state.vmf_dir = '';
            state.atm_load_dir = '';
            
            state.check();
            this.updateUI();
        end
        
        function checkFlag(this)
            Core_UI.checkFlag(this.flag_list)
        end
        
        function updateSessionFromState(this, caller, event)
            state = Core.getCurrentSettings();
            this.ui_sss_start.Children(2).JavaPeer.setDate(java.util.Date(state.sss_date_start.toString('yyyy/mm/dd')));
            this.ui_sss_start.Children(1).Children(1).String = state.sss_date_start.toString('HH:MM:SS');
            %this.ui_sss_start.setDate(java.util.Date(state.sss_date_start.toString('yyyy/mm/dd')));
            this.ui_sss_stop.Children(2).JavaPeer.setDate(java.util.Date(state.sss_date_stop.toString('yyyy/mm/dd')));
            this.ui_sss_stop.Children(1).Children(1).String = state.sss_date_stop.toString('HH:MM:SS');
            %this.ui_sss_stop.setDate(java.util.Date(state.sss_date_stop.toString('yyyy/mm/dd')));
        end
                
        function updateCheckBoxFromState(this)
            for i = 1 : length(this.check_boxes)
                value = Core.getCurrentSettings.getProperty(this.check_boxes{i}.UserData);
                if ~isempty(value)
                    this.check_boxes{i}.Value = double(value(1));
                end
            end
        end
        
        function updateEditFromState(this)
            for i = 1 : length(this.edit_texts)
                value = Core.getCurrentSettings.getProperty(this.edit_texts{i}.UserData);
                if ~isempty(value)
                    this.edit_texts{i}.String = value;
                end
            end
        end
        
        function updateEditArrayFromState(this, array_box)
            name_prop = array_box.UserData;
            array_value = Core.getCurrentSettings.getProperty(name_prop);
            n_child = length(array_box.Children);
            n_val = length(array_value);
            j = 1;
            for i = n_child : -1 : 1
                child =  array_box.Children(i);
                if strcmp(child.Style, 'edit')
                    child.String = array_value(min(j,n_val));
                    j = j+1;
                end
            end
        end
        
        function updateEditArraysFromState(this)
            for i = 1 : length(this.edit_texts_array)
                this.updateEditArrayFromState(this.edit_texts_array{i});
            end
        end
        
        function updatePopUpsState(this)
            state = Core.getCurrentSettings;
            for i = 1 : length(this.pop_ups)
                value = state.getProperty(this.pop_ups{i}.UserData);
                if ~isempty(value)
                    if  isprop(state,[upper(this.pop_ups{i}.UserData) '_UI2INI'])
                        this.pop_ups{i}.Value = find(state.([upper(this.pop_ups{i}.UserData) '_UI2INI']) == value);
                    elseif strcmpi(this.pop_ups{i}.UserData(1:6),'tparam')  && isempty(strfind(this.pop_ups{i}.UserData,'ztd')) && isempty(strfind(this.pop_ups{i}.UserData,'grad'))
                        this.pop_ups{i}.Value = find(state.(['TIME_PARAMETRIZATION_UI2INI']) == value);
                    else
                        this.pop_ups{i}.Value = value;
                    end
                end
            end
        end
        
        function updateResourcePopUpsState(this)
            % Getting current remote resource manager
            r_man = Remote_Resource_Manager.getInstance();
            
            state = Core.getCurrentSettings;
            % read current orbit center
            [center_list, center_ss] = r_man.getCenterList();
            cur_center = state.getCurCenter;
            if isempty(cur_center)
                cur_center = {'default'};
            end
            value = 1;
            while (value < numel(center_list)) && ~strcmp(center_list{value}, cur_center)
                value = value + 1;
            end
            
            % display resources tree of the current center
            if ~isempty(value)
                this.rpop_up{1}.Value = value;
            end
            
            % Update constellation Available for the center
            this.rpop_up{1}.Parent.Children(2).String = sprintf('Supported satellites: "%s"', center_ss{value});
            
            % Update Orbit Preferences
            available_orbit = r_man.getOrbitType(cur_center{1});
            for i = 1 : 5
                this.ropref{i}.Enable = iif(available_orbit(i), 'on', 'off');
            end
            flag_preferred_orbit = state.getPreferredOrbit();
            for i = 1 : 5
                if available_orbit(i)
                    this.ropref{i}.Value = this.ropref{i}.Value | flag_preferred_orbit(i);
                end
            end
            
             % read current iono center
            [center_list, center_ss] = r_man.getCenterList(1);
            cur_center = state.getCurIonoCenter;
            if isempty(cur_center)
                cur_center = {'default'};
            end
            value = 1;
            while (value < numel(center_list)) && ~strcmp(center_list{value}, cur_center)
                value = value + 1;
            end
            
            % display resources tree of the current center
            if ~isempty(value)
                this.rpop_up{2}.Value = value;
            end
%             
%             % Update constellation Available for the center
%             this.rpop_up{2}.Parent.Children(2).String = sprintf('Supported satellites: "%s"', center_ss{value});
%             
            
            % Update Iono Preferences
            available_iono = r_man.getIonoType(cur_center{1});
            for i = 1 :  numel(this.ripref)
                this.ripref{i}.Enable = iif(available_iono(i), 'on', 'off');
            end
            flag_preferred_iono = state.getPreferredIono();            
            for i = 1 :  numel(this.ripref)
                if available_iono(i)
                    this.ripref{i}.Value = this.ripref{i}.Value | flag_preferred_iono(i);
                end
            end
            % Update VMF source Preferences
            available_vmf = r_man.getVMFSourceType();
            for i = 1 : 3
                this.rv2pref{i}.Enable = iif(available_vmf(i), 'on', 'off');
            end
            flag_preferred_vmf = state.getPreferredVMFSource();            
            for i = 1 : 3
                if available_vmf(i)
                    this.rv2pref{i}.Value = this.rv2pref{i}.Value | flag_preferred_vmf(i);
                end
            end
            
        end
                
        function openInspector(this, caller, event)
            % Create a new project            
            goInspector;
        end
        
        function about(this, caller, event)
            % Show About window
            new = GUI_About(this);
        end
             
        function runGoGPS(this, caller, event)
            % Run goGPS
            this.close;
            goGPS();
        end
        
        function loadState(this, caller, event)
            % Load state settings
            
            state = Core.getCurrentSettings;
            config_dir = state.getHomeDir();
            if exist([config_dir filesep 'config'], 'dir')
                config_dir = [config_dir filesep 'config'];
            end
            % On MacOS this doesn't work anymore: [file_name, pathname] = uigetfile({'*.ini;','INI configuration file (*.ini)'; '*.mat;','state file goGPS < 0.5 (*.mat)'}, 'Choose file with saved settings', config_dir);
            [file_name, path_name] = uigetfile('*.ini', 'Choose file with saved settings', config_dir);
            
            if path_name ~= 0 % if the user pressed cancelled, then we exit this callback
                % get the extension (mat/ini):
                [~, ~, ext] = fileparts(file_name);
                
                % build the path name of the file to be loaded
                settings_file = fullfile(path_name, file_name);
                if strcmp(ext, '.ini')
                    state.importIniFile(settings_file);
                    Core.getReferenceFrame.init();
                    this.updateUI();
                    this.updateRecList();
                else
                    Core.getLogger.addError('Unrecognized input file format!');
                end
            end
        end
        
        function close(this, caller, event)
            
            % This is closing definitively, I prefer to hide this interface, it takes a while to restore it
            %delete(this.win);
            if isvalid(this.win)
                this.win.Visible = 'off';
                %uiresume(this.win);
            end            
        end
        
        function go(this, caller, event)
            
            %if isvalid(this.win)
            %    this.win.Visible = 'off';
            %end
            core = Core.getCurrentCore();
            core.log.addMarkedMessage('Starting download!');
            fw = File_Wizard;
            [sss_start, sss_stop] = this.getSessionsLimits();
            fw.downloadResource(this.getDownloadableItems(),sss_start, sss_stop);
            this.ok_go = true;
        end
        
        function updateUI(this)
            if isvalid(this.win)
                this.ini_path.String = Core.getCurrentSettings.getIniPath();
                this.updateSessionGUI();
                this.updateSessionFromState();
                this.updateCheckBoxFromState();
                this.updateEditFromState();
                this.updateEditArraysFromState();
                this.updatePopUpsState();
                this.updateResourcePopUpsState();
                this.checkFlag();
                
                % Enable VMF if used
                this.check_boxes_dwn{6}.Enable = iif(Core.getState.isVMF,'on', 'off');
            end
        end
                
        function updateSessionGUI(this)
            % enable disable fields
            ui_tspan = findobj(this.win, 'Tag', 'sss_duration');
            ui_buffer = findobj(this.win, 'Tag', 'sss_buffer');
            ui_smooth_tropo = findobj(this.win, 'Tag', 'sss_smooth');
            if Core.getCurrentSettings.isRinexSession()
                Core_UI.disableElement(ui_tspan);
                Core_UI.disableElement(ui_buffer);
                Core_UI.disableElement(ui_smooth_tropo);
            else
                Core_UI.enableElement(ui_tspan);
                Core_UI.enableElement(ui_buffer);                
                Core_UI.enableElement(ui_smooth_tropo);                
            end            
        end
        
        function setCheckBox(this, name_prop, value)
            for i = 1 : length(this.check_boxes)
                if this.check_boxes{i}.isvalid && strcmp(name_prop, this.check_boxes{i}.UserData)
                    this.check_boxes{i}.Value = double(value);
                end
            end
        end
    end
    
    methods
        function addGoMenu(this)
            this.menu.goGPS = uimenu(this.win, 'Label', 'goGPS');
            uimenu(this.menu.goGPS, ...
                'Label', 'About', ...
                'Callback', @this.about);
            uimenu(this.menu.goGPS, ...
                'Label', 'Open Inspector', ...
                'Callback', @this.openInspector);
            uimenu(this.menu.goGPS, ...
                'Label', 'Open goGPS', ...
                'Callback', @this.runGoGPS);
        end
    end
    
    methods (Static)
        function clr = getColorTrck(sys_c, band, parity)
            Gcol = '#a6cee3';
            Rcol = '#1f78b4';
            Ecol = '#b2df8a';
            Jcol = '#33a02c';
            Icol = '#fb9a99';
            Scol = '#e31a1c';
            Ccol = '#fdbf6f';
            if sys_c == 'G'
                clr = Gcol;
                rgb = hex2rgb(clr);
                
                if band == '1'
                    
                elseif band == '2'
                    rgb = max(rgb -20/255,0);
                    
                elseif band == '5'
                    rgb = max(rgb -10/255,0);
                    
                end
            elseif sys_c == 'R'
                clr = Rcol;
                rgb = hex2rgb(clr);
                
                if band == '1'
                    
                elseif band == '2'
                    rgb = max(rgb -20/255,0);
                    
                elseif band == '3'
                    rgb = max(rgb -10/255,0);
                    
                end
            elseif sys_c == 'E'
                
                clr = Ecol;
                rgb = hex2rgb(clr);
                
                if band == '1'
                    
                elseif band == '5'
                    rgb = max(rgb -20/255,0);
                    
                elseif band == '7'
                    rgb = max(rgb -10/255,0);
                elseif band == '6'
                    rgb = min(rgb +20/255,1);
                end
            elseif sys_c == 'C'
                clr = Ccol;
                rgb = hex2rgb(clr);
                
                if band == '2'
                elseif band == '7'
                    rgb = max(rgb -20/255,0);
                elseif band == '6'
                    rgb = max(rgb -10/255,0);
                end
            elseif sys_c == 'J'
                clr = Jcol;
                rgb = hex2rgb(clr);
                
                if band == '1'
                elseif band == '2'
                    rgb = max(rgb -20/255,0);
                elseif band == '5'
                    rgb = max(rgb -10/255,0);
                elseif band == '6'
                    rgb = min(rgb +20/255,1);
                end
            elseif sys_c == 'I'
                clr = Icol;
                rgb = hex2rgb(clr);
                
                if band == '5'
                elseif band == '9'
                    rgb = max(rgb -20/255,0);
                end
            elseif sys_c == 'S'
                clr = Icol;
                rgb = hex2rgb(clr);
                
                if band == '1'
                elseif band == '5'
                    rgb = max(rgb -20/255,0);
                end
            else
                clr = '#fcfcfc';
                rgb = hex2rgb(clr);
                
            end
            if parity
                rgb = min(rgb +6/255,1);
            else
                rgb = max(rgb -6/255,0);
                
            end
            clr = rgb2hex(rgb);
        end
    end
end
