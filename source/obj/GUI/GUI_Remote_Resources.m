%   CLASS GUI_Remote_Resources
% =========================================================================
%
% DESCRIPTION
%   class to manage the about window of goGPSz
%
% EXAMPLE
%   ui = GUI_Remote_Resources.getInstance();
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
%  Copyright (C) 2021 (GReD srl) Andrea Gatti
%  Written by:        Andrea Gatti
%  Contributors:      Andrea Gatti, ...
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

classdef GUI_Remote_Resources < GUI_Unique_Win
    
    properties (Constant)
        WIN_NAME = 'goGPS_Inspect_RR';
    end
    
    properties (Constant, Access = 'protected')
        BG_COLOR = Core_UI.DARK_GREY_BG;
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
        w_edt       % Handle of the main window
        win         % Handle to this window
        j_rrini     % ini resources file java pointer
    end
    
    %% PROPERTIES STATUS
    % ==================================================================================================================================================
    properties (GetAccess = private, SetAccess = private)
    end
    
    %% METHOD CREATOR
    % ==================================================================================================================================================
    
    methods (Static, Access = private)
        function this = GUI_Remote_Resources(w_edt)
            % GUI_MAIN object creator
            this.init();
            this.openGUI();
            if nargin == 1
                this.w_edt = w_edt;
            end
        end
    end
    
    methods (Static, Access = public)
        function this = getInstance(w_edt)
            % Get the persistent instance of the class
            persistent unique_instance_gui_rr__
            
            if nargin == 0
                w_edt = [];
            end
            
            if isempty(unique_instance_gui_rr__)
                this = GUI_Remote_Resources(w_edt);
                unique_instance_gui_rr__ = this;
            else
                this = unique_instance_gui_rr__;
                this.init();
                this.openGUI();
            end
        end
        
        function closeGUI()
            fh_list = get(groot, 'Children');
            fig_handle = [];
            
            % bad code writing style but fast
            for f = 1 : numel(fh_list)
                try
                    if isfield(fh_list(f).UserData, 'name') && strcmp(fh_list(f).UserData.name, GUI_Remote_Resources.WIN_NAME)
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
            this.log = Core.getLogger();
            this.state = Core.getState();
        end
        
        function openGUI(this)
            % Main Window ----------------------------------------------------------------------------------------------
            
            old_win = this.getUniqueWinHandle();
            if ~isempty(old_win)
                log = Core.getLogger();
                log.addMarkedMessage('Resetting the old Inspector window');
                win = old_win;
                
                if strcmp(this.win.Visible, 'off')
                    if isunix && not(ismac())
                        win.Position(1) = round((win.Parent.ScreenSize(3) - win.Position(3)) / 2);
                        win.Position(2) = round((win.Parent.ScreenSize(4) - win.Position(4)) / 2);
                    else
                        win.OuterPosition(1) = round((win.Parent.ScreenSize(3) - win.OuterPosition(3)) / 2);
                        win.OuterPosition(2) = round((win.Parent.ScreenSize(4) - win.OuterPosition(4)) / 2);
                    end
                end
                this.win = win;               
            else
                win = figure( 'Name', 'Remote Resources Inspector', ...
                    'Visible', 'on', ...
                    'DockControls', 'off', ...
                    'MenuBar', 'none', ...
                    'ToolBar', 'none', ...
                    'NumberTitle', 'off', ...
                    'Position', [0 0 1040, 640], ...
                    'Resize', 'on');
                win.UserData.name = this.WIN_NAME;

                this.win = win;

                if isunix && not(ismac())
                    win.Position(1) = round((win.Parent.ScreenSize(3) - win.Position(3)) / 2);
                    win.Position(2) = round((win.Parent.ScreenSize(4) - win.Position(4)) / 2);
                else
                    win.OuterPosition(1) = round((win.Parent.ScreenSize(3) - win.OuterPosition(3)) / 2);
                    win.OuterPosition(2) = round((win.Parent.ScreenSize(4) - win.OuterPosition(4)) / 2);
                end

                try
                    main_vb = uix.VBox('Parent', win, ...
                        'Padding', 5, ...
                        'BackgroundColor', Core_UI.DARKER_GREY_BG);
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
                top_bh = uix.HBox('Parent', main_vb);

                logo_GUI.BG_COLOR = Core_UI.DARK_GREY_BG;
                left_tbv = uix.VBox('Parent', top_bh, ...
                    'BackgroundColor', logo_GUI.BG_COLOR, ...
                    'Padding', 5);

                % Logo/title box -------------------------------------------------------------------------------------------

                logo_g = uix.Grid('Parent', left_tbv, ...
                    'Padding', 5, ...
                    'BackgroundColor', logo_GUI.BG_COLOR);

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

                Core_UI.insertEmpty(left_tbv, logo_GUI.BG_COLOR);
                left_tbv.Heights = [82 -1];

                % Title Panel -----------------------------------------------------------------------------------------------
                right_tvb = uix.VBox('Parent', top_bh, ...
                    'Padding', 5, ...
                    'BackgroundColor', logo_GUI.BG_COLOR);

                top_bh.Widths = [106 -1];

                title = uix.HBox('Parent', right_tvb, ...
                    'BackgroundColor', logo_GUI.BG_COLOR);
                txt = this.insertBoldText(title, 'Remote Resources in use', 10, Core_UI.LBLUE, 'left');
                txt.BackgroundColor = logo_GUI.BG_COLOR;

                Core_UI.insertEmpty(right_tvb, logo_GUI.BG_COLOR);

                % Disclaimer Panel -----------------------------------------------------------------------------------------------
                txt = this.insertText(right_tvb, {'Inspect the remoute resources locations', '(read-only)'},  8, [], 'left');

                right_tvb.Heights = [20 5 -1];

                string_bh = uix.VBox('Parent', main_vb, ...
                    'Padding', 10, ...
                    'BackgroundColor', GUI_Remote_Resources.BG_COLOR);

                %Core_UI.insertEmpty(string_bh);

                this.j_rrini = this.insertRRBox(string_bh);
                this.j_rrini.setText('RR');
                

                % Manage dimension -------------------------------------------------------------------------------------------
                main_vb.Heights = [80 -1];
            end
            this.updateResourcePopUpsState();
            this.win.Visible = 'on';
            this.bringOnTop();
        end
        
        function updateResourcePopUpsState(this)
            % Getting current remote resource manager
            r_man = Remote_Resource_Manager.getInstance();
            
            state = Core.getCurrentSettings;
            % read current center
            [center_list, center_ss] = r_man.getCenterList();
            cur_center = Core.getState.getCurCenter;
            if isempty(cur_center)
                cur_center = {'default'};
            end
            value = 1;
            while (value < numel(center_list)) && ~strcmp(center_list{value}, cur_center)
                value = value + 1;
            end
            
            % display resources tree of the current center
            if ~isempty(value)
                try
                    str = r_man.centerToString(state.getRemoteCenter());
                    str = strrep(['% ' str], char(10), [char(10) ' ']);
                catch
                    str = sprintf('[!!] Resource file missing:\n"%s"\nnot found\n\ngoGPS may not work properly', state.getRemoteSourceFile);
                end
                this.j_rrini.setText(str);
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
                'BackgroundColor', GUI_Remote_Resources.BG_COLOR);
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
                'BackgroundColor', GUI_Remote_Resources.BG_COLOR);
        end
        
        function insertLink(parent, url, prefix, prefix_color, link_label, font_size, color, alignment)
            if nargin < 4 || isempty(color)
                color = Core_UI.WHITE;
            end
            
            % Create and display the text label
            label_str = sprintf('<html><span style="color: rgb(%d,%d,%d);">%s </span><a href="" style="color: rgb(%d,%d,%d);">%s</a></html>', ...
                round(prefix_color(1) * 255), round(prefix_color(2) * 255), round(prefix_color(3) * 255), prefix, ...
                round(color(1) * 255), round(color(2) * 255), round(color(3) * 255), link_label);
            j_label = javaObjectEDT('javax.swing.JLabel', label_str);
            bg_color = num2cell(GUI_Remote_Resources.BG_COLOR);
            j_label.setBackground(java.awt.Color(bg_color{:}));
            %% DEPRECATE!!!
            warning on
            [hj_label, h_container] = javacomponent(j_label, [10,10,250,20], parent);
            warning off
            
            % Modify the mouse cursor when hovering on the label
            hj_label.setCursor(java.awt.Cursor.getPredefinedCursor(java.awt.Cursor.HAND_CURSOR));
            
            % Set the label's tooltip
            hj_label.setToolTipText(['Visit the ' url ' website']);
            
            % Set the mouse-click callback
            set(hj_label, 'MouseClickedCallback', @(h,e)web(['http://' url], '-browser'))
        end
        
        function j_rrini = insertRRBox(container)
            j_rrini = com.mathworks.widgets.SyntaxTextPane;
            %codeType = j_rrini.M_MIME_TYPE;  % j_rrini.contentType='text/m-MATLAB'
            %j_rrini.setContentType(codeType);
            str = ' Remote Resources in use: cd ';
            j_rrini.setText(str);
            % Create the ScrollPanel containing the widget
            j_scroll_settings = com.mathworks.mwswing.MJScrollPane(j_rrini);
            % Inject edit box with the Java Scroll Pane into the main_window
            %% DEPRECATE!!!
            warning off
            [panel_j, panel_h] = javacomponent(j_scroll_settings, [1 1 1 1], container);
            warning on
            j_rrini.setEditable(0);
        end
        
    end
    
    %% METHODS EVENTS
    % ==================================================================================================================================================
    methods (Access = public)
    end
end
