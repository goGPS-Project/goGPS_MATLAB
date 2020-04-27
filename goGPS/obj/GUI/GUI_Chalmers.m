%   CLASS GUI_Chalmers
% =========================================================================
%
% DESCRIPTION
%   class to manages the about window of goGPSz
%
% EXAMPLE
%   ui = GUI_Chalmers.getInstance();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Core_UI


%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b7
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

classdef GUI_Chalmers < handle
    
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
        w_main      % Handle of the main window
        win         % Handle to this window
    end
    
    %% PROPERTIES STATUS
    % ==================================================================================================================================================
    properties (GetAccess = private, SetAccess = private)
    end
    
    %% METHOD CREATOR
    % ==================================================================================================================================================
    methods (Static)
        function this = GUI_Chalmers(w_main)
            % GUI_MAIN object creator
            this.init();
            this.openGUI();
            if nargin == 1
                this.w_main = w_main;
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
            
            win = figure( 'Name', 'Chalmers Ocean Loading', ...
                'Visible', 'on', ...
                'DockControls', 'off', ...
                'MenuBar', 'none', ...
                'ToolBar', 'none', ...
                'NumberTitle', 'off', ...
                'Position', [0 0 640 400], ...
                'Resize', 'on');
            
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
            
            txt = this.insertBoldText(title, 'goGPS', 10, Core_UI.LBLUE, 'left');
            txt.BackgroundColor = logo_GUI.BG_COLOR;
            title_l = uix.VBox('Parent', title, 'BackgroundColor', GUI_Chalmers.BG_COLOR);
            title.Widths = [54 -1];
            Core_UI.insertEmpty(title_l, logo_GUI.BG_COLOR)
            txt = this.insertBoldText(title_l, ['- software V' Core.GO_GPS_VERSION], 8, [], 'left');
            txt.BackgroundColor = logo_GUI.BG_COLOR;
            title_l.Heights = [2, -1];
            
            % Disclaimer Panel -----------------------------------------------------------------------------------------------
            Core_UI.insertEmpty(right_tvb, logo_GUI.BG_COLOR)
            txt = this.insertText(right_tvb, {'Ocean loading computation must be required manually:', ...
                'go to Chalmers website and request a BLQ file using ocean tide model FES2004', ...
                ... % 'select also to compensate the values for the motion'...
                },  8, [], 'left');
            this.insertLink(right_tvb, 'holt.oso.chalmers.se/loading/', ...
                'Link:', Core_UI.WHITE, ...
                'http://holt.oso.chalmers.se/loading/', 8, Core_UI.LBLUE, 'left');
            
            right_tvb.Heights = [20 3 -1 20];
            
            string_bh = uix.VBox('Parent', main_vb, ...
                'Padding', 10, ...
                'BackgroundColor', GUI_Chalmers.BG_COLOR);
            
            txt = this.insertText(string_bh, {'Use the following string for the station locations:'},  9, [], 'left');
            txt.BackgroundColor = logo_GUI.BG_COLOR;
            
            %Core_UI.insertEmpty(string_bh);
            
            j_chalmers = this.insertChalmersBox(string_bh);
            j_chalmers.setText(this.getChalmersString);
            
            txt = this.insertText(string_bh, {'Save the ocean loading parameter values into a BLQ file and pass it to goGPS'},  9, [], 'left');
            txt.BackgroundColor = logo_GUI.BG_COLOR;
            
            % IMPORTANT NOTE:
            % If you are a contributor and your name is not in this list feel free to add your name
            
            % Manage dimension -------------------------------------------------------------------------------------------
            main_vb.Heights = [100 -1];
            string_bh.Heights = [23 -1 23];
            
            this.win.Visible = 'on';
        end
    end
    
    %% METHODS INSERT
    % ==================================================================================================================================================
    methods (Static)
        function str = getChalmersString(mode)
            if nargin == 0 || isempty(mode)
                mode = 'rinex';
            end
            
            str = sprintf('//------------------------------------------------------------------------');
            str = sprintf('%s\n%s', str, '// Station list for ocean loading computation');
            str = sprintf('%s\n%s',str, '//------------------------------------------------------------------------');
            
            new_only = true;
            if new_only
                station_code = GUI_Chalmers.getLoadingStationCode();
            else
                has_blq = false;
            end
            
            switch mode
                case {'rinex', 'RINEX'}
                    rec_path = Core.getState.getRecPath();
                    for r = 1 : numel(rec_path)
                        file_list = rec_path{r};
                        i = 0;
                        has_no_coo = true;
                        while (i < numel(file_list) && has_no_coo)
                            i = i + 1;
                            fr = File_Rinex(file_list(i), 100, true);
                            has_no_coo = isempty(fr.coo.getXYZ) || all(fr.coo.getXYZ == 0);
                        end
                        if fr.isValid()
                            name = fr.marker_name{1};
                            name = upper(name(1:min(4, numel(name))));
                            
                            if new_only
                                has_blq = ~isempty(find(station_code == Core_Utils.code4Char2Num(name), 1, 'first'));
                            end
                            if ~has_blq
                                xyz = median(fr.coo.getXYZ,1,'omitnan');
                                if ~isempty(xyz)
                                    str = sprintf('%s\n%s', str, sprintf('%-24s %16.4f%16.4f%16.4f', name, xyz(1), xyz(2),xyz(3)));
                                end
                            end
                        end
                    end
                case {'rec'}
                    core = Core.getCurrentCore;
                    sta_list = core.rec;
                    for r = 1 : size(sta_list, 2)
                        rec = sta_list(~sta_list(:,r).isEmpty, r);
                        if new_only
                            has_blq = ~isempty(find(station_code == Core_Utils.code4Char2Num(upper(rec.getMarkerName4Ch)), 1, 'first'));
                        end
                        if ~has_blq
                            if ~isempty(rec)
                                xyz = rec.out.getMedianPosXYZ();
                                if isempty(xyz)
                                    xyz = rec.work.getMedianPosXYZ();
                                end
                                str = sprintf('%s\n%s', str, sprintf('%-24s %16.4f%16.4f%16.4f', rec(1).getMarkerName4Ch, xyz(1), xyz(2),xyz(3)));
                            end
                        end
                    end
            end
            str = sprintf('%s\n%s', str,  '//------------------------------------------------------------------------');
        end
        
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
                'BackgroundColor', GUI_Chalmers.BG_COLOR);
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
                'BackgroundColor', GUI_Chalmers.BG_COLOR);
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
            bg_color = num2cell(GUI_Chalmers.BG_COLOR);
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
        
        function j_chalmers = insertChalmersBox(container)
            j_chalmers = com.mathworks.widgets.SyntaxTextPane;
            %codeType = j_chalmers.M_MIME_TYPE;  % j_chalmers.contentType='text/m-MATLAB'
            %j_chalmers.setContentType(codeType);
            str = '// Chalmers station list';
            j_chalmers.setText(str);
            % Create the ScrollPanel containing the widget
            j_scroll_settings = com.mathworks.mwswing.MJScrollPane(j_chalmers);
            % Inject edit box with the Java Scroll Pane into the main_window
            %% DEPRECATE!!!
            warning off
            [panel_j, panel_h] = javacomponent(j_scroll_settings, [1 1 1 1], container);
            warning on
            j_chalmers.setEditable(0);
        end
        
        function station_code = getLoadingStationCode(this)
            % Get numeric 4char marker names present in the ocean loading file
            %
            % SYNTAX
            %   station_code = this.getLoadingStationCode();
            
            ocean_file = Core.getState.getOceanFile;
            %ocean_file = '/Users/Andrea/Repositories/goGPS_MATLAB/data/project/Japan/station/ocean/ocean_loading.blq';
            
            % open RINEX observation file
            fid = fopen(ocean_file,'rt');
            if fid > 0
                txt = fread(fid,'*char')';
                % try to see if carriage return is present in the file (Windows stupid standard)
                % On Windows file lines ends with char(13) char(10)
                % instead of just using char(10)
                if ~isempty(find(txt(1:min(1000,numel(txt))) == 13, 1, 'first'))
                    has_cr = true;  % The file has carriage return - I hate you Bill!
                else
                    has_cr = false;  % The file is UNIX standard
                end
                % txt = txt(txt ~= 13);  % remove carriage return - I hate you Bill!
                fclose(fid);
                
                % get new line separators                
                nl = regexp(txt, '\n')';
                if isempty(nl)
                    station_code = [];
                else
                    if nl(end) <  (numel(txt) - double(has_cr))
                        nl = [nl; numel(txt)];
                    end
                    lim = [[1; nl(1 : end - 1) + 1] (nl - 1 - double(has_cr))];
                    lim = [lim (lim(:,2) - lim(:,1) + 1)];
                    while lim(end,3) < 3
                        lim(end,:) = [];
                    end
                    
                    % removing empty lines at end of file
                    min_line_width = 6;
                    lim(lim(1:end,3) < min_line_width,:) = [];
                    
                    % removing commented lines
                    lim(txt(lim(:,1)) == '$', :) = [];
                    
                    %txt(repmat(lim(:,1), 1, 4) + repmat(2 : 5, size(lim, 1), 1))
                    
                    % Convert 4ch marker to num
                    station_code = Core_Utils.code4Char2Num(upper(txt(repmat(lim(:,1), 1, 4) + repmat(2 : 5, size(lim, 1), 1))));
                end
            else
                station_code = [];
            end
        end
    end
    %% METHODS getters
    % ==================================================================================================================================================
    methods
    end
    
    %% METHODS EVENTS
    % ==================================================================================================================================================
    methods (Access = public)
    end
end
