%   CLASS GUI_Inspector
% =========================================================================
%
% DESCRIPTION
%   class to manages the about window of goGPSz
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
        end
        
        function openGUI(this)
            % Main Window ----------------------------------------------------------------------------------------------
            
            win = figure( 'Name', 'goGPS inspector', ...
                'Visible', 'on', ...
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
            
            % Logo/title box -------------------------------------------------------------------------------------------
            
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
            
            % Title Panel -----------------------------------------------------------------------------------------------
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
            
            % Disclaimer Panel -----------------------------------------------------------------------------------------------
            Core_UI.insertEmpty(right_tvb, logo_GUI_Inspector.BG_COLOR)
            txt = this.insertText(right_tvb, {'A GNSS processing software powered by GReD'}, 9, [], 'left');
            txt.BackgroundColor = logo_GUI_Inspector.BG_COLOR;
            right_tvb.Heights = [25 3 -1];
            
            
            Core_UI.insertEmpty(main_vb, Core_UI.DARKER_GREY_BG)
            
            % Logging Panel --------------------------------------------------------------------------------------------------
            
            main_hb = uix.HBox('Parent', main_vb, ...
                'Padding', 5, ...
                'BackgroundColor', logo_GUI_Inspector.BG_COLOR);
            main_vb.Heights = [84 5 -1];
            
            this.insertRecList(main_hb);
            Core_UI.insertEmpty(main_hb, Core_UI.DARK_GREY_BG)
            Core_UI.insertEmpty(main_hb, Core_UI.DARKER_GREY_BG)
            main_hb.Widths = [175 5 -1];
            
            % Manage dimension -------------------------------------------------------------------------------------------
            
                        
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
                'Callback', @this.selectAll); %#ok<NASGU>
            clear_but = uicontrol( 'Parent', list_but, ...
                'String', 'Clear All', ...
                'Callback', @this.unselectAll); %#ok<NASGU>           
                        
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
            for r = 1 : numel(core.rec)
                name = core.rec(r).getMarkerName4Ch;                                
                
                this.rec_tbl.Data{r,1} = true;
                this.rec_tbl.Data{r,2} = sprintf('%s style="font-weight: bold; font-size: 9px; color: #6666FF; ">%d', '<html><tr><td width=9999 align=center ', r);
                this.rec_tbl.Data{r,3} = sprintf('%s style="font-weight: bold; font-size: 9px; color: #6666FF; ">%s', '<html><tr><td width=9999 align=center ', upper(name(1:4)));
                this.rec_tbl.Data{r,4} = sprintf('%s style="font-weight: bold; font-size: 9px; color: #444444; ">%s', '<html><tr><td width=9999 align=center ', core.rec(r).getActiveSys());
            end
        end
        
        function selectAll(this, caller, event)
            % Select all the receivers
            for r = 1 : size(this.rec_tbl.Data, 1)
                this.rec_tbl.Data{r,1} = true;
            end
        end
        
        function unselectAll(this, caller, event)
            % Select all the receivers
            for r = 1 : size(this.rec_tbl.Data, 1)
                this.rec_tbl.Data{r,1} = false;
            end
        end

    end
    
    %% METHODS EVENTS
    % ==================================================================================================================================================
    methods (Access = public)         
    end
end
