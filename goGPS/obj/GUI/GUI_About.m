%   CLASS GUI_About
% =========================================================================
%
% DESCRIPTION
%   class to manages the about window of goGPSz
%
% EXAMPLE
%   ui = GUI_About.getInstance();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Core_UI


%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b6
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

classdef GUI_About < handle
    
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
        function this = GUI_About(w_main)
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
            
            win = figure( 'Name', 'About', ...
                'Visible', 'on', ...
                'DockControls', 'off', ...
                'MenuBar', 'none', ...
                'ToolBar', 'none', ...
                'NumberTitle', 'off', ...
                'Position', [0 0 680 540], ...
                'Resize', 'off');
            
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
            
            logo_GUI_About.BG_COLOR = Core_UI.DARK_GREY_BG;
            left_tbv = uix.VBox('Parent', top_bh, ...
                'BackgroundColor', logo_GUI_About.BG_COLOR, ...
                'Padding', 5);
            
            % Logo/title box -------------------------------------------------------------------------------------------
            
            logo_g = uix.Grid('Parent', left_tbv, ...
                'Padding', 5, ...
                'BackgroundColor', logo_GUI_About.BG_COLOR);
            
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
                        
            Core_UI.insertEmpty(left_tbv, logo_GUI_About.BG_COLOR);
            left_tbv.Heights = [82 -1];
            
            % Title Panel -----------------------------------------------------------------------------------------------
            right_tvb = uix.VBox('Parent', top_bh, ...
                'Padding', 5, ...
                'BackgroundColor', logo_GUI_About.BG_COLOR);

            top_bh.Widths = [106 -1];
            
            title = uix.HBox('Parent', right_tvb, ...
                'BackgroundColor', logo_GUI_About.BG_COLOR);
            
            txt = this.insertBoldText(title, 'goGPS', 10, Core_UI.LBLUE, 'left');
            txt.BackgroundColor = logo_GUI_About.BG_COLOR;
            title_l = uix.VBox('Parent', title, 'BackgroundColor', GUI_About.BG_COLOR);
            title.Widths = [54 -1];
            Core_UI.insertEmpty(title_l, logo_GUI_About.BG_COLOR)
            txt = this.insertBoldText(title_l, ['- software V' Core.GO_GPS_VERSION], 8, [], 'left');
            txt.BackgroundColor = logo_GUI_About.BG_COLOR;
            title_l.Heights = [2, -1];
            
            % Disclaimer Panel -----------------------------------------------------------------------------------------------
            Core_UI.insertEmpty(right_tvb, logo_GUI_About.BG_COLOR)
            txt = this.insertText(right_tvb, {['This release is loosely based on the original goGPS software, most of the code has '...
                'been rewritten for better performances, automation, flexibility, but with the primary goal of quasi static processing. ' ...
                'Legacy version can still be found at:'], ...
                'https://github.com/goGPS-Project/goGPS_MATLAB/tree/legacy'}, 7, [], 'left');
            txt.BackgroundColor = logo_GUI_About.BG_COLOR;
            
            right_tvb.Heights = [20 3 -1];
            
            credits_bh = uix.HBox('Parent', main_vb, ...
                'Padding', 5, ...
                'BackgroundColor', GUI_About.BG_COLOR);
            credits_legacy = uix.VBox('Parent', credits_bh, ...
                'BackgroundColor', GUI_About.BG_COLOR);
            credits_new = uix.VBox('Parent', credits_bh, ...
                'BackgroundColor', GUI_About.BG_COLOR);
            
            % Old Credits -----------------------------------------------------------------------------------------------
            
            this.insertBoldText(credits_legacy, 'Legacy version', 9);
            bar_cont = uix.HBox('Parent', credits_legacy, ...
                'BackgroundColor', GUI_About.BG_COLOR);
            
            Core_UI.insertEmpty(bar_cont, GUI_About.BG_COLOR)
            Core_UI.insertHBarDark(bar_cont);
            Core_UI.insertEmpty(bar_cont, GUI_About.BG_COLOR)
            bar_cont.Widths = [-1 220 -1];

            this.insertBoldText(credits_legacy, 'Developers of the original core', 8);
            this.insertText(credits_legacy, {'Mirko Reguzzoni', 'Eugenio Realini'}, 7);
            Core_UI.insertEmpty(credits_legacy, GUI_About.BG_COLOR)
            
            this.insertBoldText(credits_legacy, 'Founding fathers', 8);
            this.insertText(credits_legacy, {'Alessio Dominioni', ...
                'Roberto Terruzzi', ...
                'Marco Colla', ...
                'Marcello Palmulli', ...
                'Michele Tettamanti'}, 7);
            
            Core_UI.insertEmpty(credits_legacy, GUI_About.BG_COLOR)
            this.insertBoldText(credits_legacy, 'Contributors (in temporal order)', 8);
            this.insertText(credits_legacy, {'Sara Lucca', ...
                'Lisa Pertusini', ...
                'Daniele Sampietro', ...
                'Luana Valentini', ...
                'Alemu Befkadu Nigussie', ...
                'Ivan Reguzzoni', ...
                'Stefano Caldera', ...
                'Damiano Triglione', ...
                'Giuliano Sirioni', ...
                'Antonio Herrera Olmo', ...
                'Andrea Gatti', ...
                'Hendy F. Suhandri', ...
                'Andrea Nardo', ...
                'Daisuke Yoshida', ...
                'Marco Negretti', ...
                'Yuya Iwaki', ...
                'Louis Osei-Poku', ...
                'Paiyuan Zhou'}, 7);
            
            credits_legacy.Heights = [22 5 18 14*2 5 18 14*5 5 18 14*17];
            
            % New Credits -----------------------------------------------------------------------------------------------

            this.insertBoldText(credits_new, 'Object oriented new version', 9);
            bar_cont = uix.HBox('Parent', credits_new, ...
                'BackgroundColor', GUI_About.BG_COLOR);
            
            Core_UI.insertEmpty(bar_cont, GUI_About.BG_COLOR)
            Core_UI.insertHBarDark(bar_cont);
            Core_UI.insertEmpty(bar_cont, GUI_About.BG_COLOR)
            bar_cont.Widths = [-1 220 -1];

            this.insertBoldText(credits_new, 'Developers of the new goGPS core', 8);
            this.insertText(credits_new, {'Andrea Gatti', 'Giulio Tagliaferro'}, 7);
            Core_UI.insertEmpty(credits_new, GUI_About.BG_COLOR)

            this.insertBoldText(credits_new, 'Other Contributors', 8);
            this.insertText(credits_new, {'Eugenio Realini', ...
                'Stefano Caldera', ...
                'Lisa Pertusini', ...
                'Andreas Krietemeyer', ...
                'Alice Bonfiglio', ...
                'Stefano Barindelli',...
                'Alessandra Mascitelli', ...
                }, 7);
            Core_UI.insertEmpty(credits_new, GUI_About.BG_COLOR)

            % IMPORTANT NOTE:
            % If you are a contributor and your name is not in this list feel free to add your name
            
            this.insertBoldText(credits_new, 'Acknoledgments', 8);
            this.insertText(credits_new, {'This release is funded by', 'GReD srl (www.g-red.eu)'}, 7);
            txt = this.insertText(credits_new, {'I dedicate my work on this software to my father, Tiziano Gatti, ', 'who taught me, among countless things, the importance of putting passion in everything you choose to do. I will never forget you,'}, 7, [], 'right'); txt.FontAngle = 'italic';
            txt = this.insertText(credits_new, {'Andrea Gatti'}, 7, [20 183 242]/255, 'right'); txt.FontAngle = 'italic';
            
            credits_new.Heights = [22 5 18 14*2 5 18 14*7 5 18 -1 14*3 28];
            
            

            % Manage dimension -------------------------------------------------------------------------------------------
            
            main_vb.Heights = [84 -1];
            
            %session_height = sum(left_bv.Children(2).Children(1).Heights);
            
            this.win.Visible = 'on';
            %uiwait(this.win);
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
                'BackgroundColor', GUI_About.BG_COLOR);
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
                'BackgroundColor', GUI_About.BG_COLOR);
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
