%   CLASS Core_UI
% =========================================================================
%
% DESCRIPTION
%   class to manages the user interface of goGPS
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

classdef Core_UI < Logos
    
    properties (Constant)        
        DEFAULT_MODE = 'dark';
        DEFAULT_EXPORT_MODE = 'light';
        
        FONT_SIZE_CONVERSION_LNX = 1;
        FONT_SIZE_CONVERSION_MAC = 1.45;
        FONT_SIZE_CONVERSION_WIN = 1;
        
        LIGHT_GREY_BG = 0.85 * ones(3, 1);
        DARK_GREY_BG = [0.2431372549 0.2470588235 0.2509803922];
        DARKER_GREY_BG = [0.1529411765 0.1529411765 0.1568627451];
        WHITE = ones(3, 1);
        BLACK = zeros(3, 1);
        GREY = [0.6 0.6 0.6];
        BLUE = [0 0 1];
        RED  = [1 0 0];
        GREEN = [0 1 0];
        ORANGE = [1 0.6 0];
        LBLUE = [0 163 222]/255;
        
        COLOR_ORDER = [ ...
            0     0.447 0.741;
            0.85  0.325 0.098;
            0.929 0.694 0.125;
            0.494 0.184 0.556;
            0.466 0.674 0.188;
            0.301 0.745 0.933;
            0.635 0.078 0.184];
        
        CMAP_51 = [      0                         0                         0
                         1         0.901960784313726                         1
                         1         0.843137254901961                         1
                         1         0.784313725490196                         1
          0.96078431372549         0.705882352941177          0.96078431372549
         0.901960784313726         0.627450980392157         0.901960784313726
         0.862745098039216         0.549019607843137         0.862745098039216
         0.803921568627451         0.470588235294118         0.803921568627451
         0.745098039215686         0.352941176470588         0.745098039215686
         0.686274509803922         0.235294117647059         0.686274509803922
         0.588235294117647         0.176470588235294         0.588235294117647
         0.549019607843137                         0         0.549019607843137
         0.549019607843137                         0         0.431372549019608
         0.549019607843137                         0         0.274509803921569
         0.549019607843137                         0                         0
         0.627450980392157                         0                         0
         0.725490196078431                         0                         0
         0.784313725490196                         0                         0
         0.843137254901961                         0                         0
         0.894117647058824                         0                         0
         0.941176470588235                         0                         0
                         1                         0                         0
                         1         0.156862745098039                         0
                         1         0.274509803921569                         0
                         1         0.392156862745098                         0
                         1         0.509803921568627                         0
                         1         0.627450980392157                         0
                         1         0.725490196078431                         0
                         1         0.823529411764706                         0
                         1         0.901960784313726                         0
                         1                         1                         0
         0.949019607843137                         1                         0
         0.901960784313726                         1                         0
         0.823529411764706                         1                         0
         0.549019607843137                         1                         0
                         0                         1                         0
                         0         0.901960784313726                         0
                         0         0.823529411764706                         0
                         0         0.745098039215686                         0
                         0         0.647058823529412                         0
                         0         0.588235294117647         0.117647058823529
                         0         0.529411764705882         0.274509803921569
                         0         0.470588235294118         0.470588235294118
                         0         0.411764705882353         0.647058823529412
                         0         0.352941176470588         0.823529411764706
                         0         0.392156862745098         0.980392156862745
                         0         0.509803921568627         0.980392156862745
                         0         0.705882352941177         0.980392156862745
                         0         0.862745098039216         0.980392156862745
         0.313725490196078         0.941176470588235         0.980392156862745
         0.705882352941177         0.980392156862745         0.980392156862745];
     
        LINE_HEIGHT = 23
    end
    %% PROPERTIES SINGLETON POINTERS
    % ==================================================================================================================================================
    properties % Utility Pointers to Singletons
        state
        w_bar
    end
    %% PROPERTIES GUI
    % ==================================================================================================================================================
    properties
        main        % Handle of the main
    end    
    %% METHOD CREATOR
    % ==================================================================================================================================================
    methods (Static, Access = private)
        % Concrete implementation.  See Singleton superclass.
        function this = Core_UI()
            % Core object creator
        end
    end
    %% METHODS UI
    % ==================================================================================================================================================
    methods (Static, Access = public)
        function this = getInstance()
            % Get the persistent instance of the class
            persistent unique_instance_core_ui__
            
            if isempty(unique_instance_core_ui__)
                this = Core_UI();
                unique_instance_core_ui__ = this;
            else
                this = unique_instance_core_ui__;
            end
            this.init();
        end
        
        function showTextHeader()
            % Display as a text the Header containing goGPS ASCII title and version
            %
            % SYNTAX:
            %   Core_UI.showTextHeader();
            log = Core.getLogger();
            if log.getColorMode()
                %cprintf([241 160 38]/255,'\n               ___ ___ ___\n     __ _ ___ / __| _ | __|\n    / _` / _ \\ (_ |  _|__ \\\n    \\__, \\___/\\___|_| |___/\n    |___/                    '); cprintf('text','v '); cprintf('text', Core.GO_GPS_VERSION); fprintf('\n');
                %cprintf([0 163 222]/255,'\n               ___ ___ ___\n     __ _ ___ / __| _ | __|\n    / _` / _ \\ (_ |  _|__ \\\n    \\__, \\___/\\___|_| |___/\n    |___/                    '); cprintf('text','v '); cprintf('text', Core.GO_GPS_VERSION); fprintf('\n');
                cprintf([0 163 222]/255,'\n               ___ ___ ___\n     __ _ ___ / __| _ | __|\n    / _` / _ \\ (_ |  _|__ \\\n    \\__, \\___/\\___|_| |___/\n    |___/   '); cprintf('text',[iif(Core.isGReD, 'GReD', 'OPEN') ' EDITION     v ']); cprintf('text', Core.GO_GPS_VERSION); fprintf('\n');
                fprintf('\n--------------------------------------------------------------------------\n');
                fprintf('    GNSS data processing powered by GReD\n');
            else
                %fprintf('\n               ___ ___ ___\n     __ _ ___ / __| _ | __|\n    / _` / _ \\ (_ |  _|__ \\\n    \\__, \\___/\\___|_| |___/\n    |___/                    v %s\n', Core.GO_GPS_VERSION);
                fprintf(['\n               ___ ___ ___\n     __ _ ___ / __| _ | __|\n    / _` / _ \\ (_ |  _|__ \\\n    \\__, \\___/\\___|_| |___/\n    |___/   ' iif(Core.isGReD, 'GReD', 'OPEN') ' EDITION     v %s\n'], Core.GO_GPS_VERSION);
                fprintf('\n');
                fprintf('\n--------------------------------------------------------------------------\n');
                fprintf('    GNSS data processing powered by GReD\n');
                fprintf('\n--------------------------------------------------------------------------\n');
            end
            % log.addWarning('This is goGPS nightly build\nSome parts (or all of it) could not work properly\nUse at your own risk!');
            log.addWarning('This is a goGPS beta build, use it at your own risk!\nSome parts of it could not work properly\nPlease open a new issue on github if you found any bug');
            log.simpleSeparator();
            fprintf('\n');
        end
        
        function str_out = getTextHeader()
            % Export as a string the Header containing goGPS ASCII title and version
            %
            % SYNTAX:
            %   str_out = Core_UI.getTextHeader();
            
            str_out = sprintf('\n               ___ ___ ___\n     __ _ ___ / __| _ | __|\n    / _` / _ \\ (_ |  _|__ \\\n    \\__, \\___/\\___|_| |___/\n    |___/                    v %s\n', Core.GO_GPS_VERSION);
            str_out = sprintf('%s\n--------------------------------------------------------------------------\n',str_out);
        end        
    end
    
    %% METHODS FIGURE MODIFIER
    methods (Static, Access = public)
        function addBeautifyMenu(fig_handle)
            % Add a menu Aspect to figure
            % with 2 sub menu for changing figure aspect 
            % to Light or Dark Mode
            %
            % SYNTAX
            %   Core_Utils.addBeautifyMenu(fig_handle)
            
            if nargin == 0 || isempty(fig_handle)
                fig_handle = gcf;
            end
            
            m = findall(fig_handle.Children, 'Type', 'uimenu', 'Label', 'Aspect');
            if ~isempty(m)
                m = m(1);
            else
                m = uimenu(fig_handle, 'Label', 'Aspect');
            end
            mitem = findall(m.Children, 'Type', 'uimenu', 'Label', '&Dark mode');
            if ~isempty(mitem)
                % Item already present
                %    mitem = mitem(1);
            else
                mitem = uimenu(m,'Label','&Dark mode');
                mitem.Accelerator = 'D';
                mitem.Callback = @toDark;
            end
            
            mitem = findall(m.Children, 'Type', 'uimenu', 'Label', '&Light mode');
            if ~isempty(mitem)
                % Item already present
                %    mitem = mitem(1);
            else
                mitem = uimenu(m,'Label','&Light mode');
                mitem.Accelerator = 'L';
                mitem.Callback = @toLight;
            end
            
            function toDark(src, event)
                Core_UI.beautifyFig(src.Parent.Parent, 'dark');
            end
            
            function toLight(src, event)
                Core_UI.beautifyFig(src.Parent.Parent, 'light');
            end
        end
        
        function beautifyFig(fig_handle, color_mode)
            % Change font size / type colors of a figure
            % 
            % INPUT:
            %   fig_handle  figure handler (e.g. gcf)
            %   color_mode  two modes are supported:
            %               'light'     classic Light mode
            %               'dark'      dark mode
            %
            % SYNTAX:
            %   Core_UI.beautifyFig(fig_handle, color_mode)
            FONT = 'Open Sans';
            %FONT = 'Helvetica';
            if nargin == 0 || isempty(fig_handle)
                fig_handle = gcf;
            end
            if nargin < 2 || isempty(color_mode)
                color_mode = Core_UI.DEFAULT_MODE;
            end
            ax_list = findall(fig_handle,'type','axes');
            for ax = ax_list(:)'
                ax.FontName = 'Open Sans';
                ax.FontWeight = 'bold';
            end
            set(fig_handle, ...
                'DefaultFigureColor', 'w', ...
                'DefaultAxesLineWidth', 0.5, ...
                'DefaultAxesXColor', 'k', ...
                'DefaultAxesYColor', 'k', ...
                'DefaultAxesFontUnits', 'points', ...
                'DefaultAxesFontSize', Core_UI.getFontSize(8), ...
                'DefaultAxesFontName', FONT, ...
                'DefaultLineLineWidth', 1, ...
                'DefaultTextFontUnits', 'Points', ...
                'DefaultTextFontSize', Core_UI.getFontSize(13), ...
                'DefaultTextFontName', FONT, ...
                'DefaultTextFontWeight', 'bold', ...
                'DefaultAxesBox', 'off', ...
                'DefaultAxesTickLength', [0.02 0.025]);
            
            set(fig_handle, 'DefaultAxesTickDir', 'out');
            set(fig_handle, 'DefaultAxesTickDirMode', 'manual');
            
            if strcmp(color_mode, 'dark') % ------------------------------------------------------------------- DARK
                fig_handle.Color = [0.15, 0.15 0.15];
                box_list = findall(fig_handle, 'Type', 'uicontainer');
                for box = box_list(:)'
                    if ~ischar(box.BackgroundColor) && all(box.BackgroundColor > 0.5)
                        box.BackgroundColor = 1 - box.BackgroundColor;
                    end
                end
                ui_list = findall(fig_handle, 'Type', 'uicontrol');
                for ui = ui_list(:)'
                    if ~ischar(ui.BackgroundColor) && all(ui.BackgroundColor > 0.5)
                        ui.ForegroundColor = 1 - ui.ForegroundColor;
                        ui.FontName = FONT;
                        ui.FontSize = iif(ui.FontSize == 12, Core_UI.getFontSize(12), Core_UI.getFontSize(13));
                        ui.BackgroundColor = 1 - ui.BackgroundColor;
                    end
                end
                
                cb_list = findall(fig_handle, 'Type', 'colorbar');
                for cb = cb_list(:)'
                    cb.Color = [0.9 0.9 0.9];
                    cb.FontName = FONT;
                    cb.FontSize = iif(cb.FontSize == 12, Core_UI.getFontSize(12), Core_UI.getFontSize(13));
                    cbt_list = findall(cb.UserData, 'Type', 'text');
                    for cbt = cbt_list(:)'
                        cbt.Color = [0.9 0.9 0.9];
                        cbt.FontName = FONT;
                        cbt.FontSize = Core_UI.getFontSize(12);
                    end
                end
                ax_list = findall(fig_handle,'type','axes');
                for ax = ax_list(:)'
                    ax.Color = [0.2 0.2 0.2];
                    ax.Title.Color = [1 1 1];
                    ax.XLabel.Color = [0.9 0.9 0.9];
                    ax.YLabel.Color = [0.9 0.9 0.9];
                    ax.ZLabel.Color = [0.9 0.9 0.9];
                    ax.XColor = [0.9 0.9 0.9];
                    ax.YColor = [0.9 0.9 0.9];
                    ax.ZColor = [0.9 0.9 0.9];
                    ax.FontName = FONT;
                    ax.FontSize = Core_UI.getFontSize(11);
                    text_label = findall(ax, 'Type', 'text');
                    for txt = text_label(:)'
                        txt.Color = [0.9 0.9 0.9];
                        if ~ischar(txt.BackgroundColor) && all(txt.BackgroundColor > 0.5)
                            txt.BackgroundColor = 1 - txt.BackgroundColor;
                        end
                        txt.FontName = FONT;
                        txt.FontSize = iif(txt.FontSize == 10, Core_UI.getFontSize(9), Core_UI.getFontSize(12));
                    end
                end
                text_label = findall(gcf,'Tag', 'm_grid_xticklabel');
                for txt = text_label(:)'
                    txt.Color = [0.9 0.9 0.9];
                    txt.FontName = FONT;
                    txt.FontSize = iif(txt.FontSize == 12, 13, 15);
                end
                text_label = findall(gcf,'Tag', 'm_grid_yticklabel');
                for txt = text_label(:)'
                    txt.Color = [0.9 0.9 0.9];
                    txt.FontName = FONT;
                    txt.FontSize = iif(txt.FontSize == 12, 13, 15);
                end
                text_label = findall(gcf,'Tag', 'm_ruler_label');
                for txt = text_label(:)'
                    txt.Color = [0.9 0.9 0.9];
                    txt.FontName = FONT;
                    txt.FontSize = iif(txt.FontSize == 12, 13, 15);
                end
                legend = findall(gcf, 'type', 'legend');
                for lg = legend(:)'
                    lg.FontName = FONT;
                    lg.Color = [1 1 1];
                    lg.Title.Color = [0 0 0];
                    lg.TextColor = 1-[0.9 0.9 0.9];
                    lg.EdgeColor = [0.5 0.5 0.5];
                end
                
            elseif strcmp(color_mode, 'light') % ------------------------------------------------------------------- LIGHT
                fig_handle.Color = [1 1 1];
                box_list = findall(fig_handle, 'Type', 'uicontainer');
                for box = box_list(:)'
                    if ~ischar(box.BackgroundColor) && all(box.BackgroundColor < 0.5)
                        box.BackgroundColor = 1 - box.BackgroundColor;
                    end
                end
                ui_list = findall(fig_handle, 'Type', 'uicontrol');
                for ui = ui_list(:)'
                    if ~ischar(ui.BackgroundColor) && all(ui.BackgroundColor < 0.5)
                        ui.ForegroundColor = 1 - ui.ForegroundColor;
                        ui.BackgroundColor = 1 - ui.BackgroundColor;
                    end
                end
                
                cb_list = findall(fig_handle, 'Type', 'colorbar');
                for cb = cb_list(:)'
                    cb.Color = 1-[0.9 0.9 0.9];
                    cb.FontName = FONT;
                    cb.FontSize = iif(cb.FontSize == 12, Core_UI.getFontSize(12), Core_UI.getFontSize(13));
                    cbt_list = findall(cb.UserData, 'Type', 'text');
                    for cbt = cbt_list(:)'
                        cbt.Color = 1-[0.9 0.9 0.9];
                        cbt.FontName = FONT;
                        cbt.FontSize = Core_UI.getFontSize(12);
                    end
                end
                ax_list = findall(fig_handle,'type','axes');
                for ax = ax_list(:)'
                    ax.Color = [1 1 1];
                    ax.Title.Color = 1-[1 1 1];
                    ax.XLabel.Color = 1-[0.9 0.9 0.9];
                    ax.YLabel.Color = 1-[0.9 0.9 0.9];
                    ax.ZLabel.Color = 1-[0.9 0.9 0.9];
                    ax.XColor = 1-[0.9 0.9 0.9];
                    ax.YColor = 1-[0.9 0.9 0.9];
                    ax.ZColor = 1-[0.9 0.9 0.9];
                    ax.FontName = FONT;
                    ax.FontSize = Core_UI.getFontSize(11);
                    text_label = findall(ax, 'Type', 'text');
                    for txt = text_label(:)'
                        txt.Color = 1-[0.9 0.9 0.9];
                        if ~ischar(txt.BackgroundColor) && all(txt.BackgroundColor < 0.5)
                            txt.BackgroundColor = 1 - txt.BackgroundColor;
                        end
                        txt.FontName = FONT;
                        txt.FontSize = iif(txt.FontSize == 10, Core_UI.getFontSize(9), Core_UI.getFontSize(12));
                    end
                end                
                text_label = findall(gcf,'Tag', 'm_grid_xticklabel');
                for txt = text_label(:)'
                    txt.Color = 1-[0.9 0.9 0.9];
                    txt.FontName = FONT;
                    txt.FontSize = iif(txt.FontSize == 12, Core_UI.getFontSize(12), Core_UI.getFontSize(14));
                end
                text_label = findall(gcf,'Tag', 'm_grid_yticklabel');
                for txt = text_label(:)'
                    txt.Color = 1-[0.9 0.9 0.9];
                    txt.FontName = FONT;
                    txt.FontSize = iif(txt.FontSize == 12, 13, 15);
                end
                text_label = findall(gcf,'Tag', 'm_ruler_label');
                for txt = text_label(:)'
                    txt.Color = 1-[0.9 0.9 0.9];
                    txt.FontName = FONT;
                    txt.FontSize = iif(txt.FontSize == 12, Core_UI.getFontSize(12), Core_UI.getFontSize(14));
                end
                legend = findall(gcf, 'type', 'legend');
                for lg = legend(:)'
                    lg.FontName = FONT;
                    lg.Color = [1 1 1];
                    lg.Title.Color = [0 0 0];
                    lg.TextColor = 1-[0.9 0.9 0.9];
                    lg.EdgeColor = [0.5 0.5 0.5];
                end                
            end
            
            % ResizeFig
            unit = fig_handle.Units;
            fig_handle.Units = 'pixels';
            drawnow
            if isprop(fig_handle,'InnerPosition')
                flag_small_fig = fig_handle.InnerPosition(3) <= (600) && fig_handle.InnerPosition(4) <= (450);
            else
                flag_small_fig = fig_handle.Position(3) <= (600) && fig_handle.Position(4) <= (450);
            end
            
            if ~strcmp(fig_handle.WindowStyle, 'docked') && flag_small_fig
                if isprop(fig_handle,'InnerPosition')
                    fig_handle.InnerPosition([3, 4]) = [1020 700];
                else
                    fig_handle.Position([3, 4]) = [1020 700];
                end
                if isunix && not(ismac())
                    fig_handle.Position(1) = round((fig_handle.Parent.ScreenSize(3) - fig_handle.Position(3)) / 2);
                    fig_handle.Position(2) = round((fig_handle.Parent.ScreenSize(4) - fig_handle.Position(4)) / 2);
                else
                    fig_handle.OuterPosition(1) = round((fig_handle.Parent.ScreenSize(3) - fig_handle.OuterPosition(3)) / 2);
                    fig_handle.OuterPosition(2) = round((fig_handle.Parent.ScreenSize(4) - fig_handle.OuterPosition(4)) / 2);
                end
            end
            fig_handle.Units = unit;
        end
        
    end
    %% METHODS INIT
    % ==================================================================================================================================================
    methods
        function init(this)
            this.state = Core.getState();
            this.w_bar = Go_Wait_Bar.getInstance(100,'Init core GUI', Core.GUI_MODE);  % 0 means text, 1 means GUI, 5 both
        end
        
        function openGUI(this, flag_wait)
            if (nargin < 1) || isempty(flag_wait)
                flag_wait = false;
            end
            this.main = GUI_Edit_Settings.getInstance(flag_wait);
        end
    end
    %% METHODS INSERT
    % ==================================================================================================================================================
    methods (Static)
        
        function insertLogoGUI(container)
            % Logo
            logo_bg_color = Core_UI.LIGHT_GREY_BG;
            
            logo_g_border0 = uix.Grid('Parent', container, ...
                'Padding', 2, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            logo_g_border1 = uix.Grid('Parent', logo_g_border0, ...
                'Padding', 2, ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
            logo_g = uix.Grid('Parent', logo_g_border1, ...
                'Padding', 2, ...
                'BackgroundColor', logo_bg_color);
            
            Core_UI.insertEmpty(logo_g);
            logo_ax = axes( 'Parent', logo_g);
            if Core.isGReD()
                [logo, transparency] = Core_UI.getLogoGReDEdition();
            else
                [logo, transparency] = Core_UI.getLogoOpenEdition();
            end
            logo(repmat(sum(logo,3) == 0,1,1,3)) = 0;
            %logo = logo - 20;
            image(logo_ax, logo, 'AlphaData', transparency);
            logo_ax.XTickLabel = [];
            logo_ax.YTickLabel = [];
            axis(logo_ax, 'equal', 'tight');
            axis(logo_ax, 'off');
            logo_ax.Position([1,2]) = [size(logo,1) size(logo,2)];
            % Title and description
            version_txt = text(220, 85, [Core.GO_GPS_VERSION ' '], ...
                'FontName', 'verdana', ...
                'FontSize', Core_UI.getFontSize(5), ...
                'FontWeight', 'bold');
            
            logo_g.Widths = 184;
            logo_g.Heights = [4 78];
        end
        
        function insertLogoGUI_legacy(container)
            % Logo
            logo_bg_color = Core_UI.LIGHT_GREY_BG;
            
            logo_g_border0 = uix.Grid('Parent', container, ...
                'Padding', 2, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            logo_g_border1 = uix.Grid('Parent', logo_g_border0, ...
                'Padding', 2, ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
            logo_g = uix.Grid('Parent', logo_g_border1, ...
                'Padding', 5, ...
                'BackgroundColor', logo_bg_color);
            
            logo_ax = axes( 'Parent', logo_g);
            [logo, transparency] = Core_UI.getLogo();
            logo(repmat(sum(logo,3) == 0,1,1,3)) = 0;
            logo = logo - 20;
            image(logo_ax, logo, 'AlphaData', transparency);
            logo_ax.XTickLabel = [];
            logo_ax.YTickLabel = [];
            axis off;
            
            % Title and description
            descr_bv = uix.VBox('Parent', logo_g, ...
                'BackgroundColor', logo_bg_color);
            title_txt = uicontrol('Parent', descr_bv, ...
                'Style', 'Text', ...
                'HorizontalAlignment', 'center', ...
                'String', '  goGPS', ...
                'BackgroundColor', logo_bg_color, ...
                'ForegroundColor', 0 * ones(3, 1), ...
                'FontName', 'verdana', ...
                'FontAngle', 'italic', ...
                'FontSize', Core_UI.getFontSize(15), ...
                'FontWeight', 'bold');
            descr_txt = uicontrol('Parent', descr_bv, ...
                'Style', 'Text', ...
                'String', 'open source positioning', ...
                'BackgroundColor', logo_bg_color, ...
                'ForegroundColor', 0 * ones(3, 1), ...
                'FontName', 'arial', ...
                'FontSize', Core_UI.getFontSize(6), ...
                'FontWeight', 'bold');
            version_txt = uicontrol('Parent', descr_bv, ...
                'Style', 'Text', ...
                'HorizontalAlignment', 'right', ...
                'String', [Core.GO_GPS_VERSION ' '], ...
                'BackgroundColor', logo_bg_color, ...
                'ForegroundColor', 0 * ones(3, 1), ...
                'FontName', 'verdana', ...
                'FontSize', Core_UI.getFontSize(5), ...
                'FontWeight', 'bold');
            uicontrol('Parent', descr_bv, 'Style', 'Text', 'BackgroundColor', logo_bg_color);
            descr_bv.Heights = [-3.5 -2 -1.3 -0.5];
            
            logo_g.Widths = [64 -1];
            logo_g.Heights = 64;
        end
        
        % Single element insert ----------------------------------------------------------------
        
        function flag_handle = insertFlag(container, properties)
            % insert a flag element
            %
            % SYNTAX
            %   flag_handle = insertFlag(container)
            flag_handle = axes( 'Parent', container);
            
            if nargin == 2
                flag_handle.UserData = properties;
                Core_UI.checkFlag(flag_handle)
            else
                Core_UI.setFlagGrey(flag_handle);
            end
            flag_handle.XTickLabel = [];
            flag_handle.YTickLabel = [];                    
            axis(flag_handle, 'off');
        end            
           
        function chxbox_handle = insertCheckBoxCC(parent, title, but_tag, callback)
            chxbox_handle = uicontrol('Parent', parent,...
                'Style', 'checkbox',...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG, ...
                'String', title, ...
                'UserData', but_tag,...
                'Callback', callback ...
                );
        end
        
        function chxbox_handle = insertCheckBoxLight(container, title, state_field, callback)
            chxbox_handle = uicontrol('Parent', container,...
                'Style', 'checkbox',...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG, ...
                'FontSize', Core_UI.getFontSize(9), ...
                'String', title, ...
                'UserData', state_field,...
                'Callback', callback ...
                );
        end
        
        function chxbox_handle = insertCheckBoxDark(container, title, state_field, callback)
            chxbox_handle = uicontrol('Parent', container,...
                'Style', 'checkbox',...
                'BackgroundColor', Core_UI.DARK_GREY_BG, ...
                'FontSize', Core_UI.getFontSize(8), ...
                'String', title, ...
                'UserData', state_field,...
                'Callback', callback ...
                );
        end
        
        function insertEmpty(container, color)
            if nargin == 1
                color = Core_UI.LIGHT_GREY_BG;
            end
            uicontrol('Parent', container, 'Style', 'Text', 'BackgroundColor', color)
        end
        
        function panel_handle = insertPanelLight(container, title)
            panel_handle = uix.Panel('Parent', container, ...
                'Title', title, ...
                'FontSize', Core_UI.getFontSize(8), ...
                'FontWeight', 'bold', ...
                'Padding', 5, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
        end
        
        function insertHBarDark(container)
            bar_v = uix.VBox('Parent', container, ...
                'Padding', 0, ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
            Core_UI.insertEmpty(bar_v, Core_UI.DARK_GREY_BG);
            bar = uix.Panel( 'Parent', bar_v, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            Core_UI.insertEmpty(bar_v, Core_UI.DARK_GREY_BG);
            bar_v.Heights = [-1 1 -1];
        end
        
        function insertVBarDark(container)
            bar_h = uix.HBox('Parent', container, ...
                'Padding', 0, ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
            Core_UI.insertEmpty(bar_h, Core_UI.DARK_GREY_BG);
            bar = uix.Panel( 'Parent', bar_h, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            Core_UI.insertEmpty(bar_h, Core_UI.DARK_GREY_BG);
            bar_h.Widths = [-1 1 -1];
        end
        
        function insertHBarLight(container)
            bar_v = uix.VBox('Parent', container, ...
                'Padding', 0, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            Core_UI.insertEmpty(bar_v, Core_UI.LIGHT_GREY_BG);
            bar = uix.Panel( 'Parent', bar_v, ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
            Core_UI.insertEmpty(bar_v, Core_UI.LIGHT_GREY_BG);
            bar_v.Heights = [-1 2 -1];
        end
        
        function insertVBarLight(container)
            bar_h = uix.HBox('Parent', container, ...
                'Padding', 0, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            Core_UI.insertEmpty(bar_h, Core_UI.LIGHT_GREY_BG);
            bar = uix.Panel( 'Parent', bar_h, ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
            Core_UI.insertEmpty(bar_h, Core_UI.LIGHT_GREY_BG);
            bar_h.Widths = [-1 1 -1];
        end
        
        function date = insertDateSpinner(container, date_in, callback)
            % Initialize JIDE's usage within Matlab
            com.mathworks.mwswing.MJUtilities.initJIDE;
            
            % Display a DateChooserPanel
            date = com.jidesoft.combobox.DateSpinnerComboBox;
            
            [h_panel, h_container] = javacomponent(date,[10,10,140,20], container);
            
            set(h_panel, 'ItemStateChangedCallback', callback);
            
            date.setBackground(java.awt.Color(Core_UI.LIGHT_GREY_BG(1), Core_UI.LIGHT_GREY_BG(2), Core_UI.LIGHT_GREY_BG(3)));
            date.setDisabledBackground(java.awt.Color(Core_UI.LIGHT_GREY_BG(1),Core_UI.LIGHT_GREY_BG(2),Core_UI.LIGHT_GREY_BG(3)));
            date.setShowWeekNumbers(false);     % Java syntax
            dateFormat = java.text.SimpleDateFormat('dd-MM-yyyy');
            date.setFormat(dateFormat);
            date.setDate(java.util.Date(date_in))
        end
        
        function date_hour = insertDateSpinnerHour(container, time, call_back)
            date_hour = uix.HBox('Parent', container, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            Core_UI.insertDateSpinner(date_hour, time.toString('yyyy/mm/dd'), call_back);              
            box_handle = uix.Grid('Parent', date_hour, ...
                'Padding', 0, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            uicontrol('Parent', box_handle, ...
                'Style', 'Text', ...
                'String', 'hh:mm:ss', ...
                'ForegroundColor', Core_UI.BLACK, ...
                'HorizontalAlignment', 'center', ...
                'FontSize', Core_UI.getFontSize(9), ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            editable_handle = uicontrol('Parent', box_handle,...
                'Style', 'edit',...
                'FontSize', Core_UI.getFontSize(9), ...
                'Callback', call_back);
        end
        
        function [box_handle, pop_up_handle] = insertPopUpLight(parent, description, possible_values, property_name, callback, widths)
            if nargin < 6
                widths  = [-1 -1];
            end
            box_handle = uix.HBox('Parent', parent, ...
                'Padding', 0, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            uicontrol('Parent', box_handle, ...
                'Style', 'Text', ...
                'String', description, ...
                'ForegroundColor', Core_UI.BLACK, ...
                'HorizontalAlignment', 'left', ...
                'FontSize', Core_UI.getFontSize(9), ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            pop_up_handle = uicontrol('Parent', box_handle,...
                'Style', 'popup',...
                'UserData', property_name,...
                'String', possible_values,...
                'Callback', callback);
            box_handle.Widths = widths;
        end
        
        function [box_handle, pop_up_handle] = insertPopUpDark(parent, description, possible_values, property_name, callback, widths)
            if nargin < 6
                widths  = [-1 -1];
            end
            box_handle = uix.HBox('Parent', parent, ...
                'Padding', 0, ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
            uicontrol('Parent', box_handle, ...
                'Style', 'Text', ...
                'String', description, ...
                'ForegroundColor', Core_UI.WHITE, ...
                'HorizontalAlignment', 'left', ...
                'FontSize', Core_UI.getFontSize(9), ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
            pop_up_handle = uicontrol('Parent', box_handle,...
                'Style', 'popup',...
                'UserData', property_name,...
                'String', possible_values,...
                'Callback', callback);
            box_handle.Widths = widths;
        end
        
        function [box_handle, editable_handle, flag_handle] = insertFileBox(parent, description, property_name, callback, widths)
            if nargin < 5
                widths  = [25 -1 -1 25];
            end
            box_handle = uix.Grid('Parent', parent, ...
                'Padding', 0, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            if numel(widths) == 4
                flag_handle = Core_UI.insertFlag(box_handle, {property_name});
            end
            uicontrol('Parent', box_handle, ...
                'Style', 'Text', ...
                'String', description, ...
                'ForegroundColor', Core_UI.BLACK, ...
                'HorizontalAlignment', 'left', ...
                'FontSize', Core_UI.getFontSize(9), ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            editable_handle = uicontrol('Parent', box_handle,...
                'Style', 'edit',...
                'FontName', 'Courier New', ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'left', ...
                'UserData', property_name, ...
                'Callback', callback);
            uicontrol('Parent', box_handle,...
                'Style', 'pushbutton', ...
                'FontSize', Core_UI.getFontSize(7), ...
                'String', '...',...
                'Callback', @Core_UI.onSearchFileBox)
            box_handle.Widths = widths;
        end
        
        function [box_handle, editable_handle, flag_handle] = insertDirBox(parent, description, property_name, callback, widths)
            if nargin < 5
                widths  = [25 -1 -1 25];
            end
            box_handle = uix.Grid('Parent', parent, ...
                'Padding', 0, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            if numel(widths) == 4
                flag_handle = Core_UI.insertFlag(box_handle, {property_name});
            end
            uicontrol('Parent', box_handle, ...
                'Style', 'Text', ...
                'String', description, ...
                'ForegroundColor', Core_UI.BLACK, ...
                'HorizontalAlignment', 'left', ...
                'FontSize', Core_UI.getFontSize(9), ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            editable_handle = uicontrol('Parent', box_handle,...
                'Style', 'edit',...
                'FontName', 'Courier New', ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'left', ...
                'UserData', property_name, ...
                'Callback', callback);
            uicontrol('Parent', box_handle,...
                'Style', 'pushbutton', ...
                'FontSize', Core_UI.getFontSize(7), ...
                'String', '...',...
                'Callback', @Core_UI.onSearchDirBox)
            box_handle.Widths = widths;
            box_handle.Heights = 23;
        end
        
        function [box_handle, editable_handle, flag_handle] = insertDirBoxDark(parent, description, property_name, callback, widths)
            if nargin < 5
                widths  = [25 -1 -1 25];
            end
            box_handle = uix.Grid('Parent', parent, ...
                'Padding', 0, ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
            if numel(widths) == 4
                flag_handle = Core_UI.insertFlag(box_handle, {property_name});
            end
            uicontrol('Parent', box_handle, ...
                'Style', 'Text', ...
                'String', description, ...
                'ForegroundColor', Core_UI.WHITE, ...
                'HorizontalAlignment', 'left', ...
                'FontSize', Core_UI.getFontSize(9), ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
            editable_handle = uicontrol('Parent', box_handle,...
                'Style', 'edit',...
                'FontName', 'Courier New', ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'left', ...
                'UserData', property_name, ...
                'Callback', callback);
            uicontrol('Parent', box_handle,...
                'Style', 'pushbutton', ...
                'FontSize', Core_UI.getFontSize(7), ...
                'String', '...',...
                'Callback', @Core_UI.onSearchDirBox)
            box_handle.Widths = widths;
            box_handle.Heights = 23;
        end
        
        function [box_handle, editable_handle_dir, editable_handle_file, flag_handle] = insertDirFileBox(parent, description, property_name_dir, property_name_file, callback, widths)
            if nargin < 6
                widths  = [25 -1 -3 5 -1 25];
            end
            box_handle = uix.Grid('Parent', parent, ...
                'Padding', 0, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            if numel(widths) == 6
                flag_handle = Core_UI.insertFlag(box_handle, {property_name_dir, property_name_file});
            end
            uicontrol('Parent', box_handle, ...
                'Style', 'Text', ...
                'String', description, ...
                'ForegroundColor', Core_UI.BLACK, ...
                'HorizontalAlignment', 'left', ...
                'FontSize', Core_UI.getFontSize(9), ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            editable_handle_dir = uicontrol('Parent', box_handle,...
                'Style', 'edit',...
                'FontName', 'Courier New', ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'left', ...
                'UserData', property_name_dir, ...
                'Callback', callback);
            Core_UI.insertEmpty(box_handle);
            editable_handle_file = uicontrol('Parent', box_handle,...
                'Style', 'edit',...
                'FontName', 'Courier New', ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'left', ...
                'UserData', property_name_file, ...
                'Callback', callback);
            uicontrol('Parent', box_handle,...
                'Style', 'pushbutton', ...
                'FontSize', Core_UI.getFontSize(7), ...
                'String', '...',...
                'Callback', @Core_UI.onSearchDirFileBox)
                                 
            box_handle.Widths = widths;
        end
                
        function [box_handle, editable_handle_dir, editable_handle_file] = insertDirFileBoxObsML(parent, description, property_name_dir, property_name_file, callback, widths)
            if nargin < 6
                widths  = {[-1 -1 25], [-1 -1 25]};
            end
            box_handle = uix.VBox('Parent', parent, ...
                'Padding', 0, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            box_top = uix.HBox('Parent', box_handle, ...
                'Padding', 0, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            box_bot = uix.HBox('Parent', box_handle, ...
                'Padding', 0, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            
            uicontrol('Parent', box_top, ...
                'Style', 'Text', ...
                'String', [description 's directory'], ...
                'ForegroundColor', Core_UI.BLACK, ...
                'HorizontalAlignment', 'left', ...
                'FontSize', Core_UI.getFontSize(9), ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            editable_handle_dir = uicontrol('Parent', box_top,...
                'Style','edit',...
                'FontName', 'Courier New', ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'left', ...
                'UserData', property_name_dir, ...
                'Callback', callback);
            uicontrol('Parent', box_top,...
                'Style', 'pushbutton', ...
                'FontSize', Core_UI.getFontSize(7), ...
                'String', '...',...
                'Callback', @Core_UI.onSearchDirBoxMLStationObs)
            
            uicontrol('Parent', box_bot, ...
                'Style', 'Text', ...
                'String', [description ' files'], ...
                'ForegroundColor', Core_UI.BLACK, ...
                'HorizontalAlignment', 'left', ...
                'FontSize', Core_UI.getFontSize(9), ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            editable_handle_file = uicontrol('Parent', box_bot,...
                'Style','edit',...
                'FontName', 'Courier New', ...
                'FontWeight', 'bold', ...
                'Min',0,'Max',2,...  % This is the key to multiline edits.
                'HorizontalAlignment', 'left', ...
                'UserData', property_name_file, ...
                'Callback', callback);
            uicontrol('Parent', box_bot,...
                'Style', 'pushbutton', ...
                'FontSize', Core_UI.getFontSize(7), ...
                'String', '...',...
                'Callback', @Core_UI.onSearchDirFileBoxML)
            
            box_top.Widths = widths{1};
            box_bot.Widths = widths{1};
            box_handle.Heights = [23 -1];
        end
        
        function [box_handle, editable_handle_dir, editable_handle_file] = insertDirFileBoxMetML(parent, description, property_name_dir, property_name_file, callback, widths)
            if nargin < 6
                widths  = {[-1 -1 25], [-1 -1 25]};
            end
            box_handle = uix.VBox('Parent', parent, ...
                'Padding', 0, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            box_top = uix.HBox('Parent', box_handle, ...
                'Padding', 0, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            box_bot = uix.HBox('Parent', box_handle, ...
                'Padding', 0, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            
            uicontrol('Parent', box_top, ...
                'Style', 'Text', ...
                'String', [description ' directory'], ...
                'ForegroundColor', Core_UI.BLACK, ...
                'HorizontalAlignment', 'left', ...
                'FontSize', Core_UI.getFontSize(9), ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            editable_handle_dir = uicontrol('Parent', box_top,...
                'Style','edit', ...
                'FontName', 'Courier New', ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'left', ...
                'UserData', property_name_dir, ...
                'Callback', callback);
            uicontrol('Parent', box_top,...
                'Style', 'pushbutton', ...
                'FontSize', Core_UI.getFontSize(7), ...
                'String', '...',...
                'Callback', @Core_UI.onSearchDirBoxMLStationMet)
            
            uicontrol('Parent', box_bot, ...
                'Style', 'Text', ...
                'String', [description ' files'], ...
                'ForegroundColor', Core_UI.BLACK, ...
                'HorizontalAlignment', 'left', ...
                'FontSize', Core_UI.getFontSize(9), ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            editable_handle_file = uicontrol('Parent', box_bot,...
                'Style','edit',...
                'FontName', 'Courier New', ...
                'FontWeight', 'bold', ...
                'Min',0,'Max',2,...  % This is the key to multiline edits.
                'HorizontalAlignment', 'left', ...
                'UserData', property_name_file, ...
                'Callback', callback);
            uicontrol('Parent', box_bot,...
                'Style', 'pushbutton', ...
                'FontSize', Core_UI.getFontSize(7), ...
                'String', '...',...
                'Callback', @Core_UI.onSearchDirFileBoxML)
            
            box_top.Widths = widths{1};
            box_bot.Widths = widths{2};
            box_handle.Heights = [23 -1];
        end
        
        function [box_handle, editable_handle] = insertFileBoxML(parent, description, property_name, callback, widths)
            if nargin < 5
                widths  = [-1 -1 25];
            end
            box_handle = uix.Grid('Parent', parent, ...
                'Padding', 0, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            uicontrol('Parent', box_handle, ...
                'Style', 'Text', ...
                'String', description, ...
                'ForegroundColor', Core_UI.BLACK, ...
                'HorizontalAlignment', 'left', ...
                'FontSize', Core_UI.getFontSize(9), ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            editable_handle = uicontrol('Parent', box_handle,...
                'Style','edit',...
                'FontName', 'Courier New', ...
                'FontWeight', 'bold', ...
                'Min',0,'Max',2,...  % This is the key to multiline edits.
                'HorizontalAlignment', 'left', ...
                'UserData', property_name, ...
                'Callback', callback);
            uicontrol('Parent', box_handle,...
                'Style', 'pushbutton', ...
                'FontSize', Core_UI.getFontSize(7), ...
                'String', '...',...
                'Callback', @Core_UI.onSearchFileBox)
            box_handle.Widths = widths;
            box_handle.Heights = -1;
        end
        
        function [box_handle, editable_handle] = insertEditBox(parent, text_left, property_name, text_right, callback, widths)
            if nargin < 6
                if nargin < 4 || isempty(text_right)
                    widths  = [-1 -1];
                else
                    widths  = [-1 -1 5 -1];
                end
            end
            
            box_handle = uix.Grid('Parent', parent, ...
                'Padding', 0, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            uicontrol('Parent', box_handle, ...
                'Style', 'Text', ...
                'String', text_left, ...
                'ForegroundColor', Core_UI.BLACK, ...
                'HorizontalAlignment', 'left', ...
                'FontSize', Core_UI.getFontSize(9), ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            editable_handle = uicontrol('Parent', box_handle,...
                'Style', 'edit',...
                'UserData', property_name, ...
                'Callback', callback);
            Core_UI.insertEmpty(box_handle, Core_UI.LIGHT_GREY_BG);
            uicontrol('Parent', box_handle, ...
                'Style', 'Text', ...
                'String', text_right, ...
                'ForegroundColor', Core_UI.BLACK, ...
                'HorizontalAlignment', 'left', ...
                'FontSize', Core_UI.getFontSize(9), ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            box_handle.Widths = widths;
            box_handle.Heights = 23;
        end
        
        function [box_handle, editable_handle] = insertEditBoxDark(parent, text_left, property_name, text_right, callback, widths)
            if nargin < 6
                if nargin < 4 || isempty(text_right)
                    widths  = [-1 -1];
                else
                    widths  = [-1 -1 5 -1];
                end
            end
            
            box_handle = uix.Grid('Parent', parent, ...
                'Padding', 0, ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
            uicontrol('Parent', box_handle, ...
                'Style', 'Text', ...
                'String', text_left, ...
                'ForegroundColor', Core_UI.WHITE, ...
                'HorizontalAlignment', 'left', ...
                'FontSize', Core_UI.getFontSize(9), ...
                'BackgroundColor', Core_UI.DARK_GREY_BG);
            editable_handle = uicontrol('Parent', box_handle,...
                'Style', 'edit',...
                'UserData', property_name, ...
                'Callback', callback);
            if nargin > 3 && ~isempty(text_right)
                Core_UI.insertEmpty(box_handle, Core_UI.DARK_GREY_BG);
                uicontrol('Parent', box_handle, ...
                    'Style', 'Text', ...
                    'String', text_right, ...
                    'ForegroundColor', Core_UI.WHITE, ...
                    'HorizontalAlignment', 'left', ...
                    'FontSize', Core_UI.getFontSize(9), ...
                    'BackgroundColor', Core_UI.DARK_GREY_BG);
            end
            box_handle.Widths = widths;
            box_handle.Heights = 23;
        end
        
        function [box_handle] = insertEditBoxArray(parent,num_boxes, text_left, property_name, text_right, callback, widths)
            if nargin < 6
                if nargin < 4 || isempty(text_right)
                    widths  = [-1 -1];
                else
                    widths  = [-1 -1 5 -1];
                end
            end
            
            box_handle = uix.Grid('Parent', parent, ...
                'Padding', 0, ...
                'UserData', property_name, ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            uicontrol('Parent', box_handle, ...
                'Style', 'Text', ...
                'String', text_left, ...
                'ForegroundColor', Core_UI.BLACK, ...
                'HorizontalAlignment', 'left', ...
                'FontSize', Core_UI.getFontSize(9), ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            for i = 1 : num_boxes
                uicontrol('Parent', box_handle,...
                    'Style', 'edit',...
                    'UserData', property_name, ...
                    'Callback', callback);
                  Core_UI.insertEmpty(box_handle, Core_UI.LIGHT_GREY_BG);
            end
          
            uicontrol('Parent', box_handle, ...
                'Style', 'Text', ...
                'String', text_right, ...
                'ForegroundColor', Core_UI.BLACK, ...
                'HorizontalAlignment', 'left', ...
                'FontSize', Core_UI.getFontSize(9), ...
                'BackgroundColor', Core_UI.LIGHT_GREY_BG);
            box_handle.Widths = [widths(1) repmat(widths(2:3), 1, num_boxes) widths(4)];
            box_handle.Heights = 23;
        end
        
        function [j_edit_box, h_log_panel] = insertLog(parent)
            % Insert a log editbox
            %
            % SYNTAX
            %   [j_edit_box, h_log_panel] = Core_UI.insertLog(parent)
            h_log_panel = uicontrol('style','edit', 'max', 5, 'Parent', parent, 'Units','norm', 'Position', [0,0.2,1,0.8], 'Background', 'w');
            % Get the underlying Java editbox, which is contained within a scroll-panel
            jScrollPanel = findjobj(h_log_panel);
            try
                jScrollPanel.setVerticalScrollBarPolicy(jScrollPanel.java.VERTICAL_SCROLLBAR_AS_NEEDED);
                jScrollPanel = jScrollPanel.getViewport;
            catch
                % may possibly already be the viewport, depending on release/platform etc.
            end
            j_edit_box = handle(jScrollPanel.getView,'CallbackProperties');
        end
    end     
    %% METHODS ELEMENT MODIFIER
    % ==================================================================================================================================================
    methods (Static)
        function checkFlag(flag_handle_list)
            % Check with the content of the properties of state if a flag should be of a certain color
            %
            % SYNTAX
            %   Core_UI.checkFlag(flag_handle_list)
            
            fnp = File_Name_Processor;
            
            if ~iscell(flag_handle_list)
                flag_handle_list = {flag_handle_list};            
            end
            for f = 1 : numel(flag_handle_list)
                flag_handle = flag_handle_list{f};
                properties = flag_handle.UserData;
                if isempty(properties)
                    Core_UI.setFlagGrey(flag_handle);
                else
                    state = Core.getState();
                    try
                        if numel(properties) == 2
                            full_file = fullfile(fnp.getFullDirPath(state.(properties{1})), state.(properties{2}));
                        else
                            full_file = File_Name_Processor.getFullDirPath(state.(properties{1}));
                        end                        
                    catch
                        % Properties are not a field of state
                        full_file = '';
                    end
                    if exist(full_file, 'file')
                        Core_UI.setFlagGreen(flag_handle);
                    else
                        if any(full_file == '$')
                            % limited check support, check only beginning and end (for speed- up)
                            full_file_sss0 = fnp.dateKeyRep(full_file, state.sss_date_start, '0');
                            full_file_sss1 = fnp.dateKeyRep(full_file, state.sss_date_stop, '0');
                            if exist(full_file_sss0, 'file') && exist(full_file_sss1, 'file')                        
                                Core_UI.setFlagGreen(flag_handle);
                            else
                                Core_UI.setFlagOrange(flag_handle);
                            end
                        else
                            if isempty(full_file)
                                Core_UI.setFlagGrey(flag_handle);
                            else
                                Core_UI.setFlagRed(flag_handle);
                            end
                        end
                    end
                end
            end                        
        end
        
        function setFlagRed(flag_handle)
            % Change flag color to red
            %
            % SYNTAX
            %   COre_UI.setFlagRed(flag_handle)
            save_properties = flag_handle.UserData;
            hold(flag_handle, 'off')
            plot(flag_handle, 0, 0, '.', 'MarkerSize', 45, 'Color', max(0, Core_UI.RED - 0.5)); 
            hold(flag_handle, 'on')
            plot(flag_handle, 0, 0, '.', 'MarkerSize', 33, 'Color', min(1, Core_UI.RED + 0.3));
            flag_handle.XTickLabel = [];
            flag_handle.YTickLabel = [];                    
            axis(flag_handle, 'off');
            flag_handle.UserData = save_properties;            
        end
        
        function setFlagOrange(flag_handle)
            % Change flag color to orange
            %
            % SYNTAX
            %   COre_UI.setFlagOrange(flag_handle)
            save_properties = flag_handle.UserData;
            hold(flag_handle, 'off')
            plot(flag_handle, 0, 0, '.', 'MarkerSize', 45, 'Color', max(0, Core_UI.ORANGE - 0.5)); 
            hold(flag_handle, 'on')
            plot(flag_handle, 0, 0, '.', 'MarkerSize', 33, 'Color', min(1, Core_UI.ORANGE + 0.1));
            flag_handle.XTickLabel = [];
            flag_handle.YTickLabel = [];                    
            axis(flag_handle, 'off');
            flag_handle.UserData = save_properties;            
        end
        
        function setFlagGreen(flag_handle)
            % Change flag color to green
            %
            % SYNTAX
            %   COre_UI.setFlagGreen(flag_handle)
            save_properties = flag_handle.UserData;
            hold(flag_handle, 'off')
            plot(flag_handle, 0, 0, '.', 'MarkerSize', 45, 'Color', max(0, Core_UI.GREEN - 0.5)); 
            hold(flag_handle, 'on')
            plot(flag_handle, 0, 0, '.', 'MarkerSize', 33, 'Color', min(1, Core_UI.GREEN + 0.3));
            flag_handle.XTickLabel = [];
            flag_handle.YTickLabel = [];                    
            axis(flag_handle, 'off');
            flag_handle.UserData = save_properties;            
        end
        
        function setFlagGrey(flag_handle)
            % Change flag color to grey
            %
            % SYNTAX
            %   COre_UI.setFlagGrey(flag_handle)
            save_properties = flag_handle.UserData;
            hold(flag_handle, 'off')
            plot(flag_handle, 0, 0, '.', 'MarkerSize', 45, 'Color', max(0, Core_UI.GREY - 0.5)); 
            hold(flag_handle, 'on')
            plot(flag_handle, 0, 0, '.', 'MarkerSize', 33, 'Color', Core_UI.GREY);
            flag_handle.XTickLabel = [];
            flag_handle.YTickLabel = [];                    
            axis(flag_handle, 'off');
            flag_handle.UserData = save_properties;            
        end

        function disableElement(el)
            Core_UI.enableElement(el, false)
        end
        
        function enableElement(el, status)
            if nargin == 1
                status = true;
            end
            for e = 1 : numel(el)
                if isprop(el, 'Enable')
                    %switch el(e).Style
                    %    case {'edit', 'text'}
                            if status
                                el(e).Enable = 'on';
                            else
                                el(e).Enable = 'off';
                            end
                    %end
                else % is Grid, VBox, HBox, ...
                    obj = findobj(el);
                    if numel(obj) > 1
                        Core_UI.enableElement(obj(2:end), status);
                    end
                end
            end
        end
        
        function guiAddMessage(j_edit_box, text, severity)
            % Log messages in a Log element
            %
            % SYNTAX
            %   Core_UI.guiAddMessage(jEditbox, text, severity)
            
            % Ensure we have an HTML-ready editbox
            HTMLclassname = 'javax.swing.text.html.HTMLEditorKit';
            if ~isa(j_edit_box.getEditorKit,HTMLclassname)
                j_edit_box.setContentType('text/html');
            end
            % Parse the severity and prepare the HTML message segment
            if nargin < 3,  severity='none';  end
            switch lower(severity(1))
                %case 'm',  icon = 'greenarrowicon.gif'; color='black';
                %case 'w',  icon = 'warning.gif';        color='orange'; text = ['WARNING: ', text];                
                %case 'e',  icon = 'demoicon.gif';       color='red';    text = ['ERROR: ', text];
                case 'm',  icon = 'info32b.png';          color='#222222';
                case 'w',  icon = 'alert-yellow32b.png';  color='orange'; text = ['WARNING: ', text];                
                case 'e',  icon = 'error32b.png';         color='red';    text = ['ERROR: ', text];
                case 'v',  icon = 'ok32b.png';            color='#222222';
                case 'x',  icon = 'ko32b.png';            color='#222222';
                case '-',  icon = 'off32b.png';           color='#222222';
                otherwise, icon = 'empty32b.png';         color='gray';
            end
                   
            msg_txt = ['<font color=', color, ' face="Helvetica, Arial, sans-serif">', text, '</font>'];
            
            %icon = fullfile(matlabroot,'toolbox/matlab/icons',icon);
            icon = fullfile(Core.getInstallDir, '../icons', icon);
            icon_txt =['<img src="file:///', icon, '" height=20 width=16 style="height: 20px; width: 16px;"/>'];
            new_text = ['<table style="width:100%;"><tr><td style="width:16px;vertical-align:top">' icon_txt '</td><td>' msg_txt '</td></tr></table>'];
            
            % Place the HTML message segment at the bottom of the editbox
            doc = j_edit_box.getDocument();
            j_edit_box.getEditorKit().read(java.io.StringReader(new_text), doc, doc.getLength());
            j_edit_box.setCaretPosition(doc.getLength());

            end_pos = j_edit_box.getDocument.getLength;
            j_edit_box.setCaretPosition(end_pos); % end of content
        end
        
        function guiAddHTML(j_edit_box, html_txt)
            % Log messages in a Log element
            %
            % SYNTAX
            %   Core_UI.guiAddMessage(jEditbox, text, severity)
            
            % Ensure we have an HTML-ready editbox
            HTMLclassname = 'javax.swing.text.html.HTMLEditorKit';
            if ~isa(j_edit_box.getEditorKit,HTMLclassname)
                j_edit_box.setContentType('text/html');
            end
                        
            % Place the HTML message segment at the bottom of the editbox
            doc = j_edit_box.getDocument();
            j_edit_box.getEditorKit().read(java.io.StringReader(html_txt), doc, doc.getLength());
            j_edit_box.setCaretPosition(doc.getLength());

            end_pos = j_edit_box.getDocument.getLength;
            j_edit_box.setCaretPosition(end_pos); % end of content
        end
        
        function guiClearLog(j_edit_box, html_txt)
            % Empty the logging Window
            %
            % SYNTAX
            %   Core_UI.guiClearLog(jEditbox, text, severity)
            
            % Ensure we have an HTML-ready editbox
            HTMLclassname = 'javax.swing.text.html.HTMLEditorKit';
            if ~isa(j_edit_box.getEditorKit,HTMLclassname)
                j_edit_box.setContentType('text/html');
            end
                        
            % Place the HTML message segment at the bottom of the editbox
            j_edit_box.setText('');% end of content
        end
    end
    %% METHODS EVENTS
    % ==================================================================================================================================================
    methods (Static, Access = public)
        function onSearchDirFileBox(caller, event)
            % [file, path] = uigetfile({'*.*',  'All Files (*.*)'}, this.state.getHomeDir);
            [file, path] = uigetfile([caller.Parent.Children(4).String filesep '*.*']);
            if file ~= 0
                caller.Parent.Children(2).String = file;
                caller.Parent.Children(4).String = path;
                % trigger the edit of the field
                callbackCell = get(caller.Parent.Children(2),'Callback');
                callbackCell(caller.Parent.Children(2));
                callbackCell = get(caller.Parent.Children(4),'Callback');
                callbackCell(caller.Parent.Children(4));
            end
        end
        
        function onSearchDirFileBoxML(caller, event)
            % [file, path] = uigetfile({'*.*',  'All Files (*.*)'}, this.state.getHomeDir);
            [file, path] = uigetfile([caller.Parent.Parent.Children(2).Children(2).String filesep '*.*'], 'MultiSelect', 'on');
            if ~iscell(file)
                if file ~= 0
                    file = {file};
                end
            end
            if ~isempty(file) && iscell(file)
                caller.Parent.Children(2).String = file;
                caller.Parent.Parent.Children(2).Children(2).String = path;
                % trigger the edit of the field
                callbackCell = get(caller.Parent.Children(2),'Callback');
                callbackCell(caller.Parent.Children(2));
                callbackCell = get(caller.Parent.Parent.Children(2).Children(2),'Callback');
                callbackCell(caller.Parent.Parent.Children(2).Children(2));
            end
        end
        
        function onSearchFileBox(caller, event)
            state = Core.getCurrentSettings();
            [file, path] = uigetfile([state.getHomeDir filesep '*.*']);
            if file ~= 0
                caller.Parent.Children(2).String = file;
                % trigger the edit of the field
                callbackCell = get(caller.Parent.Children(2),'Callback');
                callbackCell(caller.Parent.Children(2));
            end
        end
        
        function onSearchDirBox(caller, event)
            dir_path = uigetdir(caller.Parent.Children(2).String);
            if dir_path ~= 0
                caller.Parent.Children(2).String = dir_path;
                % trigger the edit of the field
                callbackCell = get(caller.Parent.Children(2),'Callback');
                callbackCell(caller.Parent.Children(2));
            end
        end
        
        function onSearchDirBoxMLStationObs(caller, event)
            dir_path = uigetdir(caller.Parent.Children(2).String);
            station_list = Core_Utils.getStationList(dir_path,'oO');
            
            if dir_path ~= 0
                caller.Parent.Children(2).String = dir_path;
                % trigger the edit of the field
                callbackCell = get(caller.Parent.Children(2),'Callback');
                callbackCell(caller.Parent.Children(2));
                if ~isempty(station_list)
                    caller.Parent.Parent.Children(1).Children(2).String = station_list;
                    % trigger the edit of the field
                    callbackCell = get(caller.Parent.Parent.Children(1).Children(2),'Callback');
                    callbackCell(caller.Parent.Parent.Children(1).Children(2));
                end
            end
        end
        
        function onGetRecursiveMarkers(caller, event)
            dir_box = caller.Parent.Parent.Children(end).Children(end).Children(2);
            list_box = caller.Parent.Parent.Children(end).Children(1).Children(2);
            
            dir_path = uigetdir(dir_box.String);
            station_list = Core_Utils.getStationList(dir_path, 'oO', true);
            
            if dir_path ~= 0
                dir_box.String = dir_path;
                % trigger the edit of the field
                callbackCell = get(dir_box,'Callback');
                callbackCell(dir_box);
                if ~isempty(station_list)
                    list_box.String = station_list;
                    % trigger the edit of the field
                    callbackCell = get(list_box,'Callback');
                    callbackCell(list_box);
                end
            end
        end
        
        function onSearchDirBoxMLStationMet(caller, event)
            dir_path = uigetdir(caller.Parent.Children(2).String);
            station_list = Core_Utils.getStationList(dir_path,'mM');
            
            if dir_path ~= 0
                caller.Parent.Children(2).String = dir_path;
                % trigger the edit of the field
                callbackCell = get(caller.Parent.Children(2),'Callback');
                callbackCell(caller.Parent.Children(2));
                if ~isempty(station_list)
                    caller.Parent.Parent.Children(1).Children(2).String = station_list;
                    % trigger the edit of the field
                    callbackCell = get(caller.Parent.Parent.Children(1).Children(2),'Callback');
                    callbackCell(caller.Parent.Parent.Children(1).Children(2));
                end
            end
        end        
    end
    %% METHODS getters
    % ==================================================================================================================================================
    methods
        function ok_go = isGo(this)
            ok_go = this.main.isGo();
        end
    end    
    %% METHODS (static) getters
    % ==================================================================================================================================================
    methods (Static)
        function os_size = getSizeConversion()
            if isunix()
                if ismac()
                    os_size = Core_UI.FONT_SIZE_CONVERSION_MAC;
                else
                    os_size = Core_UI.FONT_SIZE_CONVERSION_LNX;
                end
            else
                os_size = Core_UI.FONT_SIZE_CONVERSION_WIN;
            end
        end
        
        function font_o = getFontSize(font_i)
            font_o = round(font_i * Core_UI.getSizeConversion());
        end
        
        function color = getColor(color_id, color_num)
            % get a color taken from COLOR ORDER CONSTANT
            %
            % SYNTAX:
            %   Core_UI.getColor(id);
            %
            if (nargin == 2) && (color_num > 7)
                color = Cmap.getColor(color_id, color_num, 'linspaced');
            else
                if nargin == 1
                    color = Cmap.getColor(color_id, 7);
                else
                    color = Cmap.getColor(color_id, color_num);
                end
            end
        end
        
        function logo_ax = insertLogo(container, location)  
            % Insert a new axis containing goGPS Logo
            %
            % SYNTAX 
            %   logo_ax = Mesonet.insertLogo(container)  
            
            if nargin < 1 || isempty(container)
                container = gcf;
            end
            
            if (nargin < 2) || isempty(location)
                location = 'SouthEast';
            end
            logo_ax = axes( 'Parent', container);
            % [logo, ~, transparency] = imread('./reserved/GReD logo_86px.png');
            [logo, transparency] = Core_UI.getLogo();
            logo(repmat(sum(logo,3) == 0,1,1,3)) = 0;
            logo = logo - 20;
            image(logo_ax, logo, 'AlphaData', transparency);
            logo_ax.XTickLabel = [];
            logo_ax.YTickLabel = [];
            
            switch location
                case 'NorthEast'
                    logo_ax.Units = 'pixels';
                    logo_ax.Position([3 4]) = [size(logo,2) size(logo,1)];            
                    cu = container.Units;
                    container.Units = 'pixels';
                    logo_ax.Position(1) = container.Position(3) - logo_ax.Position(3) - 10;
                    logo_ax.Position(2) = container.Position(4) - logo_ax.Position(4) - 10;
                    logo_ax.Units = 'normalized';            
                    container.Units = cu;
                case 'NorthWest'
                    logo_ax.Units = 'pixels';
                    logo_ax.Position([3 4]) = [size(logo,2) size(logo,1)];            
                    cu = container.Units;
                    container.Units = 'pixels';
                    logo_ax.Position(1) = 10;
                    logo_ax.Position(2) = container.Position(4) - logo_ax.Position(4) - 10;
                    logo_ax.Units = 'normalized';            
                    container.Units = cu;
                case 'SouthEast'
                    logo_ax.Units = 'pixels';
                    logo_ax.Position([3 4]) = [size(logo,2) size(logo,1)];            
                    logo_ax.Position(2) = 10;
                    cu = container.Units;
                    container.Units = 'pixels';
                    logo_ax.Position(1) = container.Position(3) - logo_ax.Position(3) - 10;
                    logo_ax.Units = 'normalized';            
                    container.Units = cu;
                case 'SouthWest'
                    logo_ax.Units = 'pixels';
                    logo_ax.Position([3 4]) = [size(logo,2) size(logo,1)];            
                    logo_ax.Position(2) = 10;
                    cu = container.Units;
                    container.Units = 'pixels';
                    logo_ax.Position(1) = 10;
                    logo_ax.Units = 'normalized';            
                    container.Units = cu;
            end            
            axis off;
            if Core.isGReD
                GReD_Utility.insertLogoGReD(container, location);
            end

        end
    end
end
