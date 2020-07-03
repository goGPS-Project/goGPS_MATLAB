%   CLASS GUI_Unique_Win
% =========================================================================
%
% DESCRIPTION
%   Parent class for goGPS windows
%
%
% FOR A LIST OF CONSTANTs and METHODS use doc GUI_Unique_Win


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
%  Contributors:     Andrea Gatti
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

classdef GUI_Unique_Win < handle   
    properties (Constant, Abstract)
        WIN_NAME
    end
    
    %% PROPERTIES SINGLETON POINTERS
    % ==================================================================================================================================================
    
    %% PROPERTIES GUI
    % ==================================================================================================================================================
    properties (Abstract)
        win       % Handle of the this window       
    end    
   
    %% METHOD CREATOR
    % ==================================================================================================================================================
    methods (Static, Access = protected)
        function this = GUI_Unique_Win()
            % GUI_Unique_Win object creator
        end
    end
    %% METHODS UTILITY
    % ==================================================================================================================================================
    methods                                
        function fig_handle = getUniqueWinHandle(this)
            
            fig_handle = [];
            if ~isempty(this.win) && isvalid(this.win)
                % if the win is open and stored in this singleton object
                fig_handle = this.win;
            end
            % clean way of doing this:
            % fig_handle = findobj(get(groot, 'Children'), 'UserData', this.WIN_NAME);
            
            % fast way of doing this:
            fh_list = get(groot, 'Children');
            fh_list = setdiff(fh_list, fig_handle);
            
            % bad code writing style but fast
            for f = 1 : numel(fh_list)
                try
                    if isfield(fh_list(f).UserData, 'name') && strcmp(fh_list(f).UserData.name, this.WIN_NAME)
                        % If there are lone Edit figures close them
                        delete(fh_list(f));
                    end
                catch ex
                    Core_Utils.printEx(ex);
                end
            end
        end         
        
        function bringOnTop(this)
            % Bring the window on top of all the others
            this.win.Visible = 'on';
            warning off
            % Ignore :-( JavaFrame property will be obsoleted in a future release
            j_win = get(handle(this.win),'JavaFrame');
            warning on
            if isHG2
                j_frame = j_win.fHG2Client;
            else
                j_frame = j_win.fHG1Client;
            end
            drawnow
            try
                j_frame.getWindow.setAlwaysOnTop(1);
                j_frame.getWindow.setAlwaysOnTop(0);
            catch
            end
        end
    end       
end
