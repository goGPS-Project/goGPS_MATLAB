%   CLASS Mode_Settings
% =========================================================================
%
% DESCRIPTION
%   Class to store all the goGPS mode
%
% EXAMPLE
%   settings = Mode_Settings();
%
% FOR A LIST OF CONSTANTS and METHODS use doc Mode_Settings

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___ 
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.9.1
% 
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Gatti Andrea
%  Contributors:     Gatti Andrea, ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011 
%--------------------------------------------------------------------------

classdef Mode_Settings < Settings_Interface
    
    properties (Constant, GetAccess = public)
        % Deprecated modes, in a future version of goGPS modes will be managed in a different way     
        MODE_RT_NAV          = 24;  % Real Time Navigation (Kalman Filter on Code and Phase Double Differences (with/without a constraint)
        MODE_RT_R_MON        = 21;  % Real Time Rover Monitor
        MODE_RT_M_MON        = 22;  % Real Time Master Monitor
        MODE_RT_RM_MON       = 23;  % Real Time Master + Rover Monitor
        
        MODE_PP_LS_C_SA      = 1;   % Post Proc Least Squares on Code Stand Alone
        MODE_PP_LS_CP_SA     = 3;   % Post Proc Least Squares on Code and Phase Stand Alone
        MODE_PP_LS_CP_VEL    = 5; % Post Proc Least Squares on Code and Phase for Velocity estimation
        MODE_PP_LS_C_DD      = 11;  % Post Proc Least Squares on Code Double Differences
        MODE_PP_LS_CP_DD_L   = 13;  % Post Proc Least Squares on Code and Phase Double Differences with LAMBDA
        MODE_PP_LS_CP_DD_MR  = 16;  % Post Proc Least Squares on Code and Phase Double Differences, Multiple Receivers
        MODE_PP_LS_C_SA_MR   = 17;  % Post Proc Least Squares on Code Stand Alone, Multiple Receivers
        
        MODE_PP_KF_C_SA          = 2;   % Post Proc Kalman Filter on Code Stand Alone
        MODE_PP_KF_C_DD          = 12;  % Post Proc Kalman Filter on Code Double Differences
        MODE_PP_KF_CP_SA         = 4;   % Post Proc Kalman Filter on Code and Phase Stand Alone (PPP)
        MODE_PP_KF_CP_DD         = 14;  % Post Proc Kalman Filter on Code and Phase Double Differences
        MODE_PP_KF_CP_DD_MR      = 15;  % Post Proc Kalman Filter on Code and Phase Double Differences, Multiple Receivers
        MODE_PP_SEID_PPP         = 18;  % SEID followed by PPP (Kalman Filter on Code and Phase Stand Alone (PPP)) it is both stand alone and DD
        
        % goGPS MODES -----------------------------------------------------
                 
        % Group of post processing modes
        GMODE_PP = [ Mode_Settings.MODE_PP_LS_C_SA ...     
            Mode_Settings.MODE_PP_LS_CP_SA ...
            Mode_Settings.MODE_PP_LS_C_DD ...
            Mode_Settings.MODE_PP_LS_CP_DD_L ...
            Mode_Settings.MODE_PP_LS_CP_VEL ...
            Mode_Settings.MODE_PP_KF_C_SA ...
            Mode_Settings.MODE_PP_KF_C_DD ...
            Mode_Settings.MODE_PP_KF_CP_SA ...
            Mode_Settings.MODE_PP_KF_CP_DD ...
            Mode_Settings.MODE_PP_LS_C_SA_MR ...
            Mode_Settings.MODE_PP_LS_CP_DD_MR ...
            Mode_Settings.MODE_PP_KF_CP_DD_MR ...
            Mode_Settings.MODE_PP_SEID_PPP];
        
        % Group of real time modes
        GMODE_RT = [ Mode_Settings.MODE_RT_NAV ...         
            Mode_Settings.MODE_RT_R_MON ...
            Mode_Settings.MODE_RT_M_MON ...
            Mode_Settings.MODE_RT_RM_MON];
        
        % Group of monitor modes
        GMODE_MON = [ Mode_Settings.MODE_RT_R_MON ...      
            Mode_Settings.MODE_RT_M_MON ...
            Mode_Settings.MODE_RT_RM_MON];

        % Group of stand alone modes
        GMODE_SA = [ Mode_Settings.MODE_PP_LS_C_SA ...     
            Mode_Settings.MODE_PP_LS_CP_SA ...
            Mode_Settings.MODE_PP_LS_CP_VEL ...
            Mode_Settings.MODE_PP_KF_C_SA ...
            Mode_Settings.MODE_PP_KF_CP_SA ...
            Mode_Settings.MODE_PP_LS_C_SA_MR];
        
        % Group of double differences modes
        GMODE_DD = [ Mode_Settings.MODE_PP_LS_C_DD ...     
            Mode_Settings.MODE_PP_LS_CP_DD_L ...
            Mode_Settings.MODE_PP_KF_C_DD ...
            Mode_Settings.MODE_PP_KF_CP_DD ...
            Mode_Settings.MODE_PP_LS_CP_DD_MR ...
            Mode_Settings.MODE_PP_KF_CP_DD_MR ...
            Mode_Settings.MODE_PP_SEID_PPP];
        
        % Group of multi-receiver modes
        GMODE_MR = [ Mode_Settings.MODE_PP_LS_C_SA_MR ...  
            Mode_Settings.MODE_PP_LS_CP_DD_MR ...
            Mode_Settings.MODE_PP_KF_CP_DD_MR ...
            Mode_Settings.MODE_PP_SEID_PPP];
       
        % Group of modes using Phase
        GMODE_PH = [ Mode_Settings.MODE_RT_NAV ...         
            Mode_Settings.MODE_RT_R_MON ...
            Mode_Settings.MODE_RT_M_MON ...
            Mode_Settings.MODE_RT_RM_MON ...
            Mode_Settings.MODE_PP_LS_CP_SA ...
            Mode_Settings.MODE_PP_LS_CP_VEL ...
            Mode_Settings.MODE_PP_LS_CP_DD_MR ...
            Mode_Settings.MODE_PP_LS_CP_DD_L ...
            Mode_Settings.MODE_PP_KF_CP_SA ...
            Mode_Settings.MODE_PP_KF_CP_DD ...
            Mode_Settings.MODE_PP_KF_CP_DD_MR ...
            Mode_Settings.MODE_PP_SEID_PPP];
        
        % Group of modes using Kalman Filter
        GMODE_KM = [ Mode_Settings.MODE_PP_KF_C_SA ...      
            Mode_Settings.MODE_PP_KF_C_DD ... 
            Mode_Settings.MODE_PP_KF_CP_SA ...
            Mode_Settings.MODE_PP_KF_CP_DD ...
            Mode_Settings.MODE_PP_KF_CP_DD_MR ...
            Mode_Settings.MODE_PP_SEID_PPP];
        
    end
    
    properties (Constant, GetAccess = protected)
        % This constant provides the connection among the numeric value of the mode with the string id and the set of 4 combobox used in the GUI to select a mode (id / GUI code / internal mode code)
        P_MODE_2_ID = [ 01, 1100, Mode_Settings.MODE_RT_NAV;
                        02, 1200, Mode_Settings.MODE_RT_R_MON;
                        03, 1300, Mode_Settings.MODE_RT_M_MON;
                        04, 1400, Mode_Settings.MODE_RT_RM_MON;
                        05, 2011, Mode_Settings.MODE_PP_LS_C_SA;
                        06, 0000, Mode_Settings.MODE_PP_LS_CP_SA;
                        07, 2014, Mode_Settings.MODE_PP_LS_CP_VEL;
                        08, 2012, Mode_Settings.MODE_PP_LS_C_DD;
                        09, 2013, Mode_Settings.MODE_PP_LS_CP_DD_L;
                        10, 0000, Mode_Settings.MODE_PP_LS_CP_DD_MR;
                        11, 0000, Mode_Settings.MODE_PP_LS_C_SA_MR;
                        12, 2021, Mode_Settings.MODE_PP_KF_C_SA;
                        13, 2022, Mode_Settings.MODE_PP_KF_C_DD;
                        14, 2023, Mode_Settings.MODE_PP_KF_CP_SA;
                        15, 2024, Mode_Settings.MODE_PP_KF_CP_DD;
                        16, 2031, Mode_Settings.MODE_PP_KF_CP_DD_MR;
                        17, 2032, Mode_Settings.MODE_PP_SEID_PPP;];                         
        
        P_MODE = { sprintf('%02d: %s', Mode_Settings.P_MODE_2_ID(1,3),'Real Time Navigation (Kalman Filter on Code and Phase Double Differences (with/without a constraint)'), ...
                       sprintf('%02d: %s', Mode_Settings.P_MODE_2_ID(2,3),'Real Time Rover Monitor'), ...
                       sprintf('%02d: %s', Mode_Settings.P_MODE_2_ID(3,3),'Real Time Master Monitor'), ...
                       sprintf('%02d: %s', Mode_Settings.P_MODE_2_ID(4,3),'Real Time Master + Rover Monitor'), ...
                       sprintf('%02d: %s', Mode_Settings.P_MODE_2_ID(5,3),'Post Proc Least Squares on Code Stand Alone'), ...
                       sprintf('%02d: %s', Mode_Settings.P_MODE_2_ID(6,3),'Post Proc Least Squares on Code and Phase Stand Alone'), ...
                       sprintf('%02d: %s', Mode_Settings.P_MODE_2_ID(7,3),'Post Proc Least Squares on Code and Phase for Velocity estimation'), ...
                       sprintf('%02d: %s', Mode_Settings.P_MODE_2_ID(8,3),'Post Proc Least Squares on Code Double Differences'), ...
                       sprintf('%02d: %s', Mode_Settings.P_MODE_2_ID(9,3),'Post Proc Least Squares on Code and Phase Double Differences with LAMBDA'), ...
                       sprintf('%02d: %s', Mode_Settings.P_MODE_2_ID(10,3),'Post Proc Least Squares on Code and Phase Double Differences, Multiple Receivers'), ...
                       sprintf('%02d: %s', Mode_Settings.P_MODE_2_ID(11,3),'Post Proc Least Squares on Code Stand Alone, Multiple Receivers'), ...
                       sprintf('%02d: %s', Mode_Settings.P_MODE_2_ID(12,3),'Post Proc Kalman Filter on Code Stand Alone'), ...
                       sprintf('%02d: %s', Mode_Settings.P_MODE_2_ID(13,3),'Post Proc Kalman Filter on Code Double Differences'), ...
                       sprintf('%02d: %s', Mode_Settings.P_MODE_2_ID(14,3),'Post Proc Kalman Filter on Code and Phase Stand Alone (PPP)'), ...
                       sprintf('%02d: %s', Mode_Settings.P_MODE_2_ID(15,3),'Post Proc Kalman Filter on Code and Phase Double Differences'), ...
                       sprintf('%02d: %s', Mode_Settings.P_MODE_2_ID(16,3),'Post Proc Kalman Filter on Code and Phase Double Differences, Multiple Receivers (SEID - only rinex writing)'), ...
                       sprintf('%02d: %s', Mode_Settings.P_MODE_2_ID(17,3),'SEID followed by PPP (Kalman Filter on Code and Phase Stand Alone (PPP)) it is both stand alone and DD') };
    end
    
    properties (SetAccess = protected, GetAccess = public)
        %------------------------------------------------------------------
        % PROCESSING PARAMETERS 
        %------------------------------------------------------------------
        
        % Processing mode
        p_mode = 1;        
    end
        
    % =========================================================================
    %  INIT
    % =========================================================================
    methods
        function this = Mode_Settings()
            % Creator of Mode_Settings
        end
    end
        
    % =========================================================================
    %  INTERFACE REQUIREMENTS
    % =========================================================================
    methods
        function import(this, settings)
            % This function import Mode (only) settings from another setting object
            if isa(settings, 'Ini_Manager')                
                this.p_mode = settings.getData('p_mode');
            else                
                this.p_mode = settings.p_mode;
            end
        end
        
        function str = toString(this, str)
            % Display the satellite system in use
            if (nargin == 1)
                str = '';
            end            
            str = [str sprintf(' Processing using %s\n\n', this.P_MODE{this.P_MODE_2_ID(this.P_MODE_2_ID(:,3) == this.p_mode, 1)})];

        end
        
        function str_cell = export(this, str_cell)            
            % Conversion to string ini format of the minimal information needed to reconstruct the this            
            if (nargin == 1)
                str_cell = {};
            end
            
            str_cell = Ini_Manager.toIniStringComment('Processing using mode:', str_cell);
            str_cell = Ini_Manager.toIniString('p_mode', this.p_mode, str_cell);
            for i = 1 : numel(this.P_MODE)
                str_cell = Ini_Manager.toIniStringComment(sprintf(' %s', this.P_MODE{i}), str_cell);
            end
        end        
    end    
    
    % =========================================================================
    %  LEGACY IMPORT
    % =========================================================================
    methods (Access = 'public')
        function legacyImport(this, state)
            % import from the state variable (saved into the old interface mat file of goGPS)
            try
                this.p_mode = this.gui2mode(state.mode, state.nav_mon, state.kalman_ls, state.code_dd_sa);
            catch ex
                this.logger.addWarning(['Legacy import "Processing mode" failed - ', ex.message])
            end
        end
    end
    
    % =========================================================================
    %  GETTERS
    % =========================================================================
    methods (Access = 'public')
        function [mode, nav_mon, ls_kalman, code_dd_sa] = getGuiMode(this)
            % Get the current mode in GUI format 
            % SYNTAX: [mode, nav_mon, ls_kalman, code_dd_sa] = this.getGuiMode();
            [mode, nav_mon, ls_kalman, code_dd_sa] = this.mode2gui(this.p_mode);
        end
        
        function mode = getMode(this)
            % Get the current mode
            % SYNTAX: mode = this.getMode();
            mode = this.p_mode;
        end
        
    end

    % =========================================================================
    %  LEGACY CONVERTERS
    % =========================================================================
    methods (Static, Access = 'protected')
        function [mode, nav_mon, ls_kalman, code_dd_sa] = mode2gui(p_mode)
            % Convert GUI mode from internal mode - In the GUI mode is splitted in 4 list box, in  goGPS the mode has a unique id
            p_mode = Mode_Settings.P_MODE_2_ID(Mode_Settings.P_MODE_2_ID(:,3) == p_mode, 2);
            mode = floor(p_mode / 1000);
            p_mode = p_mode - mode * 1000;
            nav_mon = floor(p_mode / 100);
            p_mode = p_mode - nav_mon * 100;
            ls_kalman = floor(p_mode / 10);
            p_mode = p_mode - ls_kalman * 10;
            code_dd_sa = p_mode;            
        end
        
        function  p_mode = gui2mode(mode, nav_mon, ls_kalman, code_dd_sa)
            % Convert internal mode from GUI mode - In the GUI mode is splitted in 4 list box, in  goGPS the mode has a unique id
            p_mode = mode * 1000 + (mode == 1) * nav_mon * 100 + (mode > 1) * ls_kalman * 10 + (mode > 1) * code_dd_sa;
            p_mode = Mode_Settings.P_MODE_2_ID(Mode_Settings.P_MODE_2_ID(:,2) == p_mode, 3);
            p_mode = p_mode(1); % this shouldn't be necessary, gui mode are never duplicated
        end
    end
    
    % =========================================================================
    %   MODE FUNCTION (STATIC)
    % =========================================================================
    % function to detect a certain kind of processing
    methods (Static, Access = 'public')
        function isPostProcessing = isPP(mode)
            % return whether or not the mode given in use is a Post Processing mode
            isPostProcessing = sum(intersect(mode, Mode_Settings.GMODE_PP));
        end
        
        function isMonitor = isMON(mode)
            % return whether or not the mode given in use is a Monitor mode
            isMonitor = sum(intersect(mode, Mode_Settings.GMODE_MON));
        end
                
        function isRealTime = isRT(mode)
            % return whether or not the mode given in use is a Real Time mode
            isRealTime = sum(intersect(mode, Mode_Settings.GMODE_RT));
        end
        
        function isDoubleDifferences = isDD(mode)
            % return whether or not the mode given in use is a Double Difference mode
            isDoubleDifferences = sum(intersect(mode, Mode_Settings.GMODE_DD));
        end
        
        function isStandAlone = isSA(mode)
            % return whether or not the mode given in use is a Stand Alone mode
            isStandAlone = sum(intersect(mode, Mode_Settings.GMODE_SA));
        end
        
        function isMultiReceiver = isMR(mode)
            % return whether or not the mode given in use is a Stand Alone mode
            isMultiReceiver = sum(intersect(mode, Mode_Settings.GMODE_MR));
        end
        
        function isUsingPhase = isPH(mode)
            % return whether or not the mode given in use uses Phase
            isUsingPhase = sum(intersect(mode, Mode_Settings.GMODE_PH));
        end
        
         function isKalman = isKM(mode)
            % return whether or not the mode given in use uses Kalman Filter
            isKalman = sum(intersect(mode, Mode_Settings.GMODE_KM));
        end       
        
    end
    
    % =========================================================================
    %  TEST
    % =========================================================================
    methods (Static, Access = 'public')
        function test()      
            % test the class
            s = Mode_Settings();
            s.testInterfaceRoutines();
        end
    end
end
