%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3 beta
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%----------------------------------------------------------------------------------------------
%
%    Code contributed by Andrea Gatti
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
%    along with this program. If not, see <http://www.gnu.org/licenses/>.
%---------------------------------------------------------------------------------------------

classdef goGUIclass < handle

    % =========================================================================
    %   PROPERTIES
    % =========================================================================
    
    properties (Constant)
        % Constellations id => naming inside the goObservation object
        isWin = 0;
        isUnix = 1;
        %isLinux = 2;
        %isMac = 3;
        
        % Colors
        disableCol = [0.502 0.502 0.502];   % Grey (disabled color)
        enableCol = [0 0 0];                % Black (enabled color)
        green = [0 0.8 0];                  % Green - for flag
        yellow = [1 0.8 0.1];               % Yellow - for flag
        red = [1 0 0];                      % Red - for flag
        blue = [0 0 1];                     % Blue - for flag
        
        % goGPS Modes
    end

    properties (GetAccess = 'private', SetAccess = 'private')

    %  HANDLERS
    % =========================================================================
    
        goWB = [];                % waitbar handle
        goh = [];                 % goGPS gui handler
        
        interfaceOS = 0;          % 0 = Windows
                                  % 1 = Linux
                                  % 2 = Mac
        
        edtINI = [];              % Handler to everything related to the editor pf the ini files
                                  
        status; % DEPRECATE: structure containing the state of the parameters of the figure
                                          
    %  INTERFACE STATUS
    % =========================================================================
        
        curState = [];      % This array [n x 1] contains the current status of abilitation of each element of the interface
        newState = [];      % This array [n x 1] contains the future status of abilitation of each element of the interface
        initialState = [];  % This array [n x 1] contains a saved status of abilitation of each element of the interface
        
        getFlag = [];       % Logical flag array [n x 1] contain 1 where a change of the value of the element is needed (from UI to current status)
        setFlag = [];       % Logical flag array [n x 1] contain 1 where a change of the value of the element is needed (from current status to UI)
        curVal = {};        % contains the value of each element of the interface, if one of the element has no value the cell is empty
        newVal = {};        % contains the value of the element of the interface that require to be channged
        
        initialized = 0;
        echoChar = '*';
        
    %  INTERFACE STATUS - PATH
    % =========================================================================

        intSettingsDir = './settings/';                 % Settings folder of goGPS, it contains default settings files
        settingsDir = './settings/';                    % Settings folder of goGPS, it contains all the other settings files
        workingDir = '../data/';                        % Working folder of goGPS, it contains data files/folders
        defaultSettingsFile = 'default_settings.mat';   % Name of the file containing default_settings
        lastSettingsFile = 'last_settings.mat';         % Name of the file containing last settings
        defaultINIFile = 'default_InputFiles.ini';      % Name of the file containing the default example of an ini file (yamatogawa case)
        defaultINIKeywordsFile = 'goGPS_iniDefaultKeywords.ini'; % File containing hint parameters for building an ini files
        
    %  POP UP MENUS
    % ======================================================================
    % Look initPopUps for the complete initialization

        % Processing Mode
        idRealTime = 1;     % id of the pop-up menu Processing mode for Real Time
        idPostProc = 2;     % id of the pop-up menu Processing mode for Post processing
        strMode = {};       % string containing the pop-up menu fields
        
        % Capture mode
        idNav = 1;           % Navigation + capture
        idRMon = 2;          % Rover monitoring (only Rover capture)
        idMMon = 3;          % Master monitoring (only Master capture)
        idRMMon = 4;         % Rover + Master monitoring (Rover and Master capture)
        strCaptureMode = {}; % string containing the pop-up menu fields
        
        % Algorithm type
        idLS = 1;          % Least Squares
        idKF = 2;          % Kalman Filter
        strAlgorithm = {}; % string containing the pop-up menu fields
        
        % Processing type
        idC_SA = 1;      % Code stand-alone (i.e. undifferenced)
        idC_DD = 2;      % Code double difference
        idCP_DD_L = 3;   % Code double difference with LAMBDA
		idCP_Vel = 4;    % Variometric approach for velocity estimation
        idCP_SA = 3;     % Code and phase stand-alone (i.e. undifferenced)    
        idCP_DD = 4;     % Code and phase double difference
        idC_SA_MR = 5;   % Code and phase stand-alone (i.e. undifferenced) for multiple receivers
        idCP_DD_MR = 6;  % Code and phase double difference for multiple receivers
        strTypeLS = {};  % string containing the pop-up menu fields
        strTypeKF = {};  % string containing the pop-up menu fields
        
        % File types
        idRin = 1;
        idBin = 2;
        
        % Integer ambiguity resolution method
        idILS_enum_old = 1;   % ILS method based enumeration in search             (LAMBDA 2.0)
        idILS_shrink = 2;     % ILS method based on search-and-shrink              (LAMBDA 3.0)
        idILS_enum = 3;       % ILS method based enumeration in search             (LAMBDA 3.0)
        idIntRound = 4;       % integer rounding method                            (LAMBDA 3.0)
        idIntBootstr = 5;     % integer bootstrapping method                       (LAMBDA 3.0)
        idPAR = 6;            % PAR with the input P0 of user-defined success rate (LAMBDA 3.0)
        strLAMBDAMethod = {}; % String containing the pop-up menu fields
        
        % Dynamic model
        idCVel = 1;          % Costant velocity
        idCAcc = 2;          % Costant acceleration
        idStatic = 3;        % Static model
        idVariable = 4;      % Variable model
        strDynModel = {};    % String containing the pop-up menu fields
        idMonConstant = 1;   % Monitor mode constant
        idMonVariable = 2;   % Monitor mode variable
        strMonDynModel = {}; % String containing the pop-up menu fields
        
        % Ports
        strPorts = {};       % String containing the pop-up menu fields
        
        % Master station
        idXYZ = 1;           % Pop up Master station (id in the pop-up for idXYZ)
        idGeodetic = 2;      % Pop up Master station (id in the pop-up for idGeodetic)

    %  OTHER IDs
    % ======================================================================
    
        ledOff = 0;          % Led status off
        ledOn  = 1;          % Led status on
        ledOk  = 3;          % Led status ok
        ledKo  = 4;          % Led status ko
        ledCk  = 5;          % Led status check
        ledOp  = 6;          % Led status optional parameter
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
    
    %  ELEMENT IDS
    % ======================================================================
    % Each element of the interface has an id, this is needed to search
    % inside the interface status its values and enable mode
    
        idUI;               % Structure containing all the UI identifiers
        idGroup;            % Structure containing groups of elements        
        
        id2handle;          % Array n x 2 containing at the idUI index the handle of the object in the GUI
    end
    
    % =========================================================================
    %   METHODS
    % =========================================================================
    
    %   CREATOR DESTRUCTOR and INIT
    % -------------------------------------------------------------------------
    % functions to be used to modify the GUI properties    
    methods
        %   CREATOR / DESTRUCTOR
        % =================================================================
        
        % Creator (Brahma)
        function obj = goGUIclass(handles, typeOS)            
            % Creator (Brahma)
            if (length(intersect(typeOS, [obj.isWin obj.isUnix])) ~= 1)
                typeOS = obj.isWin;
            end
            obj.interfaceOS = typeOS;
            
            obj.init(handles);
        end        
        
        function isOS = typeOS(obj, typeOS)
            isOS = obj.interfaceOS == typeOS;
        end
        
        % Destructor (Shiva)
        % function delete(obj); end
    end
    
    % Internal Init functions
    methods(Access = 'private')
        
        % Init all the UI
        function init(obj, handles)
            % Init all the UI
            clear global goIni

            tic;
            obj.goWB = goWaitBar(5, 'Initializing goGPS GUI...');
            obj.goWB.titleUpdate('Init GUI');
            obj.goh = handles;  % Save the handle of the figure
            % Init pup-up strings
            obj.initPopUps();
            
            % Choose default command line output for gui_goGPS_unix
            obj.goh.output = obj.goh.main_panel;
            
            set(obj.goh.main_panel,'CloseRequestFcn',@obj.closeGUI);            
              
            % Update handles structure
            guidata(obj.goh.main_panel, obj.goh);
                
            obj.goWB.goMsg('Centering interface...');
            
            %pixels
            set(obj.goh.main_panel, 'Units', 'pixels' );
            
            %get display size
            screenSize = get(0, 'ScreenSize');
            
            %calculate the center of the display
            position = get(obj.goh.main_panel, 'Position');
            position(1) = (screenSize(3)-position(3))/2;
            position(2) = (screenSize(4)-position(4))/2;
            
            %center the window
            set(obj.goh.main_panel, 'Position', position);

            % Init elements ids
            obj.initInterface();
            obj.goWB.close();
            t0 = toc;
            fprintf('goGPS GUI initialization completed in %.2f seconds\n', t0);
        end

        % Fill all the Pop-up menus
        function initPopUps(obj)
            % Fill all the Pop-up menus

            % Processing Mode
            obj.strMode{obj.idRealTime} = 'Real-time';
            obj.strMode{obj.idPostProc} = 'Post-processing';
                        
            obj.strCaptureMode{obj.idNav} = 'Navigation';
            obj.strCaptureMode{obj.idRMon} = 'Rover Monitor';
            obj.strCaptureMode{obj.idMMon} = 'Master Monitor';
            obj.strCaptureMode{obj.idRMMon} = 'Rover and Master Monitor';

            obj.strAlgorithm{obj.idLS} = 'Least squares';
            obj.strAlgorithm{obj.idKF} = 'Kalman filter';
            
            obj.strTypeLS{obj.idC_SA} = 'Code undifferenced';
            obj.strTypeLS{obj.idC_DD} = 'Code double difference';
            obj.strTypeLS{obj.idCP_DD_L} = 'Code and phase double difference (for LAMBDA)';
			obj.strTypeLS{obj.idCP_Vel} = 'Variometric approach for velocity estimation';
            %obj.strTypeLS{obj.idC_SA_MR} = 'Code undifferenced for multiple receivers';
            %obj.strTypeLS{obj.idCP_DD_MR} = 'Code and phase double difference for multiple receivers';

            obj.strTypeKF{obj.idC_SA} = 'Code undifferenced';
            obj.strTypeKF{obj.idC_DD} = 'Code double difference';
            obj.strTypeKF{obj.idCP_SA} = 'Code and phase undifferenced';
            obj.strTypeKF{obj.idCP_DD} = 'Code and phase double difference';
            %obj.strTypeKF{obj.idCP_DD_MR} = 'Code and phase double difference for multiple receivers';
            
            obj.strLAMBDAMethod{obj.idILS_enum_old} = 'LAMBDA 2.0 - ILS, enumeration';
            obj.strLAMBDAMethod{obj.idILS_shrink}   = 'LAMBDA 3.0 - ILS, search-and-shrink';
            obj.strLAMBDAMethod{obj.idILS_enum}     = 'LAMBDA 3.0 - ILS, enumeration';
            obj.strLAMBDAMethod{obj.idIntRound}     = 'LAMBDA 3.0 - Integer rounding';
            obj.strLAMBDAMethod{obj.idIntBootstr}   = 'LAMBDA 3.0 - Integer bootstrapping';
            obj.strLAMBDAMethod{obj.idPAR}          = 'LAMBDA 3.0 - Partial ambiguity resolution';
            
            obj.strDynModel{obj.idCVel} = 'Const. velocity';
            obj.strDynModel{obj.idCAcc} = 'Const. acceleration';
            obj.strDynModel{obj.idStatic} = 'Static';
            obj.strDynModel{obj.idVariable} = 'Variable';
            
            obj.strMonDynModel{obj.idMonConstant} = 'Constant';
            obj.strMonDynModel{obj.idMonVariable} = 'Variable';
            
            obj.strPorts{1} = '1';
            obj.strPorts{2} = '2';
            obj.strPorts{3} = '3';
            obj.strPorts{4} = '4';
            
            obj.initPorts(obj.strPorts);
        end
        
        % Fill the CaptureMode pop-up (Navigation, Monitor...)
        function initCaptureMode(obj, str)
            % Fill the CaptureMode pop-up (Navigation, Monitor...)
            
            if nargin < 2
                str = obj.strCaptureMode;
            end
            value = get(obj.goh.nav_mon,'Value');
            value = min(1,max(length(str), value));
            set(obj.goh.nav_mon,'Value', value);
            set(obj.goh.nav_mon,'String', str);
        end
        
        % Fill the AlgorithmType pop-up (LS, KF...)
        function initAlgorithmType(obj, str)
            % Fill the AlgorithmType pop-up (LS, KF...)
            if nargin < 2
                str = obj.strAlgorithm;
            end
            value = get(obj.goh.kalman_ls,'Value');
            value = min(1,max(length(str), value));
            set(obj.goh.kalman_ls,'Value', value);
            set(obj.goh.kalman_ls,'String', str);
        end

        % Fill the ProcessingType pop-up (C_SA, CP_DD, ...)
        function initProcessingType(obj, str)
            % Fill the ProcessingType pop-up (C_SA, CP_DD, ...)
            if nargin < 2
                if get(obj.goh.kalman_ls,'Value') == obj.idLS
                    str = obj.strTypeLS;
                else
                    str = obj.strTypeKF;
                end
            end
            value = get(obj.goh.code_dd_sa,'Value');
            value = max(1,min(length(str), value));
            set(obj.goh.code_dd_sa,'Value', value);
            set(obj.goh.code_dd_sa,'String', str);
        end
        
        % Fill the LAMBDA method pop-up (ILS_shrink, ILS_enum, ...)
        function initLAMBDAMethod(obj, str)
            % Fill the LAMBDA method pop-up (ILS_shrink, ILS_enum, ...)
            if nargin < 2
                str = obj.strLAMBDAMethod;
            end
            value = get(obj.goh.lLAMBDAMethod,'Value');
            value = max(1,min(length(str), value));
            set(obj.goh.lLAMBDAMethod,'Value', value);
            set(obj.goh.lLAMBDAMethod,'String', str);
        end
        
        % Fill the Dynamic Model pop-up (constant velocity, acceleration...)
        function initDynModel(obj, str)
            % Fill the Dynamic Model pop-up (constant velocity, acceleration...)        
            if nargin < 2
                str = obj.strDynModel;
            end
            value = get(obj.goh.dyn_mod,'Value');
            value = max(1,min(length(str), value));
            set(obj.goh.dyn_mod, 'String', str);
            set(obj.goh.dyn_mod,'Value', value);
        end
        
        % Dyn model changes wrt the mode
        function resetDynModel(obj)
            if obj.isRealTime()
                switch obj.getElVal(obj.idUI.lCaptMode)
                    case obj.idNav
                        % Full: cv, ca, static, var
                        obj.initDynModel(obj.strDynModel);
                    otherwise
                        % Constant / variable 
                        obj.initDynModel(obj.strMonDynModel);
                end
            else
                switch obj.getElVal(obj.idUI.lProcType)
                    case {obj.idCP_DD, obj.idCP_DD_MR}
                        % Full: cv, ca, static, var
                        obj.initDynModel(obj.strDynModel);
                    otherwise 
                        % cv, ca, static
                        obj.initDynModel(obj.strDynModel(1:3));
                end
            end
        end
        
        % Fill the num ports pop-up list
        function initPorts(obj, str)
            if nargin < 2
                str = obj.strPorts;
            end
            % Set list box with the number of receivers
            value = get(obj.goh.num_receivers,'Value');
            value = min(1,max(length(str), value));
            set(obj.goh.num_receivers,'Value', value);

            % Close old connections that are still open:
            try
                oldConnections = instrfind('Type', 'serial');
                if ~(isempty(oldConnections))
                    fclose(oldConnections);
                end
            catch e
                fprintf('There has been a problem with old serial connections: %s\n', e.message);
            end

            % Detect available ports
            try
                serialInfo = instrhwinfo('serial');
                num_ports = size(serialInfo.AvailableSerialPorts,1);
            catch e
                num_ports = 0;
            end
            
            if num_ports == 0
                set(obj.goh.num_receivers,'String','1');
            elseif num_ports <= size(str,1);
                set(obj.goh.num_receivers,'String',str(1:num_ports));
            else
                set(obj.goh.num_receivers,'String',str);
            end
            
            availablePorts = [];
            try
                availablePorts = serialInfo.AvailableSerialPorts;
            catch e
            end
            if (isempty(availablePorts))
                availablePorts = {'NA'};
            end
            
            set(obj.goh.com_select_0,'String', availablePorts);
            set(obj.goh.com_select_1,'String', availablePorts);
            set(obj.goh.com_select_2,'String', availablePorts);
            set(obj.goh.com_select_3,'String', availablePorts);
        end

        % Fill the available ports pop-ups
        function [select_0 select_1 select_2 select_3] = getPortValues(obj, s0, s1, s2, s3)
            contents = get(obj.goh.com_select_0,'String');
            select_0 = 1; select_1 = 1; select_2 = 1; select_3 = 1;
            for i = 1 : numel(contents)
                if (strcmp(contents{i},s0))
                    select_0 = i;
                elseif (strcmp(contents{i},s1))
                    select_1 = i;
                elseif (strcmp(contents{i},s2))
                    select_2 = i;
                elseif (strcmp(contents{i},s3))
                    select_3 = i;
                end
            end

        end
        
        % Get working directory
        function wdir = getWorkingDir(obj)
            wdir = obj.workingDir;
            if isempty(wdir)
                wdir = './';
            end
            if wdir == 0
                wdir = './';
            end
        end
        
        % Get settings directory
        function sdir = getSettingsDir(obj)
            sdir = obj.intSettingsDir;
            if isempty(sdir)
                sdir = './';
            end
            if sdir == 0
                sdir = './';
            end
        end
    end
    
    %   INTERNAL OBJ STRUCTURE    < = >    GUI
    % -------------------------------------------------------------------------
    % Functions to match the obj internal structure and the goGPS GUI
    % All the functions that could modify directly the object status
    % (acting on obj.goh handle) should be here.
    methods (Access = 'private');
      
        % Set an id value for the each User Interface element
        % This is the most important function containing the entire 
        % configuration of the interface.
        % To add an element in the interface (and to be able to manage it
        % from this object follow these steps:
        %
        %  1. Add an id (with incremental i value)
        %      i=i+1;  id.<name> = i;
        %
        %  2. Associate to this id the correct handle
        %      id2h(i) = obj.goh.<handle>;
        %
        function initUIids(obj)
            
            % Rough pre-allocation
            id2h = zeros(177,1);	% rough estimation of the handle array size
                                    % at the end of the function it will contain the exact size

            % Id naming convections:
            %  p    panels
            %  t    text fields (passive)
            %  s    text fields (active)
            %  n    text fields containing numbers
            %  c    check boxes
            %  r    radio buttons
            %  l    list boxes
            %  f    LED
            %  g    groups
                        
          %   MAIN FIG
          % --------------------------------------------------------------- 
            
            i=1;   id.Fig           = i;    id2h(i) = obj.goh.main_panel;
            
            i=i+1; id.tTitle        = i;    id2h(i) = obj.goh.txtTitle;
            i=i+1; id.tDesc         = i;    id2h(i) = obj.goh.txtDescription;
            
          %   PANELS
          % --------------------------------------------------------------- 
            
            i=i+1; id.pMode         = i;    id2h(i) = obj.goh.uipMode;
            i=i+1; id.pOptions      = i;    id2h(i) = obj.goh.uipOptions;
            i=i+1; id.pConstellations = i;  id2h(i) = obj.goh.pConstellations;
            i=i+1; id.pIntAmb       = i;    id2h(i) = obj.goh.pIntAmb;
            i=i+1; id.pIFiles       = i;    id2h(i) = obj.goh.file_type;
            i=i+1; id.pIOFiles      = i;    id2h(i) = obj.goh.uipInOutFiles;
            i=i+1; id.pSettings     = i;    id2h(i) = obj.goh.uipSettings;
            i=i+1; id.pKF           = i;    id2h(i) = obj.goh.uipKF;
            i=i+1; id.pEStD         = i;    id2h(i) = obj.goh.uipEStD;
            i=i+1; id.pW            = i;    id2h(i) = obj.goh.weight_select;
            i=i+1; id.pDynModel     = i;    id2h(i) = obj.goh.dynamic_model_panel;
            i=i+1; id.pARAA         = i;    id2h(i) = obj.goh.ARAA_panel;
            i=i+1; id.pMSt          = i;    id2h(i) = obj.goh.uipMaster;
            i=i+1; id.pPorts        = i;    id2h(i) = obj.goh.uipPorts;
            i=i+1; id.pMS_NTRIP     = i;    id2h(i) = obj.goh.uipMS_NTRIP;
                        
          %   MODE
          % --------------------------------------------------------------- 
            
            i=i+1; id.lProcMode     = i;    id2h(i) = obj.goh.mode;  
            i=i+1; id.lCaptMode     = i;    id2h(i) = obj.goh.nav_mon;
            i=i+1; id.lAlgType      = i;    id2h(i) = obj.goh.kalman_ls;
            i=i+1; id.lProcType     = i;    id2h(i) = obj.goh.code_dd_sa;
            
            idG.pMode = [id.pMode id.lProcMode:id.lProcType];

          %   INPUT FILE TYPE
          % --------------------------------------------------------------- 

            i=i+1; id.rRin          = i;    id2h(i) = obj.goh.rinex_files;
            i=i+1; id.rBin          = i;    id2h(i) = obj.goh.gogps_data;
            
            % Group of ids in the panel pIFiles
            idG.gIFiles = [id.rRin id.rBin]; 
            idG.pIFiles = [id.pIFiles idG.gIFiles]; 
            
          %   OPTIONS
          % --------------------------------------------------------------- 

            i=i+1; id.cConstraint   = i;    id2h(i) = obj.goh.constraint;
            i=i+1; id.cRefPath      = i;    id2h(i) = obj.goh.ref_path;
            i=i+1; id.cPlotProc     = i;    id2h(i) = obj.goh.plotproc;
            i=i+1; id.cSkyPlot      = i;    id2h(i) = obj.goh.no_skyplot_snr;
            i=i+1; id.cGEarth       = i;    id2h(i) = obj.goh.google_earth;
            i=i+1; id.cErrEllipse   = i;    id2h(i) = obj.goh.err_ellipse;
            i=i+1; id.cPlotMaster   = i;    id2h(i) = obj.goh.plot_master;
            i=i+1; id.cPlotAmb      = i;    id2h(i) = obj.goh.plot_amb;
            i=i+1; id.cUseNTRIP     = i;    id2h(i) = obj.goh.use_ntrip;
            i=i+1; id.cDoppler      = i;    id2h(i) = obj.goh.flag_doppler;
            i=i+1; id.cUse_SBAS     = i;    id2h(i) = obj.goh.use_SBAS;
            
            % Group of ids in the panel pOptions
            idG.gPlotProc = [id.cSkyPlot id.cGEarth id.cErrEllipse id.cPlotMaster id.cPlotAmb];
            idG.pOptions = [id.pOptions id.cConstraint:id.cUse_SBAS];
            
          %   CONSTELLATIONS
          % --------------------------------------------------------------- 
            i=i+1; id.cGPS          = i;    id2h(i) = obj.goh.cGPS;
            i=i+1; id.cGLONASS      = i;    id2h(i) = obj.goh.cGLONASS;
            i=i+1; id.cGalileo      = i;    id2h(i) = obj.goh.cGalileo;
            i=i+1; id.cBeiDou       = i;    id2h(i) = obj.goh.cBeiDou;
            i=i+1; id.cQZSS         = i;    id2h(i) = obj.goh.cQZSS;
            i=i+1; id.cSBAS         = i;    id2h(i) = obj.goh.cSBAS;
            
            % Group of ids in the panel pConstellations
            idG.pGNSS = [id.pConstellations id.cGPS id.cGLONASS id.cGalileo id.cBeiDou id.cQZSS id.cSBAS];
            
            % Constellation of satellites currently supported
            idG.pAvailableGNSSCode = [id.cGPS id.cGLONASS id.cGalileo id.cBeiDou id.cQZSS];
            idG.pAvailableGNSSPhase = [id.cGPS id.cGLONASS id.cGalileo id.cBeiDou id.cQZSS];
            
          %   INTEGER AMBIGUITY RESOLUTION
          % ---------------------------------------------------------------
          
            i=i+1; id.cLAMBDA       = i;    id2h(i) = obj.goh.cLAMBDA;
            i=i+1; id.tLAMBDAMethod = i;    id2h(i) = obj.goh.tLAMBDAMethod;
            i=i+1; id.lLAMBDAMethod = i;    id2h(i) = obj.goh.lLAMBDAMethod;
            i=i+1; id.tP0           = i;    id2h(i) = obj.goh.tP0;
            i=i+1; id.nP0           = i;    id2h(i) = obj.goh.nP0;
            i=i+1; id.cP0           = i;    id2h(i) = obj.goh.cP0;
            i=i+1; id.tMu           = i;    id2h(i) = obj.goh.tMu;
            i=i+1; id.nMu           = i;    id2h(i) = obj.goh.nMu;
            i=i+1; id.cMu           = i;    id2h(i) = obj.goh.cMu;
            
            idG.gLAMBDAMethod = [id.tLAMBDAMethod id.lLAMBDAMethod];
            idG.gLAMBDA3   = [id.tP0 id.nP0 id.cP0 id.cMu];
            idG.gLAMBDAILS = [idG.gLAMBDA3 id.tMu id.nMu];
            idG.gLAMBDA    = [idG.gLAMBDAMethod idG.gLAMBDAILS];
            idG.pIntAmb    = [id.pIntAmb id.cLAMBDA idG.gLAMBDA];

          %   INPUT/OUTPUT FILE AND FOLDERS
          % --------------------------------------------------------------- 

            % INI ---------------------------------------------------------
            i=i+1; id.tINI          = i;    id2h(i) = obj.goh.tINI;
            i=i+1; id.fINI          = i;    id2h(i) = obj.goh.fINI;
            i=i+1; id.sINI          = i;    id2h(i) = obj.goh.sINI;
            i=i+1; id.bINI          = i;    id2h(i) = obj.goh.bINI;
            i=i+1; id.bEditINI      = i;    id2h(i) = obj.goh.bEditINI;
            
            idG.gINI = [id.tINI id.fINI id.sINI id.bINI id.bEditINI];
                        
            % Rover -------------------------------------------------------
            i=i+1; id.tRinRover     = i;    id2h(i) = obj.goh.tRinRover;
            i=i+1; id.fRinRover     = i;    id2h(i) = obj.goh.fRinRover;
            i=i+1; id.tNumRec       = i;    id2h(i) = obj.goh.tNumRec;
            
            idG.RinRover = [id.tRinRover id.tNumRec id.fRinRover];
                        
            % Master ------------------------------------------------------
            i=i+1; id.tRinMaster    = i;    id2h(i) = obj.goh.tRinMaster;
            i=i+1; id.fRinMaster    = i;    id2h(i) = obj.goh.fRinMaster;
            
            idG.RinMaster = [id.tRinMaster id.fRinMaster];
            
            % Navigation --------------------------------------------------
            i=i+1; id.tRinNav       = i;    id2h(i) = obj.goh.tRinNav;
            i=i+1; id.fRinNav       = i;    id2h(i) = obj.goh.fRinNav;
            
            idG.RinNav = [id.tRinNav id.fRinNav ];
            
            % Binary In ---------------------------------------------------
            i=i+1; id.tBinGoIn      = i;    id2h(i) = obj.goh.tBinGoIn;
            i=i+1; id.fBinGoIn      = i;    id2h(i) = obj.goh.fBinGoIn;
            
            idG.BinGoIn = [id.tBinGoIn id.fBinGoIn];            

            % Output ------------------------------------------------------
            i=i+1; id.tDirGoOut     = i;    id2h(i) = obj.goh.tDirGoOut;
            i=i+1; id.fDirGoOut     = i;    id2h(i) = obj.goh.fDirGoOut;
            i=i+1; id.sDirGoOut     = i;    id2h(i) = obj.goh.sDirGoOut;
            i=i+1; id.bDirGoOut     = i;    id2h(i) = obj.goh.bDirGoOut;
            i=i+1; id.tPrefixGoOut  = i;    id2h(i) = obj.goh.tPrefixGoOut;
            i=i+1; id.sPrefixGoOut  = i;    id2h(i) = obj.goh.sPrefixGoOut;
            
            idG.DirGoOut = [id.tDirGoOut id.fDirGoOut id.sDirGoOut id.bDirGoOut];
            idG.GoOut = id.tDirGoOut:id.sPrefixGoOut;

            % DTM ---------------------------------------------------------
            i=i+1; id.tDTM          = i;    id2h(i) = obj.goh.tDTM;
            i=i+1; id.fDTM          = i;    id2h(i) = obj.goh.fDTM;
            
            idG.DTM = [id.tDTM id.fDTM];
            
            % Reference Path ----------------------------------------------
            i=i+1; id.tRefPath      = i;    id2h(i) = obj.goh.tRefPath;
            i=i+1; id.fRefPath      = i;    id2h(i) = obj.goh.fRefPath;
            
            idG.RefPath = [id.tRefPath id.fRefPath];
            
            % PCO/PCV file ------------------------------------------------
            i=i+1; id.tPCO          = i;    id2h(i) = obj.goh.tPCO;
            i=i+1; id.fPCO          = i;    id2h(i) = obj.goh.fPCO;
            
            idG.PCO = [id.tPCO id.fPCO];
            
            % Group of ids in the panel pIOFiles
            idG.pIOFiles = [id.pIOFiles idG.RinRover idG.RinMaster idG.RinNav idG.BinGoIn idG.GoOut idG.DTM idG.RefPath idG.PCO];
            
            % For a correct LED management these following id groups must be synchronized 
            idG.gFileLED = [id.fINI id.fRinRover id.fRinMaster id.fRinNav id.fDTM id.fRefPath id.fPCO];
            idG.gBinLED =  [id.fBinGoIn];
            idG.gInINILED = [id.fRinRover id.fRinMaster id.fRinNav id.fDTM id.fRefPath id.fPCO id.fBinGoIn];
            idG.gDirLED =  [id.fDirGoOut];
            idG.gLED = [idG.gFileLED idG.gBinLED idG.gDirLED];

          %   SETTINGS - KALMAN FILTER - STD
          % --------------------------------------------------------------- 

            % East --------------------------------------------------------            
            i=i+1; id.tStdE         = i;    id2h(i) = obj.goh.text_std_X;
            i=i+1; id.nStdE         = i;    id2h(i) = obj.goh.std_X;
            i=i+1; id.uStdE         = i;    id2h(i) = obj.goh.text_std_X_unit;
            
            idG.StdE = [id.tStdE id.nStdE id.uStdE];
            
            % North -------------------------------------------------------            
            i=i+1; id.tStdN         = i;    id2h(i) = obj.goh.text_std_Y;
            i=i+1; id.nStdN         = i;    id2h(i) = obj.goh.std_Y;
            i=i+1; id.uStdN         = i;    id2h(i) = obj.goh.text_std_Y_unit;
            
            idG.StdN = [id.tStdN id.nStdN id.uStdN];
            
            % Up ----------------------------------------------------------
            i=i+1; id.tStdU         = i;    id2h(i) = obj.goh.text_std_Z;
            i=i+1; id.nStdU         = i;    id2h(i) = obj.goh.std_Z;
            i=i+1; id.uStdU         = i;    id2h(i) = obj.goh.text_std_Z_unit;
            
            idG.StdU = [id.tStdU id.nStdU id.uStdU];
            
            % Code --------------------------------------------------------
            i=i+1; id.tStdCode      = i;    id2h(i) = obj.goh.text_std_code;
            i=i+1; id.nStdCode      = i;    id2h(i) = obj.goh.std_code;
            i=i+1; id.uStdCode      = i;    id2h(i) = obj.goh.text_std_code_unit;
            
            idG.StdCode = [id.tStdCode id.nStdCode id.uStdCode];
            
            % Phase -------------------------------------------------------
            i=i+1; id.bStdPhase     = i;    id2h(i) = obj.goh.toggle_std_phase;
            i=i+1; id.nStdPhase     = i;    id2h(i) = obj.goh.std_phase;
            i=i+1; id.uStdPhase     = i;    id2h(i) = obj.goh.text_std_phase_unit;
            
            idG.StdPhase = [id.bStdPhase id.nStdPhase id.uStdPhase];
            
            % Initial State -----------------------------------------------
            i=i+1; id.tStdT0        = i;    id2h(i) = obj.goh.text_std_init;
            i=i+1; id.nStdT0        = i;    id2h(i) = obj.goh.std_init;
            i=i+1; id.uStdT0        = i;    id2h(i) = obj.goh.text_std_init_unit;
            
            idG.StdT0 = [id.tStdT0 id.nStdT0 id.uStdT0];
            
            % DTM ---------------------------------------------------------
            i=i+1; id.bStdDTM       = i;    id2h(i) = obj.goh.toggle_std_dtm;
            i=i+1; id.nStdDTM       = i;    id2h(i) = obj.goh.std_dtm;
            i=i+1; id.uStdDTM       = i;    id2h(i) = obj.goh.text_std_dtm_unit;
            
            idG.StdDTM = [id.bStdDTM id.nStdDTM id.uStdDTM];
            
            % Velocity ----------------------------------------------------
            i=i+1; id.tStdVel       = i;    id2h(i) = obj.goh.text_std_vel;
            i=i+1; id.nStdVel       = i;    id2h(i) = obj.goh.std_vel;
            i=i+1; id.uStdVel       = i;    id2h(i) = obj.goh.text_std_vel_unit;
            
            idG.StdVel = [id.tStdVel id.nStdVel id.uStdVel];
            
            % Group of ids in the panel pKF
            idG.pKF_ENU = [idG.StdE idG.StdN idG.StdU];
            idG.pKF = [id.pKF idG.StdE idG.StdN idG.StdU idG.StdCode idG.StdPhase idG.StdT0 idG.StdDTM idG.StdVel];
                               
          %   SETTINGS - KALMAN FILTER - WEIGHT MODEL
          % --------------------------------------------------------------- 
             
            i=i+1; id.rW0           = i;    id2h(i) = obj.goh.weight_0;
            i=i+1; id.rW1           = i;    id2h(i) = obj.goh.weight_1;
            i=i+1; id.rW2           = i;    id2h(i) = obj.goh.weight_2;
            i=i+1; id.rW3           = i;    id2h(i) = obj.goh.weight_3;
            i=i+1; id.rW4           = i;    id2h(i) = obj.goh.weight_4;
            
            idG.pW = [id.pW id.rW0 id.rW1 id.rW2 id.rW3 id.rW4];
            
          %   SETTINGS - KALMAN FILTER
          % --------------------------------------------------------------- 

            % Cut-off thr -------------------------------------------------
            i=i+1; id.tCutOff       = i;    id2h(i) = obj.goh.text_cut_off;
            i=i+1; id.nCutOff       = i;    id2h(i) = obj.goh.cut_off;
            i=i+1; id.uCutOff       = i;    id2h(i) = obj.goh.text_cut_off_unit;
            
            idG.CutOff = [id.tCutOff id.nCutOff id.uCutOff];
            
            % SNR thr -----------------------------------------------------
            i=i+1; id.tSNR          = i;    id2h(i) = obj.goh.text_snr_thres;
            i=i+1; id.nSNR          = i;    id2h(i) = obj.goh.snr_thres;
            i=i+1; id.uSNR          = i;    id2h(i) = obj.goh.text_snr_thres_unit;
            
            idG.SNR = [id.tSNR id.nSNR id.uSNR];
            
            % CS thr ------------------------------------------------------
            i=i+1; id.tCS           = i;    id2h(i) = obj.goh.text_cs_thresh;
            i=i+1; id.nCS           = i;    id2h(i) = obj.goh.cs_thresh;
            i=i+1; id.uCS           = i;    id2h(i) = obj.goh.text_cs_thresh_unit;

            idG.CS = [id.tCS id.nCS id.uCS];
            
            % Min sat number ----------------------------------------------
            i=i+1; id.tMinNSat      = i;    id2h(i) = obj.goh.text_min_sat;
            i=i+1; id.nMinNSat      = i;    id2h(i) = obj.goh.min_sat;
            
            idG.MaxNumSat = [id.tMinNSat id.nMinNSat];
            
            % Antenna height ----------------------------------------------
            i=i+1; id.tHAntenna     = i;    id2h(i) = obj.goh.text_antenna_h;
            i=i+1; id.nHAntenna     = i;    id2h(i) = obj.goh.antenna_h;
            i=i+1; id.uHAntenna     = i;    id2h(i) = obj.goh.text_antenna_h_unit;
            
            idG.HAntenna = [id.tHAntenna id.nHAntenna id.uHAntenna];
            
            % Stop Go Stop ------------------------------------------------
            i=i+1; id.cStopGoStop   = i;    id2h(i) = obj.goh.stopGOstop;
            i=i+1; id.tStopGoStop   = i;    id2h(i) = obj.goh.text_stopGOstop;

            idG.StopGoStop = [id.cStopGoStop id.tStopGoStop];
                        
          %   SETTINGS - KALMAN FILTER - DYNAMIC MODEL
          % ---------------------------------------------------------------
            
            i=i+1; id.lDynModel     = i;    id2h(i) = obj.goh.dyn_mod;
            
            idG.pDynModel = [id.pDynModel id.lDynModel];
            
          %   SETTINGS - KALMAN FILTER - ARAA
          % ---------------------------------------------------------------
            
            i=i+1; id.lARAA         = i;    id2h(i) = obj.goh.amb_select;
            
            idG.pARAA = [id.pARAA id.lARAA];
            
          %   SETTINGS - MASTER STATION
          % ---------------------------------------------------------------
            
            i=i+1; id.cMPos         = i;    id2h(i) = obj.goh.master_pos;
            i=i+1; id.lCRS          = i;    id2h(i) = obj.goh.crs;
            
            i=i+1; id.tMX           = i;    id2h(i) = obj.goh.text_master_X;
            i=i+1; id.nMX           = i;    id2h(i) = obj.goh.master_X;
            i=i+1; id.uMX           = i;    id2h(i) = obj.goh.text_master_X_unit;
            i=i+1; id.tMY           = i;    id2h(i) = obj.goh.text_master_Y;
            i=i+1; id.nMY           = i;    id2h(i) = obj.goh.master_Y;
            i=i+1; id.uMY           = i;    id2h(i) = obj.goh.text_master_Y_unit;
            i=i+1; id.tMZ           = i;    id2h(i) = obj.goh.text_master_Z;
            i=i+1; id.nMZ           = i;    id2h(i) = obj.goh.master_Z;
            i=i+1; id.uMZ           = i;    id2h(i) = obj.goh.text_master_Z_unit;
            
            idG.gMXYZ = id.tMX:id.uMZ;
            
            i=i+1; id.tMLat         = i;    id2h(i) = obj.goh.text_master_lat;
            i=i+1; id.nMLat         = i;    id2h(i) = obj.goh.master_lat;
            i=i+1; id.uMLat         = i;    id2h(i) = obj.goh.text_master_lat_unit;
            i=i+1; id.tMLon         = i;    id2h(i) = obj.goh.text_master_lon;
            i=i+1; id.nMLon         = i;    id2h(i) = obj.goh.master_lon;
            i=i+1; id.uMLon         = i;    id2h(i) = obj.goh.text_master_lon_unit;
            i=i+1; id.tMh           = i;    id2h(i) = obj.goh.text_master_h;
            i=i+1; id.nMh           = i;    id2h(i) = obj.goh.master_h;
            i=i+1; id.uMh           = i;    id2h(i) = obj.goh.text_master_h_unit;
            
            idG.gMGeodetic = id.tMLat:id.uMh;
                        
            idG.pMSt = [id.pMSt id.cMPos id.lCRS idG.gMXYZ idG.gMGeodetic];
                        
          %   SETTINGS - PORTS
          % ---------------------------------------------------------------
            
            i=i+1; id.tnPorts       = i;    id2h(i) = obj.goh.text_num_receivers;
            i=i+1; id.lnPorts       = i;    id2h(i) = obj.goh.num_receivers;
            i=i+1; id.lRate         = i;    id2h(i) = obj.goh.pumCaptureRate;            
                    
            idG.nPorts = [id.tnPorts id.lnPorts];
            
            i=i+1; id.lPort0        = i;    id2h(i) = obj.goh.com_select_0;
            i=i+1; id.lProt0        = i;    id2h(i) = obj.goh.protocol_select_0;
            
            idG.lPort0 = [id.lPort0 id.lProt0];
            
            i=i+1; id.lPort1        = i;    id2h(i) = obj.goh.com_select_1;
            i=i+1; id.lProt1        = i;    id2h(i) = obj.goh.protocol_select_1;
            
            idG.lPort1 = [id.lPort1 id.lProt1];
            
            i=i+1; id.lPort2        = i;    id2h(i) = obj.goh.com_select_2;
            i=i+1; id.lProt2        = i;    id2h(i) = obj.goh.protocol_select_2;
            
            idG.lPort2 = [id.lPort2 id.lProt2];
            
            i=i+1; id.lPort3        = i;    id2h(i) = obj.goh.com_select_3;
            i=i+1; id.lProt3        = i;    id2h(i) = obj.goh.protocol_select_3;
            
            idG.lPort3 = [id.lPort3 id.lProt3];
            
            idG.pPorts = [id.pPorts id.tnPorts:id.lPort3];
            
          %   SETTINGS - MASTER SERVER
          % ---------------------------------------------------------------
            
            i=i+1; id.tIPaddr       = i;    id2h(i) = obj.goh.text_IP_address;
            i=i+1; id.sIPaddr       = i;    id2h(i) = obj.goh.IP_address;
            i=i+1; id.tIPport       = i;    id2h(i) = obj.goh.text_port;
            i=i+1; id.sIPport       = i;    id2h(i) = obj.goh.port;
            
            idG.pMS = [id.pMS_NTRIP id.tIPaddr:id.sIPport];
            
            i=i+1; id.tMnt          = i;    id2h(i) = obj.goh.text_mountpoint;
            i=i+1; id.sMnt          = i;    id2h(i) = obj.goh.mountpoint;
            i=i+1; id.tUName        = i;    id2h(i) = obj.goh.text_username;
            i=i+1; id.sUName        = i;    id2h(i) = obj.goh.username;
            i=i+1; id.tUPass        = i;    id2h(i) = obj.goh.text_password;
            i=i+1; id.sUPass        = i;    id2h(i) = obj.goh.password; %
            % Moved in first position => it's a JAVA component
            i=i+1; id.bUPass        = i;    id2h(i) = obj.goh.show_password;

            idG.gNTRIP = id.tMnt:id.bUPass;
            
            i=i+1; id.tApproxPos    = i;    id2h(i) = obj.goh.text_approx_pos;
            i=i+1; id.tVLat         = i;    id2h(i) = obj.goh.text_approx_lat;
            i=i+1; id.nVLat         = i;    id2h(i) = obj.goh.approx_lat;
            i=i+1; id.uVLat         = i;    id2h(i) = obj.goh.text_approx_lat_unit;
            i=i+1; id.tVLon         = i;    id2h(i) = obj.goh.text_approx_lon;
            i=i+1; id.nVLon         = i;    id2h(i) = obj.goh.approx_lon;
            i=i+1; id.uVLon         = i;    id2h(i) = obj.goh.text_approx_lon_unit;
            i=i+1; id.tVH           = i;    id2h(i) = obj.goh.text_approx_h;
            i=i+1; id.nVH           = i;    id2h(i) = obj.goh.approx_h;
            i=i+1; id.uVH           = i;    id2h(i) = obj.goh.text_approx_h_unit;
            
            idG.gApproxPos = id.tApproxPos:id.uVH;
            
            idG.pMS_NTRIP = [idG.pMS idG.gNTRIP idG.gApproxPos];
            
          %   BUTTONS
          % ---------------------------------------------------------------

            i=i+1; id.bExit         = i;    id2h(i) = obj.goh.exit;
            i=i+1; id.bLoad         = i;    id2h(i) = obj.goh.load_button;
            i=i+1; id.bSave         = i;    id2h(i) = obj.goh.save_button;
            i=i+1; id.bGo           = i;    id2h(i) = obj.goh.go_button;                 
          
            idG.gBut = id.bExit:id.bGo;
            
            
            
            idG.Fig = id.pMode:id.bGo;
            idG.ResetStatus = [id.tTitle id.tDesc id.lProcMode id.bExit id.bLoad];
            idG.gDTM = [id.nStdDTM id.uStdDTM idG.DTM idG.HAntenna];
                        
          % ---------------------------------------------------------------
          %   EVENTS GROUPS
          % ---------------------------------------------------------------
          
            % On Real Time
            idG.onRealTime = [idG.ResetStatus id.pConstellations ...
                              id.pMode id.lCaptMode ...
                              id.pIOFiles idG.GoOut ...
                              id.pSettings];
                          
            % On Real Time => Navigation Mode            
            idG.onRT_Nav = [idG.onRealTime idG.pAvailableGNSSPhase ...
                            id.pOptions id.cConstraint id.cRefPath id.cPlotProc id.cUseNTRIP ...
                            id.pKF idG.pKF_ENU idG.StdCode id.bStdPhase idG.StdT0 id.bStdDTM ...
                            idG.pW idG.CutOff idG.SNR idG.CS idG.MaxNumSat idG.StopGoStop ...
                            idG.pDynModel idG.pARAA ...
                            idG.pMSt idG.pMS id.pPorts idG.lPort0];
          
            % On Real Time => Rover Monitor
            idG.onRT_RMon = [idG.onRealTime idG.pAvailableGNSSCode ...
                             id.lRate ...
                             idG.pDynModel ...
                             id.pOptions ...
                             id.pPorts idG.nPorts idG.lPort0];

            % On Real Time => Master Monitor
            idG.onRT_MMon = [idG.onRealTime ...
                             id.pOptions id.cUseNTRIP ...
                             idG.pMS_NTRIP];

            % On Real Time => Master + Rover Monitor
            idG.onRT_RMMon = [idG.onRealTime idG.pAvailableGNSSPhase ...
                              idG.pDynModel ...
                              id.pOptions id.cUseNTRIP ...
                              id.pPorts idG.nPorts idG.lPort0...
                              idG.pMS_NTRIP];
                                                    
          % ---------------------------------------------------------------
            
            % On Post Proc
            idG.onPostProc = [idG.ResetStatus ...
                              id.pMode id.lAlgType id.lProcType ...
                              id.pIFiles idG.gINI...
                              id.pIOFiles id.pConstellations idG.GoOut ...
                              id.pOptions ...
                              id.pSettings idG.CutOff idG.pW];
            
            % On Post Proc => Least Squares
            idG.onPP_LS = [idG.onPostProc  ...
                           id.rBin id.rRin];
                          
            % On Post Proc => Least Squares => Code Stand Alone
            idG.onPP_LS_C_SA = [idG.onPP_LS id.cPlotProc idG.pAvailableGNSSCode ...
                                id.cUse_SBAS];
            
            % On Post Proc => Least Squares => Code Double Differences
            idG.onPP_LS_C_DD = [idG.onPP_LS id.cPlotProc idG.pAvailableGNSSCode ...
                                id.pMSt id.cMPos idG.SNR];
                          
            % On Post Proc => Least Squares => Code and Phase Double Differences
            idG.onPP_LS_CP_DD_L = [idG.onPP_LS id.cPlotProc idG.pAvailableGNSSPhase ...
                                   idG.StdCode idG.StdPhase ...
                                   id.pMSt id.cMPos idG.pIntAmb];

            % On Post Proc => Least Squares => Code and Phase Velocity estimation
            idG.onPP_LS_CP_Vel = [idG.onPP_LS idG.pAvailableGNSSPhase];
            
            % On Post Proc => Least Squares => Code Stand Alone 
            % => Multi Receivers Mode
            idG.onPP_LS_C_SA_MR = [idG.onPP_LS id.cPlotProc idG.pAvailableGNSSCode ...
                                   id.cUse_SBAS];
            
            % On Post Proc => Least Squares => Code and Phase Double 
            % => Multi Receivers Mode
            idG.onPP_LS_CP_DD_MR = [idG.onPP_LS id.cPlotProc idG.pAvailableGNSSPhase ...
                                   idG.StdCode idG.StdPhase ...
                                   id.pMSt id.cMPos idG.pIntAmb];

            % On Post Proc => On Kalman Filter
            idG.onPP_KF = [idG.onPostProc id.cPlotProc ...
                           id.rRin idG.pDynModel ...
                           id.pKF id.pEStD idG.pKF_ENU idG.StdCode idG.StdT0 ...
                           idG.SNR idG.MaxNumSat];
                          
            % On Post Proc => On Kalman Filter => Code Stand Alone
            idG.onPP_KF_C_SA = [idG.onPP_KF id.rBin idG.pAvailableGNSSCode ...
                               id.cUse_SBAS];
            
            % On Post Proc => On Kalman Filter => Code Double Differences
            idG.onPP_KF_C_DD = [idG.onPP_KF id.rBin idG.pAvailableGNSSCode ...
                                id.pMSt id.cMPos];

            % On Post Proc => On Kalman Filter => Code and Phase Stand Alone
            idG.onPP_KF_CP_SA = [idG.onPP_KF id.rBin idG.pAvailableGNSSPhase ...
                                 idG.StdPhase idG.CS ...
                                 id.cDoppler id.cUse_SBAS];
            
            % On Post Proc => On Kalman Filter => Code and Phase Double Differences
            idG.onPP_KF_CP_DD = [idG.onPP_KF id.rBin id.cConstraint idG.pAvailableGNSSPhase ...
                                 idG.StdPhase id.bStdDTM id.cRefPath ...
                                 idG.CS idG.StopGoStop idG.pARAA... 
                                 id.pMSt id.cMPos id.cDoppler idG.pIntAmb];

            % On Post Proc => On Kalman Filter => Code and Phase Double Differences
            % => Multi Receivers Mode
            idG.onPP_KF_CP_DD_MR = [idG.onPP_KF ...
                                    idG.StdPhase ...
                                    idG.CS idG.StopGoStop idG.pARAA... 
                                    id.pMSt id.cMPos id.cDoppler idG.pIntAmb];

            % ---------------------------------------------------------------
          
            
            % On RINEX / BIN
            idG.onRin = [idG.RinRover idG.RinMaster idG.RinNav idG.PCO];
            idG.onBin = [idG.BinGoIn];

            [idG.gPanels idG.strEl idG.valEl] = obj.autoElClassification(id2h);
            
            % Save in object
            obj.idUI = id;
            obj.idGroup = idG;
            obj.id2handle = id2h(1:i);
            
            obj.curState = false(i,1);
            obj.newState = false(i,1);
            obj.setFlag = false(i,1);
            obj.getFlag = false(i,1);
        end
        
        % Create different groups for each type of element
        % Element that are panels
        % Element with the main content stored in:
        %  - 'Value' field
        %  - 'String' field
        function [panel strEl valEl] = autoElClassification(obj,id2h)
            panel = [];
            strEl = [];
            valEl = [];
            for i=1:length(id2h)
                % A panel doesn't have a Style field
                el = get(id2h(i));
                if isfield(el,'Style')                    
                    if sum(strcmp(el.Style, {'edit' 'text'})) > 0
                        strEl = [strEl i];
                    elseif sum(strcmp(el.Style, {'radiobutton', 'checkbox', 'pushbuttons', 'popupmenu'})) > 0
                        valEl = [valEl i];
                    end
                else
                    % A panel doesn't have the field Style
                    % but has the field BorderType
                    if isfield(el,'BorderType')
                        panel = [panel i];
                    end
                end
            end
        end
        
        % Get enable / disable status of the element of the interface
        % This function should be called only once, later in the code
        % the status is kept updated
        function getAllElStatus(obj)            
            panels = false(length(obj.id2handle),1);
            panels(obj.idGroup.gPanels) = true;          % logical indexes of the panels
            idEl = 1:length(obj.id2handle);
            idEl = idEl(~panels);						% id of elements that are not panels
            idEl = idEl(obj.idUI.pMode:end);			% The only elements to be considered starts 
            											% from the principal panel mode
            
            % For each panel
            for i=1:length(obj.idGroup.gPanels)
            	obj.curState(obj.idGroup.gPanels(i)) = obj.isGuiPanelOn(obj.id2handle(obj.idGroup.gPanels(i)));
           	end
           	
           	% For all the other elements
           	for i=1:length(idEl)
            	obj.curState(idEl(i)) = obj.isGuiElOn(obj.id2handle(idEl(i)));
            end
            obj.newState = obj.curState;
        end
        
        % Enable / Disable various elements in the interface
        % the variable newStatus will decide the future status of the
        % interface (function also known as setGuiElStatus)
        function onoffUIEl(obj)
            % Detect modified state (logical array)
            idModified = xor(obj.curState, obj.newState);
            idModified(obj.idUI.Fig) = false; % the figure doesn't change its status;
            
            if (sum(idModified) > 0)
                obj.curState = obj.newState;
                
                % ids of the modified elements
                % state id of the modified element
                listState = uint8(obj.newState)+1;
                state = {'off','on'};
                
                idEl = 1:length(obj.id2handle);              % All the elements
                
                % Sets of panels
                panels = false(length(obj.id2handle),1);     % init logical panels group
                panels(obj.idGroup.gPanels) = true;           % set logical indexes of the panels
                % modified panels
                mPanel = idEl(panels & idModified);
                
                mIdEl = idEl(~panels & idModified);          % id of elements that are not panels that have been modified
                
                % For each modified panel
                for i=1:length(mPanel)
                    obj.onoffGuiPanel(obj.id2handle(mPanel(i)), state{listState(mPanel(i))});
                end
                
                % For all the other modified elements
                for i=1:length(mIdEl)
                    obj.onoffGuiEl(obj.id2handle(mIdEl(i)), state{listState(mIdEl(i))});
                end                
            end
        end

        % Get content of the element of the interface 
        function getAllElContent(obj)
            obj.getFlag(obj.idUI.Fig) = false; % the figure doesn't change its status;
            if sum(obj.getFlag > 0)                
                idEl = 1:length(obj.id2handle);              % All the elements
                
                % Sets of panels
                panels = false(length(obj.id2handle),1);     % init logical panels group
                panels(obj.idGroup.gPanels) = true;           % set logical indexes of the panels
                % modified panels
                mPanel = idEl(panels & obj.getFlag);
                
                % Sets of text elements                
                textEl = false(length(obj.id2handle),1);     % init logical text elements group
                textEl(obj.idGroup.strEl) = true;              % set logical indexes of the text elements
                % Modified LED elements
                mLED = obj.idGroup.gLED(obj.getFlag(obj.idGroup.gLED));
                % Modified text elements
                mTextEl = setdiff(idEl(textEl & obj.getFlag), mLED);
                
                mIdEl = setdiff(idEl(~panels & ~textEl & obj.getFlag), mLED); % id of elements that are not panels nor text elements that have been modified
                
                % For each modified panel
                for i=1:length(mPanel)
                    obj.curVal{mPanel(i)} = obj.getGuiElTitle(obj.id2handle(mPanel(i)));
                end
                
                % For each modified txt Element
                for i=1:length(mTextEl)
                    obj.curVal{mTextEl(i)} = obj.getGuiElStr(obj.id2handle(mTextEl(i)));
                end

                % For all the LEDs
                for i=1:length(mLED)
                    obj.curVal{mLED(i)} = obj.getGuiElColor(obj.id2handle(mLED(i)));
                end

                % For all the other modified elements
                for i=1:length(mIdEl)
                    obj.curVal{mIdEl(i)} = obj.getGuiElVal(obj.id2handle(mIdEl(i)));
                end
                
                % Update newVal if the size of it is different from curVal
                % It may appen in the initialization process
                if sum(size(obj.newVal) == size(obj.curVal)) < 2
                    obj.newVal = obj.curVal;
                end
                obj.newVal(obj.getFlag) = obj.curVal(obj.getFlag);
                obj.getFlag = false(size(obj.getFlag()));
            end

        end

        % Set content of the element of the interface
        function setAllElContent(obj)
            obj.setFlag(1) = false; % the figure doesn't change its status;
            
            if (sum(obj.setFlag) > 0)
                if (obj.setFlag(obj.idUI.sUPass))
                    % Password field is the only one to be managed independently
                    obj.setFlag(obj.idUI.sUPass) = false;
                    obj.showPassword();
                end
                
                idEl = 1:length(obj.id2handle);              % All the elements
                
                % Sets of panels
                panels = false(length(obj.id2handle),1);     % init logical panels group
                panels(obj.idGroup.gPanels) = true;           % set logical indexes of the panels
                % modified panels
                mPanel = idEl(panels & obj.setFlag);
                
                % Sets of text elements
                textEl = false(length(obj.id2handle),1);     % init logical text elements group
                textEl(obj.idGroup.strEl) = true;              % set logical indexes of the text elements
                % Modified LED elements
                mLED = obj.idGroup.gLED(obj.setFlag(obj.idGroup.gLED));
                % Modified text elements
                mTextEl = setdiff(idEl(textEl & obj.setFlag), mLED);
                
                mIdEl = setdiff(idEl(~panels & ~textEl & obj.setFlag), mLED); % id of elements that are not panels nor text elements that have been modified
                
                % For each modified panel
                for i=1:length(mPanel)
                    obj.setGuiElTitle(obj.id2handle(mPanel(i)), obj.newVal{mPanel(i)});
                end
                
                % For each modified txt element
                for i=1:length(mTextEl)
                    obj.setGuiElStr(obj.id2handle(mTextEl(i)), obj.newVal{mTextEl(i)});
                end
                
                % For all the LEDs
                for i=1:length(mLED)
                    obj.setGuiElColor(obj.id2handle(mLED(i)), obj.newVal{mLED(i)});
                end
                
                % For all the other modified elements
                for i=1:length(mIdEl)
                    obj.setGuiElVal(obj.id2handle(mIdEl(i)), obj.newVal{mIdEl(i)});
                end
                
                % Save the status in the object
                obj.curVal(obj.setFlag) = obj.newVal(obj.setFlag);
                obj.setFlag = false(size(obj.setFlag()));
            end
        end        
    end  
    
    %   INTERFACE MANAGEMENT FROM INSIDE THE OBJ
    % -------------------------------------------------------------------------
    % functions to be used to modify the GUI properties
    methods (Access = 'public');
        % Get the status of the interface and prepare the Obj for the
        % management of the GUI
        function initInterface(obj)
            obj.goWB.goMsg('Loading GUI manager object...');

            % Set value for elements ids
            obj.initUIids();            

            % Read interface status (initialize structures
            obj.getAllElStatus();            

            % Create the password field
            obj.initPasswordField('password');
            drawnow;
            
            % If the working folder does not exist
            if isempty(dir(obj.getSettingsDir()))
                waitfor(msgbox('WARNING: The folder containing the settings is not in the working directory, please chose settings directory!'));
                dname = uigetdir('../','Choose the directory that contains goGPS settings');
                if (dname ~= 0)
                    obj.intSettingsDir = [dname '/'];
                    if isempty(dir(obj.getWorkingDir()))
                        obj.workingDir = [obj.getSettingsDir() '../'];
                    end
                end
            end
            
            obj.goWB.goMsg('Loading GUI manager object...');
            obj.goWB.titleUpdate('Import Settings');
            
            % Fill pop up menus
            obj.initPopUp(); % Popup are also modified / reloaded in importStateMatlab
            
            obj.goWB.goMsg('Importing the last used settings...');
            if exist([obj.getSettingsDir() obj.lastSettingsFile],'file')
                obj.importStateMatlab([obj.getSettingsDir() obj.lastSettingsFile]);
            elseif exist([obj.getSettingsDir() obj.defaultSettingsFile],'file')
                obj.importStateMatlab([obj.getSettingsDir() obj.defaultSettingsFile]);
            else
                waitfor(msgbox('No settings file has been found, goGPS may not work properly!'));
            end
            obj.goWB.goMsg('Completing syncing phase...');
            obj.goWB.titleUpdate('Finishing');
            
            % Read interface status as get from file
            obj.initialState = obj.curState;            
            
            % Read current ids content
            obj.getFlag = true(size(obj.curState));
            obj.getAllElContent();
                        
            % Disable interface
            % obj.disableAll();
            obj.checkUIdependencies();

            obj.initialized = 1;
        end

        % Add an undocumented password box
        function initPasswordField(obj, pwd)
            % Undocumented password box for a better management of a password field
            % Create the widget containing the text
            jPwdINI = javax.swing.JPasswordField;
            jPwdINI.setText(pwd);
            % Substitute the eINI edit box with the Java Scroll Pane
            set(obj.goh.uipSettings, 'Units', 'pixels');
            set(obj.goh.uipMS_NTRIP, 'Units', 'pixels');
            set(obj.goh.password, 'Units', 'pixels');
            pos = get(obj.goh.password,'Position');
            posOff = get(obj.goh.uipSettings,'Position');
            pos(1:2) = pos(1:2) + posOff(1:2);
            posOff = get(obj.goh.uipMS_NTRIP,'Position');
            pos(1:2) = pos(1:2) + posOff(1:2);
            pos(1) = pos(1) + 4;
            pos(3) = pos(3) - 4;
            
            [jPwd, hPwd] = javacomponent(jPwdINI, pos, obj.goh.main_panel);
                        
            %obj.id2handle((obj.id2handle==obj.goh.password)) = hPwd;
            %delete(obj.goh.password); % I should delete this but it'll generate some problems
            set(obj.goh.password,'Visible','off');
            obj.goh.jPassword.jpwd = jPwd;
            obj.goh.jPassword.hpwd = hPwd;
            obj.echoChar = obj.goh.jPassword.jpwd.getEchoChar;
            
            obj.setPassword(get(obj.goh.password,'String'));
        end
        
        % Set new enable / disable status
        % Show all the new values stored in the internal state on the GUI
        % Get new values from the GUI
        function updateGUI(obj)
            obj.onoffUIEl();
            obj.getAllElContent();
            obj.checkUIdependencies();
            obj.setAllElContent();
        end
        
        function initPopUp(obj)
            obj.initCaptureMode();
            obj.initAlgorithmType();
            obj.initProcessingType();
            obj.initLAMBDAMethod();
            obj.initDynModel();
        end
        
        % Set the value of an element
        %  - idEl           is the identifier of the object 
        %                   (see idUI for the list of ids)
        %  - value          could be of any type string/boolean/number
        %  - <autoapply>    if set obj.setAllElContent() is automatically
        %                   called after to show immediatly the modification 
        %                   in the GUI
        function setElVal(obj, idEl, value, autoapply)
            if nargin == 3
                autoapply = true;
            end
            obj.setFlag(idEl) = true;
            obj.newVal{idEl} = value;
            if autoapply
                obj.setAllElContent();
            end
        end
        
        % Get the value of an element (from the internal copy of the object)
        function value = getElVal(obj, idEl)
            value = obj.newVal{idEl};
        end
                
        % Set the value of an element
        %  - idEl           is the identifier of the object (or group of
        %                   objects (see idUI/idGroup for the list of ids)
        %  - value          could be 0/1 'on'/'off'
        %  - <autoapply>    if set obj.setAllElContent() is automatically
        %                   called after to show immediatly the modification 
        %                   in the GUI
        function setElStatus(obj, idEl, status, autoapply)
            if nargin == 3
                autoapply = true;
            end
            if ischar(status)
                if strcmp(status,'on')
                    status = true;
                else 
                    status = false;
                end
            end
            obj.newState(idEl) = logical(status);
            if autoapply
                obj.onoffUIEl();
            end
        end
        
        % Return the status of abilitation of the element with id = idEl
        function isOn = isEnabled(obj, idEl)
            isOn = obj.newState(idEl);
        end
        
        % Return 1 if the color of the element with id = idEl is the same
        % of color (This works only for LEDs)
        function isCol = isColor(obj, idEl, color)
            tmp = obj.getElVal(idEl);
            if length(tmp) == 3
                isCol = sum(obj.getElVal(idEl) == color) == 3;
            else
                isCol = false;
            end
        end
        
        % Return the status of activation of the element with id = idEl
        function isOn = isActive(obj, idEl)
            if ischar(obj.getElVal(idEl))
                if strcmp(obj.getElVal(idEl),'on')
                    isOn = true;
                else
                    isOn = false;
                end
            else
                isOn = logical(obj.getElVal(idEl));
            end
            isOn = isOn & obj.isEnabled(idEl);  % To be active an elment must be Enabled and its value = 'on'
        end
        
        function isOk  = okGo(obj, idEl)
            isOk = ~obj.isEnabled(idEl);
            if ~isOk
                isOk = obj.isColor(idEl, obj.green) || obj.isColor(idEl, obj.blue);
            end
        end
                       
        % Disable all the UI elementss but the mode and exit buttons
        function disableAll(obj)
            obj.setElStatus(obj.idGroup.Fig, 0, 0);
            obj.setElStatus(obj.idGroup.ResetStatus, 1, 1);
        end
        
        % TMP: Testing function on the enabler function (onoffEl)
        function testOnOff(obj)
            obj.initialState = obj.curState;            

            obj.getFlag = true(size(obj.curState));
            obj.getAllElContent();

            obj.newState = true(size(obj.curState));
            obj.onoffUIEl(); drawnow;
            for i=1:8
                obj.newState = randn(length(obj.newState),1)>=0;
                obj.onoffUIEl();
                drawnow;
            end
            tic;
            obj.newState = false(size(obj.curState));
            for i=1:length(obj.curState)
               obj.newState(i) = true;
               obj.onoffUIEl();
               drawnow;
            end
            toc;
            
            tic;
            obj.newState = true(size(obj.curState));
            for i=1:length(obj.curState)
               obj.newState(i) = false;
               obj.onoffUIEl();
               drawnow;
            end
            toc;
            obj.newState = false(size(obj.curState));
            obj.newState([obj.idUI.pMode obj.idUI.lProcMode obj.idUI.bExit]) = true;
            obj.onoffUIEl();
            % pause(1);
            
            obj.newState = obj.initialState;
            obj.onoffUIEl();
            obj.setFlag = true(size(obj.curState));
            obj.setAllElContent();
            %pause(0.5);
            %obj.disableAll();
        end
    end      
    
    %   GUI GETTERS
    % -------------------------------------------------------------------------
    % Functions that get specific statuses
    methods
        %   GETTERS
        % =================================================================

        function isI = isInitialized(obj)
            isI = obj.initialized;
        end
            
        % Get the mode as integer value (note that this mode is the global
        % identifier of the goGPS algorithm to be used
        function mode = getgoGPSMode(obj)            
            if obj.isRealTime()
                switch obj.getElVal(obj.idUI.lCaptMode)
                    case obj.idNav
                        mode = goGNSS.MODE_RT_NAV;
                    case obj.idRMon
                        mode = goGNSS.MODE_RT_R_MON;
                    case obj.idMMon
                        mode = goGNSS.MODE_RT_M_MON;
                    case obj.idRMMon
                        mode = goGNSS.MODE_RT_RM_MON;
                end
            elseif obj.isPostProc()
                if obj.isLS()
                    switch obj.getElVal(obj.idUI.lProcType)
                        case obj.idC_SA
                            mode = goGNSS.MODE_PP_LS_C_SA;
                        case obj.idC_DD
                            mode = goGNSS.MODE_PP_LS_C_DD;
                        case obj.idCP_DD_L
                            mode = goGNSS.MODE_PP_LS_CP_DD_L;
						case obj.idCP_Vel
                            mode = goGNSS.MODE_PP_LS_CP_VEL;
                        case obj.idC_SA_MR
                            mode = goGNSS.MODE_PP_LS_C_SA_MR;
                        case obj.idCP_DD_MR
                            mode = goGNSS.MODE_PP_LS_CP_DD_MR;
                    end
                end
                if obj.isKF()
                    switch obj.getElVal(obj.idUI.lProcType)
                        case obj.idC_SA
                            mode = goGNSS.MODE_PP_KF_C_SA;
                        case obj.idC_DD
                            mode = goGNSS.MODE_PP_KF_C_DD;
                        case obj.idCP_SA
                            mode = goGNSS.MODE_PP_KF_CP_SA;
                        case obj.idCP_DD
                            mode = goGNSS.MODE_PP_KF_CP_DD;
                        case obj.idCP_DD_MR
                            mode = goGNSS.MODE_PP_KF_CP_DD_MR;
                    end
                end
            end
        end

        %   INTERFACE GETTERS - MODE
        % =================================================================
        
        % Get processing mode
        function isRT = isRealTime(obj)
            isOn = obj.isEnabled(obj.idUI.lProcMode);
            isRT = isOn && (obj.getElVal(obj.idUI.lProcMode) == obj.idRealTime);
        end
        
        function isPP = isPostProc(obj)
            isOn = obj.isEnabled(obj.idUI.lProcMode);
            isPP = isOn && (obj.getElVal(obj.idUI.lProcMode) == obj.idPostProc);
        end
                
        % Get capture mode
        function isCM = isCaptureMode(obj, idCaptureMode)
            isOn = obj.isEnabled(obj.idUI.lCaptMode);
            isCM = isOn && obj.isRealTime() && (obj.getElVal(obj.idUI.lCaptMode) == idCaptureMode);
        end
        
        % Get algorithm type
        function isLeastSquares = isLS(obj)
            isOn = obj.isEnabled(obj.idUI.lAlgType);
            isLeastSquares = isOn  && obj.isPostProc() && (obj.getElVal(obj.idUI.lAlgType) == obj.idLS);
        end        
        function isKalman = isKF(obj)
            isOn = obj.isEnabled(obj.idUI.lAlgType);
            isKalman = isOn && obj.isPostProc() && (obj.getElVal(obj.idUI.lAlgType) == obj.idKF);
        end

        % Get processing type
        function isPT = isProcessingType(obj, idProcessingType)
            isOn = obj.isEnabled(obj.idUI.lProcType);
            isPT = isOn && obj.isPostProc() && (obj.getElVal(obj.idUI.lProcType) == idProcessingType);
        end
        
        function isSA = isStandAlone(obj)
            isSA = obj.isProcessingType(obj.idC_SA) || (obj.isLS() && obj.isProcessingType(obj.idCP_Vel)) || (obj.isKF() && obj.isProcessingType(obj.idCP_SA));
        end
        
        function isMR = isMultiReceiver(obj)
            isMR = obj.isProcessingType(obj.idC_SA_MR) || obj.isProcessingType(obj.idCP_DD_MR);
        end

        %   INTERFACE GETTERS - INPUT FILE TYPE
        % =================================================================

        % Get file type
        function isRin = isRinex(obj)
            isRin = obj.isActive(obj.idUI.rRin);
        end        
        function isBin = isBin(obj)
            isBin = obj.isActive(obj.idUI.rBin);
        end        
        
        % Set a new password
        function pwd = getPassword(obj)
            if isfield(obj.goh,'jPassword')
                pwd = obj.goh.jPassword.jpwd.getPassword();
            else
                pwd = '';
            end
        end
        
        % Get LAMBDA version
        function isLambda2 = isLambda2(obj)
            isOn = obj.isEnabled(obj.idUI.lLAMBDAMethod);
            isLambda2 = isOn && (obj.getElVal(obj.idUI.lLAMBDAMethod) == obj.idILS_enum_old);
        end
        function isLambda3Par = isLambda3Par(obj)
            isOn = obj.isEnabled(obj.idUI.lLAMBDAMethod);
            isLambda3Par = isOn && (obj.getElVal(obj.idUI.lLAMBDAMethod) == obj.idPAR);
        end
        function isLambdaIls = isLambdaIls(obj)
            isOn = obj.isEnabled(obj.idUI.lLAMBDAMethod);
            isLambdaIls = isOn && ((obj.getElVal(obj.idUI.lLAMBDAMethod) == obj.idILS_enum_old) || ...
                                   (obj.getElVal(obj.idUI.lLAMBDAMethod) == obj.idILS_shrink)   || ...
                                   (obj.getElVal(obj.idUI.lLAMBDAMethod) == obj.idILS_enum));
        end
        
        % Get dynamic model
        function isStatic = isDynModelStatic(obj)
            isOn = obj.isEnabled(obj.idUI.lDynModel);
            isStatic = isOn && (obj.getElVal(obj.idUI.lDynModel) == obj.idStatic);
        end
    end
    
    %   GUI SETTERS
    % -------------------------------------------------------------------------
    % Functions that set specific statuses
    methods 
        % Select one of the two possible file type
        function toggleFileType(obj, value)
            if (nargin == 1)
                if obj.isRinex()
                    value = obj.idBin;
                else
                    value = obj.idRin;
                end
            end
            if (ismember(value, [obj.idRin obj.idBin]))
                obj.setElStatus(obj.idUI.pIFiles, 1);
                if value == obj.idRin
                    obj.setElStatus(obj.idUI.rRin, 1);
                else
                    obj.setElStatus(obj.idUI.rBin, 1);
                end
                obj.setElVal(obj.idUI.rRin, (value == obj.idRin),0);
                obj.setElVal(obj.idUI.rBin,(value == obj.idBin),1);
            end            
        end
        function setRinex(obj)
            obj.toggleFileType(obj.idRin);
        end        
        function setBin(obj)
            obj.toggleFileType(obj.idBin);
        end        
        
        % Set a new password
        function setPassword(obj, password)
            
            if isfield(obj.goh,'jPassword')
                try
                eval(sprintf('obj.goh.jPassword.jpwd.setText(''%s'')',password));
                catch
                end
            elseif (obj.isInitialized())
                set(obj.goh.password,'String',password);
            end
        end
                
        % Modify the password
        function modifyPassword(obj, newkey, newchar)
%             password = obj.getPassword();
%             switch newkey
%                 case 'backspace'
%                     password = password(1:end-1); % Delete the last character in the password
%                 case 'delete'
%                     password = password(2:end); % Delete the first character in the password
%                 otherwise
%                     % If pressed key produces a printable character
%                     if (uint8(newchar) > 32)
%                         password = [password newchar]; % Add the typed character to the password
%                         pause(0.001)%to avoid unwanted character output before the cursor
%                     end
%             end
%             
%             obj.setPassword(password); % Store the password in its current state
        end
        
        % Show the password in the password field
        function showPassword(obj)
            if isfield(obj.goh,'jPassword')
                if obj.isActive(obj.idUI.bUPass)
                    obj.goh.jPassword.jpwd.setEchoChar(char(0));
                else
                    obj.goh.jPassword.jpwd.setEchoChar(obj.echoChar);
                end
            end
        end

    end
    
    %   GUI EVENT MENAGEMENT
    % -------------------------------------------------------------------------
    % Functions that manage events main callers
    methods
        % Function that runs inside onoffUIEl
        % Test every logical dependence in the GUI
        % E.g. a flag that activate other fields
        function checkUIdependencies(obj)

%             % Check Input file type integrity
%             if ~xor(obj.isRinex(), obj.isBin())
%                 obj.setRinex();
%             end
                        
          %   INPUT FILE TYPE
          % --------------------------------------------------------------- 

            % Check File input dependencies
            obj.setElStatus([obj.idGroup.onRin], ~obj.isBin(), 0);
            if obj.isRinex()
                if obj.isStandAlone()
                    obj.setElStatus([obj.idGroup.RinMaster], 0, 0);
                end
            end
            obj.setElStatus([obj.idGroup.onBin], obj.isBin() && ~obj.isRealTime(), 0);
            
          %   OPTIONS
          % --------------------------------------------------------------- 

            % Reference path file
            isOn = obj.isActive(obj.idUI.cRefPath);
            obj.setElStatus([obj.idGroup.RefPath], isOn, 0);
            obj.setElStatus([obj.idUI.cConstraint], isOn, 0);
            
            % Plot while processing flag
            isOn = obj.isActive(obj.idUI.cPlotProc);
            obj.setElStatus([obj.idGroup.gPlotProc], isOn, 0);
            
            % If the master is not available disable the flag to plot it 
            if ~obj.isEnabled(obj.idUI.fRinMaster)
                obj.setElStatus([obj.idUI.cPlotMaster], 0, 0);
            end
            
            % If only code is used, no ambiguities are available
            if obj.isProcessingType(obj.idC_SA) || obj.isProcessingType(obj.idC_DD)
                obj.setElStatus([obj.idUI.cPlotAmb], 0, 0);
            end

            % Only for the variometric approach google Earth plot is disabled
            if obj.isLS() && obj.isProcessingType(obj.idCP_Vel)
                obj.setElStatus([obj.idUI.cGEarth], 0, 0);
            end
            
            % NTRIP flag
            isOn = obj.isActive(obj.idUI.cUseNTRIP);
            obj.setElStatus([obj.idGroup.gNTRIP], isOn, 0);
            
          %   INTEGER AMBIGUITY RESOLUTION
          % --------------------------------------------------------------- 
          
            % LAMBDA flag
            isOn = obj.isActive(obj.idUI.cLAMBDA);
            obj.setElStatus([obj.idGroup.gLAMBDA], isOn, 0);

            % LAMBDA version check
            if (obj.isLambda2())
                obj.setElStatus([obj.idGroup.gLAMBDA3], 0, 0);
            end
            if (~obj.isLambdaIls())
                obj.setElStatus([obj.idGroup.gLAMBDAILS], 0, 0);
                if (obj.isLambda3Par())
                    obj.setElStatus([obj.idUI.tP0], 1, 0);
                    obj.setElStatus([obj.idUI.nP0], 1, 0);
                    obj.setElStatus([obj.idUI.cP0], 1, 0);
                end
            end
            
            if (obj.isLambda3Par())
                set(obj.goh.tP0,'String','Min. success rate (P0):');
            else
                set(obj.goh.tP0,'String','Fixed failure rate (P0):');
            end

            % Automatic mu flag
            isOn = obj.isEnabled(obj.idUI.cMu) && obj.isActive(obj.idUI.cMu);
            if (isOn)
                obj.setElStatus([obj.idUI.nMu], ~isOn, 0);
            end
            
            % Default P0 flag
            isOn = obj.isEnabled(obj.idUI.cP0) && obj.isActive(obj.idUI.cP0);
            if (isOn)
                obj.setElStatus([obj.idUI.nP0], ~isOn, 0);
            end

          %   SETTINGS - KALMAN FILTER - STD
          % --------------------------------------------------------------- 

            % Error standard deviation Phase
            isOn = obj.isActive(obj.idUI.bStdPhase);
            obj.setElStatus([obj.idUI.nStdPhase], isOn, 0);
            obj.setElStatus([obj.idUI.uStdPhase], isOn, 0);
            
            % DTM toggle
            isOn = obj.isActive(obj.idUI.bStdDTM);
            obj.setElStatus([obj.idGroup.gDTM], isOn, 0);

          %   SETTINGS - KALMAN FILTER
          % --------------------------------------------------------------- 

            % Stop Go Stop
            if obj.isEnabled(obj.idUI.cStopGoStop)
                isOn = obj.isActive(obj.idUI.cStopGoStop);
                obj.setElStatus(obj.idGroup.pDynModel, ~isOn, 0);
            end
                        
          %   SETTINGS - PORTS
          % ---------------------------------------------------------------

            % Ports
            if obj.isEnabled(obj.idUI.lnPorts)
                nPorts = obj.getElVal(obj.idUI.lnPorts);
                obj.setElStatus([obj.idGroup.lPort0], nPorts > 0, 0); nPorts = nPorts-1;
                obj.setElStatus([obj.idGroup.lPort1], nPorts > 0, 0); nPorts = nPorts-1;
                obj.setElStatus([obj.idGroup.lPort2], nPorts > 0, 0); nPorts = nPorts-1;
                obj.setElStatus([obj.idGroup.lPort3], nPorts > 0, 0);
            end
            
          %   SETTINGS - KALMAN FILTER - DYNAMIC MODEL
          % ---------------------------------------------------------------

            obj.resetDynModel();
            
            % Check if static dynamic model
            if (obj.isDynModelStatic())
                obj.setElStatus([obj.idGroup.pKF_ENU], 0, 0);
            else
                if (obj.isKF())
                    obj.setElStatus([obj.idGroup.pKF_ENU], 1, 0);
                end
            end
                
          %   SETTINGS - MASTER STATION
          % ---------------------------------------------------------------
            if obj.isEnabled(obj.idUI.cMPos)
                % If the panel is enabled
                isOn = obj.isActive(obj.idUI.cMPos);
                % if isOn the rest should be disabled
                obj.setElStatus(obj.idUI.lCRS, ~isOn, 0);                
                isXYZ = obj.getElVal(obj.idUI.lCRS) == obj.idXYZ;
                obj.setElStatus([obj.idGroup.gMXYZ], isXYZ && ~isOn, 0);
                obj.setElStatus([obj.idGroup.gMGeodetic], ~isXYZ && ~isOn, 0);                
            else
                obj.setElStatus([obj.idGroup.pMSt], 0, 0);
            end

            
          %   SETTINGS - MASTER SERVER
          % ---------------------------------------------------------------
            
            % Password
            obj.showPassword();
            
          %   MODE
          % --------------------------------------------------------------- 
            
          %  % Check list boxes
          %  isOn = obj.isRealTime();
          %  obj.setElStatus([obj.idUI.lCaptMode], isOn, 0);
          %  obj.setElStatus([obj.idUI.lAlgType obj.idUI.lProcType], ~isOn, 1);
            
          %   GO BUTTON AND LEDS
          % --------------------------------------------------------------- 
          % For each file field enabled, I have to check the existence of
          % the folder / file to enable the go Button
            obj.onoffUIEl();
            obj.updateLEDstate();
            goOk = obj.test4Go();
            obj.setElStatus([obj.idUI.bSave obj.idUI.bGo] , goOk, 1);
        end
        
        % Force INI update
        function forceINIupdate(obj)
            obj.updateLEDstate(true);
            goOk = obj.test4Go();
            obj.setElStatus([obj.idUI.bSave obj.idUI.bGo] , goOk, 1);
        end
        
        % EVENT MANAGER
        % When an element is modified (and launch a callback function in
        % the GUI) this function must be called!
        % Sync the internal copy of the interface with the GUI
        function syncFromGUI(obj, idEl)
            if (nargin == 1)
                idEl = obj.idUI.lProcMode;
            end
            obj.getFlag(idEl) = true;
            % Read all the values of the elements
            obj.getAllElContent();
            obj.setAllElContent();

          %   MODE
          % --------------------------------------------------------------- 
            if sum(intersect(idEl, obj.idUI.lProcMode)) > 0
                if obj.isRealTime()
                    obj.setElStatus(obj.idUI.lCaptMode, 1, 0);
                    % Trigger the next popup menu
                    idEl = unique([idEl obj.idUI.lCaptMode]);
                else
                    obj.setElStatus([obj.idUI.lAlgType], 1, 0);
                    % Trigger the next popup menu
                    idEl = unique([idEl obj.idUI.lAlgType]);
                end
            end
            
            % I'm in real time, and this popup menu is active
            if sum(intersect(idEl, obj.idUI.lCaptMode)) > 0
                if obj.isRealTime()
                    obj.setElStatus(obj.idGroup.Fig, 0, 0);
                    switch obj.getElVal(obj.idUI.lCaptMode)
                        case obj.idNav
                            obj.setElStatus(obj.idGroup.onRT_Nav, 1, 0);
                        case obj.idRMon
                            obj.setElStatus(obj.idGroup.onRT_RMon, 1, 0);
                        case obj.idMMon
                            obj.setElStatus(obj.idGroup.onRT_MMon, 1, 0);
                        case obj.idRMMon
                            obj.setElStatus(obj.idGroup.onRT_RMMon, 1, 0);
                    end
                end
                obj.updateGUI();
            end
            
            % I'm in post processing and these popup menus are active
            if sum(intersect(idEl, obj.idUI.lAlgType)) > 0
                obj.initProcessingType();
                obj.getFlag(obj.idUI.lProcType) = true;
                obj.setElStatus([obj.idUI.lProcType], 1, 0);
                % Trigger the next popup menu
                idEl = unique([idEl obj.idUI.lProcType]);
                obj.getAllElContent();
            end        
            
            if sum(intersect(idEl, obj.idUI.lProcType)) > 0
                if obj.isKF() && (obj.getElVal(obj.idUI.lProcType) == obj.idCP_DD_MR)
                    obj.setRinex();
                end

                % Enable / Disable elements
                if obj.isPostProc()
                    if obj.isLS()
                        obj.setElStatus(obj.idGroup.Fig, 0, 0);
                        switch obj.getElVal(obj.idUI.lProcType)
                            case obj.idC_SA
                                obj.setElStatus(obj.idGroup.onPP_LS_C_SA, 1, 0);
                            case obj.idC_DD
                                obj.setElStatus(obj.idGroup.onPP_LS_C_DD, 1, 0);
                            case obj.idCP_DD_L
                                obj.setElStatus(obj.idGroup.onPP_LS_CP_DD_L, 1, 0);
                            case obj.idCP_Vel
                                obj.setElStatus(obj.idGroup.onPP_LS_CP_Vel, 1, 0);
                            case obj.idC_SA_MR
                                obj.setElStatus(obj.idGroup.onPP_LS_C_SA_MR, 1, 0);
                            case obj.idCP_DD_MR
                                obj.setElStatus(obj.idGroup.onPP_LS_CP_DD_MR, 1, 0);
                        end
                    end
                    if obj.isKF()
                        obj.setElStatus(obj.idGroup.Fig, 0, 0);
                        switch obj.getElVal(obj.idUI.lProcType)
                            case obj.idC_SA
                                obj.setElStatus(obj.idGroup.onPP_KF_C_SA, 1, 0);
                            case obj.idC_DD
                                obj.setElStatus(obj.idGroup.onPP_KF_C_DD, 1, 0);
                            case obj.idCP_SA
                                obj.setElStatus(obj.idGroup.onPP_KF_CP_SA, 1, 0);
                            case obj.idCP_DD
                                obj.setElStatus(obj.idGroup.onPP_KF_CP_DD, 1, 0);
                            case obj.idCP_DD_MR                                
                                obj.setElStatus(obj.idGroup.onPP_KF_CP_DD_MR, 1, 0);
                        end
                    end
                    obj.getFlag(idEl) = true;
                    obj.updateGUI();
                end
                
                % Verify if there's still something to enable/disable
                obj.onoffUIEl();
                % Check dependencies
                obj.checkUIdependencies();
            end         
            
          %   INPUT/OUTPUT FILE AND FOLDERS
          % ---------------------------------------------------------------

            % Browse for rover file
            if sum(intersect(idEl, obj.idUI.bINI)) > 0
                obj.browseINIFile();
            end
            if sum(intersect(idEl, obj.idGroup.gINI)) > 0
                obj.forceINIupdate();
            end
                
            % Browse output foder fo binary data
            if sum(intersect(idEl, obj.idUI.bDirGoOut)) > 0
                obj.browseOutDir()
            end

          %   SETTINGS - MASTER STATION
          % ---------------------------------------------------------------

            % Toggle show password
            if sum(intersect(idEl, obj.idUI.bUPass)) > 0
                obj.showPassword();
            end

            
          %   BUTTONS
          % ---------------------------------------------------------------
                        
            % on Load Settings
            if sum(intersect(idEl, obj.idUI.bLoad)) > 0
                obj.loadState();
            end
            
            % on Save Settings
            if sum(intersect(idEl, obj.idUI.bSave)) > 0
                obj.saveState();
            end
            
            % on GO
            if sum(intersect(idEl, obj.idUI.bGo)) > 0
                obj.go();
            end
            
            obj.onoffUIEl();
            obj.checkUIdependencies();
            
            % on Exit
            if sum(intersect(idEl, obj.idUI.bExit)) > 0
                if isfield(obj.goh, 'main_panel');
                    close(obj.goh.main_panel);
                else % Something went wrong in the initialization phase
                    close(gcf);
                end
            end
        end
    end
    
    %   GUI EVENT MENAGEMENT - AUXILIAR FUNCTIONS
    % -------------------------------------------------------------------------
    % Functions that manage events - auxiliar functions
    methods
        % Test if the active file/dir paths
        % contain valid file/dir
        function updateLEDstate(obj, force)
            global goIni
            if nargin == 1
                force = 0;      % force file update
            end            
            
            % Check INI file
            if obj.isEnabled(obj.idUI.sINI);
                filename = obj.getElVal(obj.idUI.sINI);
                if isempty(filename)
                    obj.setGUILedStatus(obj.idUI.fINI, obj.ledKo, 0);
                    % If I do not have an INI file, every LED should be RED
                    for  i = 1:length(obj.idGroup.gInINILED)
                        obj.setGUILedStatus(obj.idGroup.gInINILED(i), obj.ledKo, 0);
                    end
                else
                    if exist(filename,'file');
                        obj.setGUILedStatus(obj.idUI.fINI, obj.ledOk, 0);
                        
                        if obj.isPostProc()
                            % If needed init INI reader
                            if isempty(goIni)
                                goIni = goIniReader('', 0);                                
                            end
                            if (~goIni.getReadStatus())
                                goIni.readFile();
                            end
                            % If I have to update the ini file
                            goIni.update(filename, force);
                            % Receivers file --------------------------------------
                            nR = goIni.getData('Receivers','nRec');
                            data_path = goIni.getData('Receivers','data_path');
                            file_name = goIni.getData('Receivers','file_name');
                            
                            if (isempty(nR))
                                if iscell(file_name)
                                    nR = length(file_name);
                                else
                                    nR = 1;
                                end
                                goIni.addKey('Receivers','nRec',nR);
                            end
                            obj.setElVal(obj.idUI.tNumRec,['x ' num2str(nR)]);
                            
                            if (isempty(data_path))
                                data_path = '';
                            end
                            if (isempty(file_name))
                                obj.setGUILedStatus(obj.idUI.fRinRover, obj.ledKo, 0);
                            else
                                % If I have more than one receiver
                                if iscell(file_name)
                                    % The number of receiver is = to the number of files?
                                    if nR ~= length((file_name))   % Declared number of file ~= number of files
                                        obj.setGUILedStatus(obj.idUI.fRinRover, obj.ledCk, 0);
                                    else
                                        % Check the presence of all the files
                                        fileOk = true;
                                        for r = 1:nR
                                            if ~exist([data_path file_name{r}],'file')
                                                fileOk = false;
                                            end
                                        end
                                        if fileOk
                                            obj.setGUILedStatus(obj.idUI.fRinRover, obj.ledOk, 0);
                                        else
                                            obj.setGUILedStatus(obj.idUI.fRinRover, obj.ledCk, 0);
                                        end
                                    end
                                else
                                    % The number of receiver is = to the number of files?
                                    if nR > 1   % Declared number of file ~= number of files
                                        obj.setGUILedStatus(obj.idUI.fRinRover, obj.ledCk, 0);
                                    else
                                        % Check the presence of all the files
                                        if exist([data_path file_name],'file')
                                            obj.setGUILedStatus(obj.idUI.fRinRover, obj.ledOk, 0);
                                        else
                                            obj.setGUILedStatus(obj.idUI.fRinRover, obj.ledCk, 0);
                                        end
                                    end
                                end
                            end
                            
                            % Master file -----------------------------------------
                            data_path = goIni.getData('Master','data_path');
                            file_name = goIni.getData('Master','file_name');
                            if (isempty(data_path))
                                data_path = '';
                            end
                            if (isempty(file_name))
                                obj.setGUILedStatus(obj.idUI.fRinMaster, obj.ledKo, 0);
                            else
                                % Check the presence of all the files
                                if exist([data_path file_name],'file')
                                    obj.setGUILedStatus(obj.idUI.fRinMaster, obj.ledOk, 0);
                                else
                                    obj.setGUILedStatus(obj.idUI.fRinMaster, obj.ledCk, 0);
                                end
                            end
                            
                            % Navigation file -------------------------------------
                            data_path = goIni.getData('Navigational','data_path');
                            file_name = goIni.getData('Navigational','file_name');
                            if (isempty(data_path))
                                data_path = '';
                            end
                            if (isempty(file_name))
                                obj.setGUILedStatus(obj.idUI.fRinNav, obj.ledKo, 0);
                            else
                                % Check the presence of all the files
                                if exist([data_path file_name],'file')
                                    obj.setGUILedStatus(obj.idUI.fRinNav, obj.ledOk, 0);
                                else
                                    obj.setGUILedStatus(obj.idUI.fRinNav, obj.ledCk, 0);
                                end
                            end
                            
                            % Bin file --------------------------------------------
                            data_path = goIni.getData('Bin','data_path');
                            file_prefix = goIni.getData('Bin','file_prefix');
                            if (isempty(data_path))
                                data_path = '';
                            end
                            if (isempty(file_prefix))
                                obj.setGUILedStatus(obj.idUI.fBinGoIn, obj.ledKo, 0);
                            else
                                % Check the presence of the files
                                if ~(isempty(dir([data_path file_prefix '_obs*.bin'])) || isempty(dir([data_path file_prefix '_eph*.bin'])))
                                    obj.setGUILedStatus(obj.idUI.fBinGoIn, obj.ledOk, 0);
                                else
                                    obj.setGUILedStatus(obj.idUI.fBinGoIn, obj.ledCk, 0);
                                end
                            end
                            
                            % DTM file -----------------------------------------------
                            data_path = goIni.getData('DTM','data_path');
                            if (isempty(data_path))
                                obj.setGUILedStatus(obj.idUI.fDTM, obj.ledKo, 0);
                            else
                                % Check the presence of the directory
                                if exist(data_path,'dir')
                                    obj.setGUILedStatus(obj.idUI.fDTM, obj.ledOk, 0);
                                else
                                    obj.setGUILedStatus(obj.idUI.fDTM, obj.ledCk, 0);
                                end
                            end
                            
                            % Reference path file ------------------------------------
                            data_path = goIni.getData('RefPath','data_path');
                            file_name = goIni.getData('RefPath','file_name');
                            if (isempty(data_path))
                                data_path = '';
                            end
                            if (isempty(file_name))
                                obj.setGUILedStatus(obj.idUI.fRefPath, obj.ledKo, 0);
                            else
                                % Check the presence of all the files
                                if exist([data_path file_name],'file')
                                    obj.setGUILedStatus(obj.idUI.fRefPath, obj.ledOk, 0);
                                else
                                    obj.setGUILedStatus(obj.idUI.fRefPath, obj.ledCk, 0);
                                end
                            end
                            
                            % PCO/PCV file ------------------------------------
                            data_path = goIni.getData('PCO_PCV_file','data_path');
                            file_name = goIni.getData('PCO_PCV_file','file_name');
                            if (isempty(data_path))
                                data_path = '';
                            end
                            if (isempty(file_name))
                                obj.setGUILedStatus(obj.idUI.fPCO, obj.ledOp, 0);
                            else
                                % Check the presence of all the files
                                if exist([data_path file_name],'file')
                                    obj.setGUILedStatus(obj.idUI.fPCO, obj.ledOk, 0);
                                else
                                    obj.setGUILedStatus(obj.idUI.fPCO, obj.ledCk, 0);
                                end
                            end
                        end
                    else
                        obj.setGUILedStatus(obj.idUI.fINI, obj.ledCk, 0);
                        % If I do not have an INI file, every LED should be RED
                        for  i = 1:length(obj.idGroup.gInINILED)
                            obj.setGUILedStatus(obj.idGroup.gInINILED(i), obj.ledKo, 0);
                        end
                    end
                end
            end
            
          % Output dir --------------------------------------------------------
            
          outDir = obj.getElVal(obj.idUI.sDirGoOut);
          if isempty(outDir)
              obj.setGUILedStatus(obj.idUI.fDirGoOut, obj.ledKo, 0);
          else
              if ~exist(outDir,'dir');
                  obj.setGUILedStatus(obj.idUI.fDirGoOut, obj.ledCk, 0);
              else
                  obj.setGUILedStatus(obj.idUI.fDirGoOut, obj.ledOk, 0);
              end
          end
            
        end        
        
        % Test if all the folder / files are ok
        % and it is possible to activate go and save buttons
        function goOk = test4Go(obj)
            goOk = 0;
            % Performs all the led check before allowing the launch... ehm
            % activation of the Go! button
            if obj.isPostProc()
                for i = 1:length(obj.idGroup.gLED)
                    goOk = goOk+obj.okGo(obj.idGroup.gLED(i));
                end
                goOk = (goOk - length(obj.idGroup.gLED)) == 0;
            else % In real time I just have to check for the output folder
                goOk = obj.okGo(obj.idUI.fDirGoOut) == 1;
            end
            
            % Chek constellations to be used
            % I need at least one of these constellation active
            activeGNSS        = [obj.isActive(obj.idUI.cGPS) ...
                                 obj.isActive(obj.idUI.cGLONASS) ...
                                 obj.isActive(obj.idUI.cGalileo) ...
                                 obj.isActive(obj.idUI.cBeiDou) ...
                                 obj.isActive(obj.idUI.cQZSS)];
            if sum(activeGNSS) == 0
                obj.setElVal(obj.idUI.cGPS, true, 0);
                % goOk = 0;
            end
        end
        
        % Browse INI file
        function browseINIFile(obj)
            % In multi receiver mode, I read from ini file
            [filename, pathname] = uigetfile( ...
                {'*.ini;','INI configuration file (*.ini)'; ...
                '*.*',  'All Files (*.*)'}, ...
                'Choose an INI configuration file',[obj.getSettingsDir()]);
            if (filename ~= 0)
                obj.setElVal(obj.idUI.sINI, fullfile(pathname, filename));
            end
            obj.updateGUI();
        end
        
        % Browse for rover file
        function browseRoverObsFile(obj)
            if obj.isPostProc() && (obj.getElVal(obj.idUI.lProcType) == obj.idCP_DD_MR)
                % In multi receiver mode, I read from ini file
                [filename, pathname] = uigetfile( ...
                    {'*.ini;','INI configuration file (*.ini)'; ...
                    '*.*',  'All Files (*.*)'}, ...
                    'Choose an INI configuration file',[obj.getSettingsDir()]);                
            else
                [filename, pathname] = uigetfile( ...
                    {'*.obs;*.??o;*.??O','RINEX observation files (*.obs,*.??o,*.??O)';
                    '*.obs','Observation files (*.obs)'; ...
                    '*.??o;*.??O','Observation files (*.??o,*.??O)'; ...
                    '*.*',  'All Files (*.*)'}, ...
                    'Choose a RINEX observation file for the rover',[obj.getWorkingDir() 'data_RINEX']);
            end
            if (filename ~= 0)
                obj.setElVal(obj.idUI.sRinRover, fullfile(pathname, filename));
            end
            obj.updateGUI();
        end
        
        % Browse for a master file
        function browseMasterObsFile(obj)
            [filename, pathname] = uigetfile( ...
                {'*.obs;*.??o','RINEX observation files (*.obs,*.??o)';
                '*.obs','Observation files (*.obs)'; ...
                '*.??o','Observation files (*.??o)'; ...
                '*.*',  'All Files (*.*)'}, ...
                'Choose a RINEX observation file for the master',[obj.getWorkingDir() 'data_RINEX']);
            
            if (filename ~= 0)
                obj.setElVal(obj.idUI.sRinMaster, fullfile(pathname, filename));
            end
            obj.updateGUI();
        end
        
        % Browse for a navigation file
        function browseNavigationFile(obj)
            [filename, pathname] = uigetfile( ...
                {'*.nav;*.??n;*.??N','RINEX navigation files (*.nav,*.??n,*.??N)';
                '*.nav','Navigation files (*.nav)'; ...
                '*.??n;*.??N','Navigation files (*.??n,*.??N)'; ...
                '*.*',  'All Files (*.*)'}, ...
                'Choose a RINEX navigation file',[obj.getWorkingDir() 'data_RINEX']);
            
            if (filename ~= 0)
                obj.setElVal(obj.idUI.sRinNav, fullfile(pathname, filename));
            end
            obj.updateGUI();
        end
        
        % Browse for a binary file
        function browseBinaryFile(obj)
            [filename, pathname] = uigetfile( ...
                {'*.bin','goGPS binary data (*.bin)'}, ...
                'Choose goGPS binary data','../data');
            
            if (filename ~= 0)
                pos = find(filename == '_');
                filename = filename(1:pos(end-1)-1);
                obj.setElVal(obj.idUI.sBinGoIn, fullfile(pathname, filename));
            end
            obj.updateGUI();            
        end
        
        % Browse output foder fo binary data
        function browseOutDir(obj)
            dname = uigetdir(obj.getWorkingDir(),'Choose a directory to store goGPS data');
            if (dname ~= 0)
                obj.setElVal(obj.idUI.sDirGoOut, dname);
            end
            obj.updateGUI();
        end
        
        % Browse for a DTM folder
        function browseDtmDir(obj)
            dname = uigetdir([obj.getWorkingDir() 'dtm'],'Choose a directory containing DTM data');
            if (dname ~= 0)
                obj.setElVal(obj.idUI.sDTM, dname);
            end
            obj.updateGUI();
        end
        
        % Browse for the path containing reference points for constrained solutions
        function browseRefFile(obj)
            [filename, pathname] = uigetfile('*.mat', 'Choose file containing reference path','../data');
            
            if (filename ~= 0)
                obj.setElVal(obj.idUI.sRefPath, fullfile(pathname, filename));
            end
            obj.updateGUI();
        end
        
        % Led status
        % Set green / red status of the UI and optionally lock the UI
        % -------------------------------------------------------------------------
        function setGUILedStatus(obj, idEl, status, autoapply)
            if nargin == 3
                autoapply = true;
            end
            
            if (status == obj.ledOk)    % Led Ok
                if ~obj.isColor(idEl, obj.green)
                    obj.setElVal(idEl,obj.green)
                end
            elseif (status == obj.ledKo) % Led Ko
                if ~obj.isColor(idEl, obj.red)
                    obj.setElVal(idEl,obj.red)
                end
            elseif (status == obj.ledCk) % Led Check
                if ~obj.isColor(idEl, obj.yellow)
                    obj.setElVal(idEl,obj.yellow)
                end
            elseif (status == obj.ledOp) % Led Optional parameter
                if ~obj.isColor(idEl, obj.blue)
                    obj.setElVal(idEl,obj.blue)
                end
            end
            obj.setAllElContent();
            %obj.setElStatus(idEl, status > 1, autoapply)
        end
        
    end
    
    %   IMPORT /EXPORT SETTINGS
    % -------------------------------------------------------------------------
    % Functions to load save settings from file        
    methods
        function loadState(obj)
            [filename, pathname] = uigetfile('*.mat', 'Choose file with saved settings',obj.settingsDir);
            
            if pathname == 0 %if the user pressed cancelled, then we exit this callback
                return
            end
            %construct the path name of the file to be loaded
            loadDataName = fullfile(pathname,filename);
            
            %load the settings, which creates a new gui
            obj.importStateMatlab(loadDataName);
        end
        
        function saveState(obj)
            [filename,pathname] = uiputfile('*.mat','Save your GUI settings',obj.settingsDir);
            
            if pathname == 0 %if the user pressed cancelled, then we exit this callback
                return
            end
            %construct the path name of the save location
            saveDataName = fullfile(pathname,filename);
            
            %saves the gui data
            obj.exportStateMatlab(saveDataName);
        end
        
        % Load the state of the gui from a matlab file.
        function importStateMatlab(obj,filename)
            load(filename); % the file contains the variable state
            obj.status = state;
            
            %   MODE
            % ===============================================================
            
            obj.setElVal(obj.idUI.lProcMode, state.mode, 0);
            
            obj.initCaptureMode();
            obj.setElVal(obj.idUI.lCaptMode, state.nav_mon, 0);
            obj.initAlgorithmType();
            obj.setElVal(obj.idUI.lAlgType, state.kalman_ls, 0);
            obj.initProcessingType();
            obj.setElVal(obj.idUI.lProcType, state.code_dd_sa, 0);
            
            %   INPUT FILE TYPE
            % ===============================================================
            
            obj.setElVal(obj.idUI.rRin, state.rinex_files, 0);
            obj.setElVal(obj.idUI.rBin, state.gogps_data, 1);
            
            %   OPTIONS
            % ===============================================================
            
            obj.setElVal(obj.idUI.cMPos, state.master_pos, 0);
            obj.setElVal(obj.idUI.cConstraint, state.constraint, 0);
            obj.setElVal(obj.idUI.cPlotProc, state.plotproc, 0);
            obj.setElVal(obj.idUI.cRefPath, state.ref_path, 0);
            obj.setElVal(obj.idUI.cSkyPlot, state.no_skyplot_snr, 0);
            obj.setElVal(obj.idUI.cGEarth, state.google_earth, 0);
            obj.setElVal(obj.idUI.cErrEllipse, state.err_ellipse, 0);
            obj.setElVal(obj.idUI.cPlotMaster, state.plot_master, 0);
            obj.setElVal(obj.idUI.cPlotAmb, state.plot_amb, 0);
            obj.setElVal(obj.idUI.cUseNTRIP, state.use_ntrip, 0);
            obj.setElVal(obj.idUI.cDoppler, state.flag_doppler, 0);
            % Temporary check in the migration to SBAS use period
            if (isfield(state,'use_sbas')) %since v0.3.2beta -> backward compatibility
                obj.setElVal(obj.idUI.cUse_SBAS, state.use_sbas, 0);
            else
                obj.setElVal(obj.idUI.cUse_SBAS, 0, 0);
            end
            
            %   INTEGER AMBIGUITY RESOLUTION
            % ===============================================================
            
            obj.setElVal(obj.idUI.cLAMBDA, state.use_lambda, 0);
            obj.setElVal(obj.idUI.lLAMBDAMethod, state.lambda_method, 0);
            obj.setElVal(obj.idUI.nP0, state.lambda_P0, 0);
            obj.setElVal(obj.idUI.cP0, state.lambda_default_P0, 0);
            obj.setElVal(obj.idUI.nMu, state.lambda_mu, 0);
            obj.setElVal(obj.idUI.cMu, state.lambda_auto_mu, 0);
            
            %   INPUT/OUTPUT FILE AND FOLDERS
            % ===============================================================
            
            if (isfield(state,'INIsettings')) % backward compatibility check
                obj.setElVal(obj.idUI.sINI, state.INIsettings, 0);
            end
            obj.setElVal(obj.idUI.sDirGoOut, state.gogps_data_output, 0);
            obj.setElVal(obj.idUI.sPrefixGoOut, state.gogps_data_output_prefix, 0);
            if (isfield(state,'activeGNSS'))
                obj.setElVal(obj.idUI.cGPS, state.activeGNSS(1), 0);
                obj.setElVal(obj.idUI.cGLONASS, state.activeGNSS(2), 0);
                obj.setElVal(obj.idUI.cGalileo, state.activeGNSS(3), 0);
                obj.setElVal(obj.idUI.cBeiDou, state.activeGNSS(4), 0);
                obj.setElVal(obj.idUI.cQZSS, state.activeGNSS(5), 0);
                obj.setElVal(obj.idUI.cSBAS, state.activeGNSS(6), 0);
            end
            
            %   SETTINGS - MASTER STATION
            % ===============================================================
            
            obj.setElVal(obj.idUI.lCRS, state.crs, 0);
            
            obj.setElVal(obj.idUI.nMX, state.master_X, 0);
            obj.setElVal(obj.idUI.nMY, state.master_Y, 0);
            obj.setElVal(obj.idUI.nMZ, state.master_Z, 0);
            obj.setElVal(obj.idUI.nMLat, state.master_lat, 0);
            obj.setElVal(obj.idUI.nMLon, state.master_lon, 0);
            obj.setElVal(obj.idUI.nMh, state.master_h, 0);
            
            %   SETTINGS - KALMAN FILTER - STD
            % ===============================================================
            
            obj.setElVal(obj.idUI.nStdE, state.std_X, 0);
            obj.setElVal(obj.idUI.nStdN, state.std_Y, 0);
            obj.setElVal(obj.idUI.nStdU, state.std_Z, 0);
            obj.setElVal(obj.idUI.nStdCode, state.std_code, 0);
            obj.setElVal(obj.idUI.nStdPhase, state.std_phase, 0);
            obj.setElVal(obj.idUI.nStdDTM, state.std_dtm, 0);
            obj.setElVal(obj.idUI.bStdPhase, state.toggle_std_phase, 0);
            obj.setElVal(obj.idUI.bStdDTM, state.toggle_std_dtm, 0);
            obj.setElVal(obj.idUI.nStdT0, state.std_init, 0);
            obj.setElVal(obj.idUI.nStdVel, state.std_vel, 0);
            
            %   SETTINGS - KALMAN FILTER - WEIGHT MODEL
            % ===============================================================
            
            obj.setElVal(obj.idUI.rW0, state.weight_0, 0);
            obj.setElVal(obj.idUI.rW1, state.weight_1, 0);
            obj.setElVal(obj.idUI.rW2, state.weight_2, 0);
            obj.setElVal(obj.idUI.rW3, state.weight_3, 0);
            % Temporary check during development
            if (isfield(state,'weight_4')) %since v0.3.2beta -> backward compatibility
                obj.setElVal(obj.idUI.rW4, state.weight_4, 0);
            else
                obj.setElVal(obj.idUI.rW4, 0, 0);
            end
                        
            %   SETTINGS - KALMAN FILTER
            % ===============================================================
            
            obj.setElVal(obj.idUI.nCS, state.cs_thresh, 0);
            obj.setElVal(obj.idUI.nCutOff, state.cut_off, 0);
            obj.setElVal(obj.idUI.nSNR, state.snr_thres, 0);
            obj.setElVal(obj.idUI.nHAntenna, state.antenna_h, 0);
            obj.setElVal(obj.idUI.nMinNSat, state.min_sat, 0);
            obj.setElVal(obj.idUI.cStopGoStop, state.stopGOstop, 0);
            
            %   SETTINGS - KALMAN FILTER - DYNAMIC MODEL
            % ===============================================================
            
            obj.resetDynModel();
            obj.setElVal(obj.idUI.lDynModel, state.dyn_mod, 0);
            
            %   SETTINGS - KALMAN FILTER - ARAA
            % ===============================================================
            
            obj.setElVal(obj.idUI.lARAA, state.amb_select, 0);
            
            %   SETTINGS - PORTS
            % ===============================================================
            
            if (isfield(state,'captureRate'))
                obj.setElVal(obj.idUI.lRate, state.captureRate, 0);
            end
            [s0 s1 s2 s3] = obj.getPortValues(state.com_select_0, state.com_select_1, state.com_select_2, state.com_select_3);
            
            obj.setElVal(obj.idUI.lnPorts, state.num_receivers, 0);
            
            obj.setElVal(obj.idUI.lPort0, s0, 0);
            obj.setElVal(obj.idUI.lPort1, s1, 0);
            obj.setElVal(obj.idUI.lPort2, s2, 0);
            obj.setElVal(obj.idUI.lPort3, s3, 0);
            obj.setElVal(obj.idUI.lProt0, state.protocol_select_0, 0);
            obj.setElVal(obj.idUI.lProt1, state.protocol_select_1, 0);
            obj.setElVal(obj.idUI.lProt2, state.protocol_select_2, 0);
            obj.setElVal(obj.idUI.lProt3, state.protocol_select_3, 0);
            
            %   SETTINGS - MASTER SERVER
            % ===============================================================
            
            obj.setElVal(obj.idUI.sIPaddr, state.IP_address, 0);
            obj.setElVal(obj.idUI.sIPport, state.port, 0);
            obj.setElVal(obj.idUI.sMnt, state.mountpoint, 0);
            
            obj.setElVal(obj.idUI.sUName, state.username, 0);
            obj.setPassword(state.password);
            
            obj.setElVal(obj.idUI.nVLat, state.approx_lat, 0);
            obj.setElVal(obj.idUI.nVLon, state.approx_lon, 0);
            obj.setElVal(obj.idUI.nVH, state.approx_h, 1);

            % Check all the dependencies
            obj.syncFromGUI(obj.idUI.lProcMode);
        end
        
        % Save the stati of the gui to a matlab file
        function exportStateMatlab(obj,filename)
             
            %   MODE
            % ===============================================================
            state.mode              = obj.getElVal(obj.idUI.lProcMode);
            state.nav_mon           = obj.getElVal(obj.idUI.lCaptMode);
            state.kalman_ls         = obj.getElVal(obj.idUI.lAlgType);
            state.code_dd_sa        = obj.getElVal(obj.idUI.lProcType);

            %   INPUT FILE TYPE
            % ===============================================================

            state.rinex_files       = obj.getElVal(obj.idUI.rRin);
            state.gogps_data        = obj.getElVal(obj.idUI.rBin);

            %   OPTIONS
            % ===============================================================

            state.master_pos        = obj.getElVal(obj.idUI.cMPos);
            state.constraint        = obj.getElVal(obj.idUI.cConstraint);
            state.plotproc          = obj.getElVal(obj.idUI.cPlotProc);
            state.ref_path          = obj.getElVal(obj.idUI.cRefPath);
            state.no_skyplot_snr    = obj.getElVal(obj.idUI.cSkyPlot);
            state.google_earth      = obj.getElVal(obj.idUI.cGEarth);
            state.err_ellipse       = obj.getElVal(obj.idUI.cErrEllipse);
            state.plot_master       = obj.getElVal(obj.idUI.cPlotMaster);
            state.plot_amb          = obj.getElVal(obj.idUI.cPlotAmb);
            state.use_ntrip         = obj.getElVal(obj.idUI.cUseNTRIP);
            state.flag_doppler      = obj.getElVal(obj.idUI.cDoppler);
            state.use_sbas          = obj.getElVal(obj.idUI.cUse_SBAS);
            
            %   INTEGER AMBIGUITY RESOLUTION
            % ===============================================================
            
            state.use_lambda        = obj.getElVal(obj.idUI.cLAMBDA);
            state.lambda_method     = obj.getElVal(obj.idUI.lLAMBDAMethod);
            state.lambda_P0         = obj.getElVal(obj.idUI.nP0);
            state.lambda_default_P0 = obj.getElVal(obj.idUI.cP0);
            state.lambda_mu         = obj.getElVal(obj.idUI.nMu);
            state.lambda_auto_mu    = obj.getElVal(obj.idUI.cMu);

            %   INPUT/OUTPUT FILE AND FOLDERS
            % ===============================================================

            state.INIsettings       = obj.getElVal(obj.idUI.sINI);
            state.gogps_data_output = obj.getElVal(obj.idUI.sDirGoOut);
            state.gogps_data_output_prefix = obj.getElVal(obj.idUI.sPrefixGoOut);
            state.activeGNSS        = [obj.isActive(obj.idUI.cGPS) ...
                                       obj.isActive(obj.idUI.cGLONASS) ...
                                       obj.isActive(obj.idUI.cGalileo) ...
                                       obj.isActive(obj.idUI.cBeiDou) ...
                                       obj.isActive(obj.idUI.cQZSS) ...
                                       obj.isActive(obj.idUI.cSBAS) ];

            %   SETTINGS - MASTER STATION
            % ===============================================================

            state.crs               = obj.getElVal(obj.idUI.lCRS);
            state.master_X          = obj.getElVal(obj.idUI.nMX);
            state.master_Y          = obj.getElVal(obj.idUI.nMY);
            state.master_Z          = obj.getElVal(obj.idUI.nMZ);
            state.master_lat        = obj.getElVal(obj.idUI.nMLat);
            state.master_lon        = obj.getElVal(obj.idUI.nMLon);
            state.master_h          = obj.getElVal(obj.idUI.nMh);
            
            %   SETTINGS - KALMAN FILTER - STD
            % ===============================================================
                        
            state.std_X             = obj.getElVal(obj.idUI.nStdE);
            state.std_Y             = obj.getElVal(obj.idUI.nStdN);
            state.std_Z             = obj.getElVal(obj.idUI.nStdU);
            state.std_code          = obj.getElVal(obj.idUI.nStdCode);
            state.std_phase         = obj.getElVal(obj.idUI.nStdPhase);
            state.std_dtm           = obj.getElVal(obj.idUI.nStdDTM);
            state.toggle_std_phase  = obj.getElVal(obj.idUI.bStdPhase);
            state.toggle_std_dtm    = obj.getElVal(obj.idUI.bStdDTM);
            state.std_init          = obj.getElVal(obj.idUI.nStdT0);
            state.std_vel           = obj.getElVal(obj.idUI.nStdVel);
            
            %   SETTINGS - KALMAN FILTER - WEIGHT MODEL
            % ===============================================================
                        
            state.weight_0          = obj.getElVal(obj.idUI.rW0);
            state.weight_1          = obj.getElVal(obj.idUI.rW1);
            state.weight_2          = obj.getElVal(obj.idUI.rW2);
            state.weight_3          = obj.getElVal(obj.idUI.rW3);
            state.weight_4          = obj.getElVal(obj.idUI.rW4);

            %   SETTINGS - KALMAN FILTER
            % ===============================================================

            state.cs_thresh         = obj.getElVal(obj.idUI.nCS);
            state.cut_off           = obj.getElVal(obj.idUI.nCutOff);
            state.snr_thres         = obj.getElVal(obj.idUI.nSNR);
            state.antenna_h         = obj.getElVal(obj.idUI.nHAntenna);
            state.min_sat           = obj.getElVal(obj.idUI.nMinNSat);
            state.stopGOstop        = obj.getElVal(obj.idUI.cStopGoStop);
            
            %   SETTINGS - KALMAN FILTER - DYNAMIC MODEL
            % ===============================================================
                        
            state.dyn_mod           = obj.getElVal(obj.idUI.lDynModel);
            
            %   SETTINGS - KALMAN FILTER - ARAA
            % ===============================================================
            
            state.amb_select        = obj.getElVal(obj.idUI.lARAA);
                
            %   SETTINGS - PORTS
            % ===============================================================
            
            state.captureRate       = get(obj.goh.pumCaptureRate,'Value');
            state.num_receivers     = obj.getElVal(obj.idUI.lnPorts);
            contents = cellstr(get(obj.goh.com_select_0,'String'));
            state.com_select_0 = contents{get(obj.goh.com_select_0,'Value')};
            contents = cellstr(get(obj.goh.com_select_1,'String'));
            state.com_select_1 = contents{get(obj.goh.com_select_1,'Value')};
            contents = cellstr(get(obj.goh.com_select_2,'String'));
            state.com_select_2 = contents{get(obj.goh.com_select_2,'Value')};
            contents = cellstr(get(obj.goh.com_select_3,'String'));
            state.com_select_3 = contents{get(obj.goh.com_select_3,'Value')};                       
            state.protocol_select_0 = obj.getElVal(obj.idUI.lProt0);
            state.protocol_select_1 = obj.getElVal(obj.idUI.lProt1);
            state.protocol_select_2 = obj.getElVal(obj.idUI.lProt2);
            state.protocol_select_3 = obj.getElVal(obj.idUI.lProt3);
            
            %   SETTINGS - MASTER SERVER
            % ===============================================================
            
            state.IP_address        = obj.getElVal(obj.idUI.sIPaddr);
            state.port              = obj.getElVal(obj.idUI.sIPport);
            state.mountpoint        = obj.getElVal(obj.idUI.sMnt);
            state.username          = obj.getElVal(obj.idUI.sUName);
            state.password          = obj.getPassword();
            state.approx_lat        = obj.getElVal(obj.idUI.nVLat);
            state.approx_lon        = obj.getElVal(obj.idUI.nVLon);
            state.approx_h          = obj.getElVal(obj.idUI.nVH);
            
            save(filename, 'state');
        end
    end
      
    %   GO FUNCTIONS (OUTPUT)
    % -------------------------------------------------------------------------
    % This part still needs to be modified (cleaned)
    methods
        % Function to return values to goGPS.m
        function go(obj)
            global goObj goIni
            obj.saveConstellations();

            %master station coordinates
            crs_contents = cellstr(get(obj.goh.crs,'String'));
            master_X = str2double(get(obj.goh.master_X,'String'));
            master_Y = str2double(get(obj.goh.master_Y,'String'));
            master_Z = str2double(get(obj.goh.master_Z,'String'));
            master_lat = str2double(get(obj.goh.master_lat,'String'));
            master_lon = str2double(get(obj.goh.master_lon,'String'));
            master_h = str2double(get(obj.goh.master_h,'String'));
            %KF parameters
            std_init = str2double(get(obj.goh.std_init,'String'));
            std_X = str2double(get(obj.goh.std_X,'String'));
            std_Y = str2double(get(obj.goh.std_Y,'String'));
            std_Z = str2double(get(obj.goh.std_Z,'String'));
            std_vel = str2double(get(obj.goh.std_vel,'String'));
            std_code = str2double(get(obj.goh.std_code,'String'));
            if (get(obj.goh.toggle_std_phase,'Value'))
                std_phase = str2double(get(obj.goh.std_phase,'String'));
            else
                std_phase = 1e30;
            end
            if (get(obj.goh.toggle_std_dtm,'Value'))
                std_dtm = str2double(get(obj.goh.std_dtm,'String'));
            else
                std_dtm = 1e30;
            end
            min_nsat = str2double(get(obj.goh.min_sat,'String'));
            cutoff = str2double(get(obj.goh.cut_off,'String'));
            snr_threshold = str2double(get(obj.goh.snr_thres,'String'));
            cs_threshold = str2double(get(obj.goh.cs_thresh,'String'));
            antenna_h = str2double(get(obj.goh.antenna_h,'String'));
            contents_dyn_mod = cellstr(get(obj.goh.dyn_mod,'String'));
            flag_stopGOstop = get(obj.goh.stopGOstop,'Value');
            %input files
            filerootOUT = [get(obj.goh.sDirGoOut,'String') '/' get(obj.goh.sPrefixGoOut,'String')];
            if obj.isPostProc()
                data_path = goIni.getData('Bin','data_path');
                file_prefix = goIni.getData('Bin','file_prefix');
                filerootIN = [data_path file_prefix];
                data_path = goIni.getData('Receivers','data_path');
                file_name = goIni.getData('Receivers','file_name');
                filename_R_obs = [data_path file_name];
                data_path = goIni.getData('Master','data_path');
                file_name = goIni.getData('Master','file_name');
                filename_M_obs = [data_path file_name];
                data_path = goIni.getData('Navigational','data_path');
                file_name = goIni.getData('Navigational','file_name');
                filename_nav = [data_path file_name];
                ref_path = get(obj.goh.ref_path, 'Value');
                data_path = goIni.getData('RefPath','data_path');
                file_name = goIni.getData('RefPath','file_name');
                filename_ref = [data_path file_name];
                data_path = goIni.getData('DTM','data_path');
                dtm_dir = data_path;
                data_path = goIni.getData('PCO_PCV_file','data_path');
                file_name = goIni.getData('PCO_PCV_file','file_name');
                filename_pco = [data_path file_name];
            end
            
            %serial communication
            % global COMportR
            contents = cellstr(get(obj.goh.com_select_0,'String'));
            COMportR0 = contents{get(obj.goh.com_select_0,'Value')};
            % contents = cellstr(get(obj.goh.com_select_1,'String'));
            % COMportR1 = contents{get(obj.goh.com_select_1,'Value')};
            % contents = cellstr(get(obj.goh.com_select_2,'String'));
            % COMportR2 = contents{get(obj.goh.com_select_2,'Value')};
            % contents = cellstr(get(obj.goh.com_select_3,'String'));
            % COMportR3 = contents{get(obj.goh.com_select_3,'Value')};
            %TCPIP / NTRIP
            flag_NTRIP = get(obj.goh.use_ntrip,'Value');
            master_ip = get(obj.goh.IP_address,'String');
            master_port = str2double(get(obj.goh.port,'String'));
            ntrip_mountpoint = get(obj.goh.mountpoint,'String');
            %functioning mode
            mode = obj.getgoGPSMode();
            
            % If I'm in a mode that uses objects instead of regular code, set goObj flag to 1
            if (mode == goGNSS.MODE_PP_KF_CP_DD_MR)
                goObj = true;
                % init the objects:
            else
                goObj = false;	% set to 1 when goGPS objects are used instead of the regular code
            end
            
            ready = 1;
                        
            if obj.isPostProc()
                %check if the dataset was surveyed with a variable dynamic model
                d = dir([filerootIN '_dyn_000.bin']);
                if (obj.isPostProc && (flag_stopGOstop || strcmp(contents_dyn_mod{get(obj.goh.dyn_mod,'Value')},'Variable')) && isempty(d))
                    msgbox('The selected dataset was not surveyed with a variable dynamic model: please select another dynamic model.'); ready = 0;
                end
            end
            
            if (mode == goGNSS.MODE_RT_R_MON || mode == goGNSS.MODE_RT_RM_MON || mode == goGNSS.MODE_RT_NAV) %if a COM connection to the rover is required
                if(strcmp(COMportR0, 'NA'))
                    msgbox('Please select an existing COM port.'); ready = 0;
                end
            end
            
            if (mode == goGNSS.MODE_RT_M_MON || mode == goGNSS.MODE_RT_RM_MON || mode == goGNSS.MODE_RT_NAV) %if a TCP/IP connection to the master is required
                if (isempty(master_ip))
                    msgbox('Please provide an IP address for the connection to the master.'); ready = 0;
                elseif (isnan(master_port) || master_port < 0 || master_port > 65535)
                    msgbox('Please provide a valid port number for the connection to the master (between 0 and 65535).'); ready = 0;
                end
                if (flag_NTRIP) %if a NTRIP connection is required
                    if (isempty(ntrip_mountpoint))
                        msgbox('Please provide a mountpoint for the NTRIP connection.'); ready = 0;
                    end
                end
            end
            
            if (obj.isEnabled(obj.idUI.cLAMBDA) && obj.isActive(obj.idUI.cLAMBDA) && ~obj.isLambda2) %if LAMBDA3.x is requested
                if (~exist('LAMBDA.m','file') || ~exist('ratiotab.mat','file'))
                    msgbox(['LAMBDA 3.x code not found. Please download it from http://gnss.curtin.edu.au/research/lambda.cfm and place it in the working path (e.g. in ./positioning/lambda/lambda_v3). Switching ambiguity resolution method to' obj.strLAMBDAMethod{obj.idILS_enum_old} '...']); ready = 0;
                    obj.setElVal(obj.idUI.lLAMBDAMethod,1);
                end
            end
            
            if (ready)
                % If the working folder does not exist
                if isempty(dir(obj.getSettingsDir()))
                    msgbox('Non existent settings folder. It has beeen probably erased! It is not possible to save the last settings.');
                else
                    obj.exportStateMatlab([obj.getSettingsDir() obj.lastSettingsFile]);
                end
                uiresume(obj.goh.main_panel);
            end
        end
        
        % Function to return values to goGPS.m
        function funout = outputFun(obj)
            global goIni;
            if isempty(goIni)
                goIni = goIniReader;
            end
            
            obj.saveConstellations();
            mode = obj.getgoGPSMode();
            mode_vinc = get(obj.goh.constraint,'Value') * obj.isActive(obj.idUI.cConstraint);
            if (get(obj.goh.file_type, 'SelectedObject') == obj.goh.rinex_files)
                mode_data = 0;
            else %goGPS data
                mode_data = 1;
            end
            contents_dyn_mod = cellstr(get(obj.goh.dyn_mod,'String'));
            if (strcmp(contents_dyn_mod{get(obj.goh.dyn_mod,'Value')},'Variable') || get(obj.goh.stopGOstop,'Value'))
                flag_var_dyn_model = 1;
            else
                flag_var_dyn_model = 0;
            end

            % atm model
            iono = goIni.getData('ATM_model','iono');
            tropo = goIni.getData('ATM_model','tropo');
            if (isempty(iono))
                iono_model = 1;
            else
                iono_model = iono;
            end
            if (isempty(tropo))
                tropo_model = 1;
            else
                tropo_model = tropo;
            end
            % mixed
            fsep = goIni.getData('Various','field_separator');
            if (isempty(fsep))
                fsep_char = 'default';
            else
                fsep_char = fsep;
            end           
            
            mode_ref = get(obj.goh.ref_path,'Value');
            flag_ms_pos = get(obj.goh.master_pos,'Value');
            flag_ms = get(obj.goh.plot_master,'Value');
            flag_ge = get(obj.goh.google_earth,'Value');
            flag_cov = get(obj.goh.err_ellipse,'Value');
            flag_NTRIP = get(obj.goh.use_ntrip,'Value');
            flag_amb = get(obj.goh.plot_amb,'Value');
            flag_skyplot = get(obj.goh.no_skyplot_snr,'Value');
            flag_plotproc = get(obj.goh.plotproc,'Value');
            flag_stopGOstop = get(obj.goh.stopGOstop,'Value');
            flag_SBAS = get(obj.goh.use_SBAS,'Value');
            flag_IAR = get(obj.goh.cLAMBDA,'Value');
            filerootOUT = [get(obj.goh.sDirGoOut,'String') '/' get(obj.goh.sPrefixGoOut,'String')];
            if (obj.isPostProc) % I need these informations only in Post Processing
                data_path = goIni.getData('Bin','data_path');
                file_prefix = goIni.getData('Bin','file_prefix');
                filerootIN = [data_path file_prefix];
                i = 1;
                j = length(filerootOUT);
                while (~isempty(dir([filerootOUT '_????_rover.bin'])) || ...
                        ~isempty(dir([filerootOUT '_master*.bin'])) || ...
                        ~isempty(dir([filerootOUT '_????_obs*.bin'])) || ...
                        ~isempty(dir([filerootOUT '_????_eph*.bin'])) || ...
                        ~isempty(dir([filerootOUT '_????_dyn*.bin'])) || ...
                        ~isempty(dir([filerootOUT '_sat*.bin'])) || ...
                        ~isempty(dir([filerootOUT '_kal*.bin'])) || ...
                        ~isempty(dir([filerootOUT '_dt*.bin'])) || ...
                        ~isempty(dir([filerootOUT '_conf*.bin'])) || ...
                        ~isempty(dir([filerootOUT '_dop*.bin'])) || ...
                        ~isempty(dir([filerootOUT '_ECEF*.txt'])) || ...
                        ~isempty(dir([filerootOUT '_geod*.txt'])) || ...
                        ~isempty(dir([filerootOUT '_plan*.txt'])) || ...
                        ~isempty(dir([filerootOUT '_????_NMEA*.txt'])) || ...
                        ~isempty(dir([filerootOUT '.kml'])) )
                    
                    filerootOUT(j+1:j+4) = ['_' num2str(i,'%03d')];
                    i = i + 1;
                end
                data_path = goIni.getData('Receivers','data_path');
                file_name = goIni.getData('Receivers','file_name');
                filename_R_obs = [data_path file_name];
                data_path = goIni.getData('Master','data_path');
                file_name = goIni.getData('Master','file_name');
                filename_M_obs = [data_path file_name];
                data_path = goIni.getData('Navigational','data_path');
                file_name = goIni.getData('Navigational','file_name');
                filename_nav = [data_path file_name];
                flag_SP3 = goIni.getData('Navigational','isSP3');
                if isempty(flag_SP3)
                    if (strcmpi(filename_nav(end-3:end),'.sp3'))
                        flag_SP3 = 1;
                    else
                        flag_SP3 = 0;
                    end
                end
                data_path = goIni.getData('RefPath','data_path');
                file_name = goIni.getData('RefPath','file_name');
                filename_ref = [data_path file_name];
                data_path = goIni.getData('PCO_PCV_file','data_path');
                file_name = goIni.getData('PCO_PCV_file','file_name');
                filename_pco = [data_path file_name];
                if(obj.isMultiReceiver)
                    [multi_antenna_rf, ~] = goIni.getGeometry();
                else
                    multi_antenna_rf = [];
                end
            else
                filerootIN = '';
                filename_R_obs = '';
                filename_M_obs = '';
                filename_nav = '';
                filename_ref = '';
                filename_pco = '';
                flag_SP3 = 0;
                rates = get(obj.goh.pumCaptureRate,'String');                
                goIni.setCaptureRate(rates{get(obj.goh.pumCaptureRate,'Value')});
                multi_antenna_rf = [];
            end
            
            contents = cellstr(get(obj.goh.crs,'String'));
            if (strcmp(contents{get(obj.goh.crs,'Value')},'ECEF (X,Y,Z)'))
                XM = str2double(get(obj.goh.master_X,'String'));
                YM = str2double(get(obj.goh.master_Y,'String'));
                ZM = str2double(get(obj.goh.master_Z,'String'));
            else
                latM = str2double(get(obj.goh.master_lat,'String'));
                lonM = str2double(get(obj.goh.master_lon,'String'));
                hM = str2double(get(obj.goh.master_h,'String'));
                [XM, YM, ZM] = geod2cart (latM*pi/180, lonM*pi/180, hM, 6378137, 1/298.257222101);
            end
            pos_M_man = [XM; YM; ZM];
            
            contents = cellstr(get(obj.goh.num_receivers,'String'));
            num_rec = str2double(contents{get(obj.goh.num_receivers,'Value')});
            protocol_idx = nan(4,1);
            
            if num_rec >= 1
                contentsProt = cellstr(get(obj.goh.protocol_select_0,'String'));
                if (strcmp(contentsProt{get(obj.goh.protocol_select_0,'Value')},'UBX (u-blox)'))
                    protocol_idx(1) = 0;
                elseif (strcmp(contentsProt{get(obj.goh.protocol_select_0,'Value')},'iTalk (Fastrax)'))
                    protocol_idx(1) = 1;
                elseif (strcmp(contentsProt{get(obj.goh.protocol_select_0,'Value')},'SkyTraq'))
                    protocol_idx(1) = 2;
                elseif (strcmp(contentsProt{get(obj.goh.protocol_select_0,'Value')},'BINR (NVS)'))
                    protocol_idx(1) = 3;
                end
                
                if num_rec >= 2
                    contentsProt = cellstr(get(obj.goh.protocol_select_1,'String'));
                    if (strcmp(contentsProt{get(obj.goh.protocol_select_1,'Value')},'UBX (u-blox)'))
                        protocol_idx(2) = 0;
                    elseif (strcmp(contentsProt{get(obj.goh.protocol_select_1,'Value')},'iTalk (Fastrax)'))
                        protocol_idx(2) = 1;
                    elseif (strcmp(contentsProt{get(obj.goh.protocol_select_1,'Value')},'SkyTraq'))
                        protocol_idx(2) = 2;
                    elseif (strcmp(contentsProt{get(obj.goh.protocol_select_1,'Value')},'BINR (NVS)'))
                        protocol_idx(2) = 3;
                    end
                    
                    if num_rec >= 3
                        contentsProt = cellstr(get(obj.goh.protocol_select_2,'String'));
                        if (strcmp(contentsProt{get(obj.goh.protocol_select_2,'Value')},'UBX (u-blox)'))
                            protocol_idx(3) = 0;
                        elseif (strcmp(contentsProt{get(obj.goh.protocol_select_2,'Value')},'iTalk (Fastrax)'))
                            protocol_idx(3) = 1;
                        elseif (strcmp(contentsProt{get(obj.goh.protocol_select_2,'Value')},'SkyTraq'))
                            protocol_idx(3) = 2;
                        elseif (strcmp(contentsProt{get(obj.goh.protocol_select_2,'Value')},'BINR (NVS)'))
                            protocol_idx(3) = 3;
                        end
                        
                        if num_rec >= 4
                            contentsProt = cellstr(get(obj.goh.protocol_select_3,'String'));
                            if (strcmp(contentsProt{get(obj.goh.protocol_select_3,'Value')},'UBX (u-blox)'))
                                protocol_idx(4) = 0;
                            elseif (strcmp(contentsProt{get(obj.goh.protocol_select_3,'Value')},'iTalk (Fastrax)'))
                                protocol_idx(4) = 1;
                            elseif (strcmp(contentsProt{get(obj.goh.protocol_select_3,'Value')},'SkyTraq'))
                                protocol_idx(4) = 2;
                            elseif (strcmp(contentsProt{get(obj.goh.protocol_select_3,'Value')},'BINR (NVS)'))
                                protocol_idx(4) = 3;
                            end
                        end
                    end
                end
            end
            protocol_idx = protocol_idx(~isnan(protocol_idx));
            
            funout = cell(30,1);
            
            funout{1} = mode;
            funout{2} = mode_vinc;
            funout{3} = mode_data;
            funout{4} = mode_ref;
            funout{5} = flag_ms_pos;
            funout{6} = flag_ms;
            funout{7} = flag_ge;
            funout{8} = flag_cov;
            funout{9} = flag_NTRIP;
            funout{10} = flag_amb;
            funout{11} = flag_skyplot;
            funout{12} = flag_plotproc;
            funout{13} = flag_var_dyn_model;
            funout{14} = flag_stopGOstop;
            funout{15} = flag_SP3;
            funout{16} = flag_SBAS;
            funout{17} = flag_IAR;
            funout{18} = filerootIN;
            funout{19} = filerootOUT;
            funout{20} = filename_R_obs;
            funout{21} = filename_M_obs;
            funout{22} = filename_nav;
            funout{23} = filename_ref;
            funout{24} = filename_pco;
            funout{25} = pos_M_man;
            funout{26} = protocol_idx;
            funout{27} = multi_antenna_rf;
            funout{28} = iono_model;
            funout{29} = tropo_model;            
            funout{30} = fsep_char;
            
            global sigmaq0 sigmaq_vE sigmaq_vN sigmaq_vU sigmaq_vel
            global sigmaq_cod1 sigmaq_cod2 sigmaq_ph sigmaq0_N sigmaq_dtm
            global min_nsat cutoff snr_threshold cs_threshold weights snr_a snr_0 snr_1 snr_A order o1 o2 o3
            global h_antenna
            global tile_header tile_georef dtm_dir
            global master_ip master_port ntrip_user ntrip_pw ntrip_mountpoint
            global nmea_init
            global flag_doppler_cs
            global COMportR
            global IAR_method P0 mu flag_auto_mu flag_default_P0

            IAR_method = get(obj.goh.lLAMBDAMethod,'Value') - 1;
            P0 = str2double(get(obj.goh.nP0,'String'));
            mu = str2double(get(obj.goh.nMu,'String'));
            flag_auto_mu = get(obj.goh.cMu,'Value') && ~obj.isLambda2;
            flag_default_P0 = get(obj.goh.cP0,'Value') && ~obj.isLambda2;
            
            contents = cellstr(get(obj.goh.com_select_0,'String'));
            COMportR0 = contents{get(obj.goh.com_select_0,'Value')};
            contents = cellstr(get(obj.goh.com_select_1,'String'));
            COMportR1 = contents{get(obj.goh.com_select_1,'Value')};
            contents = cellstr(get(obj.goh.com_select_2,'String'));
            COMportR2 = contents{get(obj.goh.com_select_2,'Value')};
            contents = cellstr(get(obj.goh.com_select_3,'String'));
            COMportR3 = contents{get(obj.goh.com_select_3,'Value')};
            
            if num_rec >= 1
                COMportR{1,1} = COMportR0;
                if num_rec >= 2
                    COMportR{2,1} = COMportR1;
                    if num_rec >= 3
                        COMportR{3,1} = COMportR2;
                        if num_rec >= 4
                            COMportR{4,1} = COMportR3;
                        end
                    end
                end
            end
            
            flag_doppler_cs = get(obj.goh.flag_doppler,'Value');
            sigmaq0 = str2double(get(obj.goh.std_init,'String'))^2;
            sigmaq_vE = str2double(get(obj.goh.std_X,'String'))^2;
            sigmaq_vN = str2double(get(obj.goh.std_Y,'String'))^2;
            sigmaq_vU = str2double(get(obj.goh.std_Z,'String'))^2;
            sigmaq_vel = str2double(get(obj.goh.std_vel,'String'))^2;
            sigmaq_cod1 = str2double(get(obj.goh.std_code,'String'))^2;
            sigmaq_cod2 = 0.16;
            if (get(obj.goh.toggle_std_phase,'Value'))
                sigmaq_ph = str2double(get(obj.goh.std_phase,'String'))^2;
            else
                sigmaq_ph = 1e30;
            end
            sigmaq0_N = 1000;
            if (get(obj.goh.toggle_std_dtm,'Value'))
                sigmaq_dtm = str2double(get(obj.goh.std_dtm,'String'))^2;
            else
                sigmaq_dtm = 1e30;
            end
            min_nsat = str2double(get(obj.goh.min_sat,'String'));
            if (mode == 2)
                disp('Minimum number of satellites is forced to 4 (for undifferenced positioning)');
                min_nsat = 4;
            end
            goIni.addSection('Generic');
            goIni.addKey('Generic','cutoff', str2double(get(obj.goh.cut_off,'String')));
            cutoff = str2double(get(obj.goh.cut_off,'String'));
            goIni.addKey('Generic','snrThr', str2double(get(obj.goh.snr_thres,'String')));
            snr_threshold = str2double(get(obj.goh.snr_thres,'String'));
            goIni.addKey('Generic','csThr', str2double(get(obj.goh.cs_thresh,'String')));
            cs_threshold = str2double(get(obj.goh.cs_thresh,'String'));
            if (get(obj.goh.weight_select, 'SelectedObject') == obj.goh.weight_0)
                weights = 0;
            elseif (get(obj.goh.weight_select, 'SelectedObject') == obj.goh.weight_1)
                weights = 1;
            elseif (get(obj.goh.weight_select, 'SelectedObject') == obj.goh.weight_2)
                weights = 2;
            elseif (get(obj.goh.weight_select, 'SelectedObject') == obj.goh.weight_3)
                weights = 3;
            elseif (get(obj.goh.weight_select, 'SelectedObject') == obj.goh.weight_4)
                weights = 4;
            end
            snr_a = 30;
            snr_0 = 10;
            snr_1 = 50;
            snr_A = 30;
            global amb_restart_method
            contents = cellstr(get(obj.id2handle(obj.idUI.lARAA),'String'));
            selection = contents{min(get(obj.id2handle(obj.idUI.lARAA),'Value'), length(contents))};
            if (strcmp(selection, 'Observed code - phase difference'))
                amb_restart_method = 0;
            elseif (strcmp(selection, 'Kalman-predicted code - phase difference'))
                amb_restart_method = 1;
            else
                amb_restart_method = 2;
            end
            contents = cellstr(get(obj.id2handle(obj.idUI.lDynModel),'String'));
            if (strcmp(contents{min(get(obj.id2handle(obj.idUI.lDynModel),'Value'), length(contents))},'Static'))
                order = 1;
            else
                if (strcmp(contents{min(get(obj.id2handle(obj.idUI.lDynModel),'Value'), length(contents))},'Const. acceleration'))
                    order = 3;
                elseif (strcmp(contents{min(get(obj.id2handle(obj.idUI.lDynModel),'Value'), length(contents))},'Const. velocity'))
                    order = 2;
                else
                    order = 1;
                end
            end
            
            o1 = order;
            o2 = order*2;
            o3 = order*3;
            h_antenna = str2double(get(obj.goh.antenna_h,'String'));
%             if (obj.isPostProc) % I need these informations only in Post Processing
                dtm_dir = goIni.getData('DTM','data_path');
                try
                    load([dtm_dir '/tiles/tile_header'], 'tile_header');
                    load([dtm_dir '/tiles/tile_georef'], 'tile_georef');
                catch e
                    tile_header.nrows = 0;
                    tile_header.ncols = 0;
                    tile_header.cellsize = 0;
                    tile_header.nodata = 0;
                    tile_georef = zeros(1,1,4);
                end
%             end
            master_ip = get(obj.goh.IP_address,'String');
            master_port = str2double(get(obj.goh.port,'String'));
            ntrip_user = get(obj.goh.username,'String');
            ntrip_pw = obj.getPassword();
            ntrip_mountpoint = get(obj.goh.mountpoint,'String');
            phiApp = str2double(get(obj.goh.approx_lat,'String'));
            lamApp = str2double(get(obj.goh.approx_lon,'String'));
            hApp = str2double(get(obj.goh.approx_h,'String'));
            [XApp,YApp,ZApp] = geod2cart (phiApp*pi/180, lamApp*pi/180, hApp, 6378137, 1/298.257222101);
            if ~isnan(XApp) && ~isnan(YApp) && ~isnan(ZApp)
                nmea_init = NMEA_GGA_gen([XApp YApp ZApp],10);
            else
                nmea_init = '';
            end
        end
        
        % Function to save in the goIni object the status of activation of
        % the various GNSS
        function saveConstellations(obj)
            global goIni
            if isempty(goIni)
                goIni = goIniReader;
            end            
            goIni.addSection('Constellations');
            goIni.addKey('Constellations','GPS',obj.isActive(obj.idUI.cGPS));
            goIni.addKey('Constellations','GLONASS',obj.isActive(obj.idUI.cGLONASS));
            goIni.addKey('Constellations','Galileo',obj.isActive(obj.idUI.cGalileo));
            goIni.addKey('Constellations','BeiDou',obj.isActive(obj.idUI.cBeiDou));
            goIni.addKey('Constellations','QZSS',obj.isActive(obj.idUI.cQZSS));
            goIni.addKey('Constellations','SBAS',obj.isActive(obj.idUI.cSBAS));
        end
    end        
    

    
    
    
    
    
    %   GO FUNCTIONS (OUTPUT)
    % -------------------------------------------------------------------------
    % This part still needs to be modified (cleaned)
    methods
        % Function to load the Edit INI window
        function openEditINI(obj)
            if (isfield(obj.edtINI,'h'))
                if ishandle(obj.edtINI.h.wEditINI)
                    close(obj.edtINI.h.wEditINI);
                end
                delete obj.edtINI.h
            end
            if (obj.interfaceOS == obj.isUnix)
                guiEditINI_unix();
            else
                guiEditINI();
            end
        end
        
        % Function to init the INI editor
        % creates, objects, load default valuesm, etc...
        function initEditINI(obj, h)
            global goIni
            
            % Save handler to the INI editor
            obj.edtINI.h = h;
            
            % Get the name of the ini file
            
            if isobject(goIni)
                fileName = goIni.getFileName();
            else
                fileName = '';
            end
            if isempty(fileName)
                fileName = [obj.getSettingsDir() obj.defaultINIFile];
                if ~exist(fileName, 'file');
                    fileName = '';
                end
            end
            
            % Get the content of the ini file
            obj.setGuiElStr(h.sINI, fileName);
            obj.setGuiElStr(h.sINIout, fileName);

            if ~isempty(fileName)
                fid = fopen(fileName,'r');
                text = fread(fid, '*char');
                text = text';
                fclose(fid);
            else
                text = '# No INI file found';
            end
            
            % Undocumented edit box => better management of a text file
            % Create the widget containing the text
            jCodePaneINI = com.mathworks.widgets.SyntaxTextPane;
            jCodePaneINI.setText(text);
            % Create the ScrollPanel containing the widget
            jScrollPaneINI = com.mathworks.mwswing.MJScrollPane(jCodePaneINI);
            % Substitute the eINI edit box with the Java Scroll Pane
            set(obj.edtINI.h.eINI, 'Units', 'pixels');
            [jhPanel, hContainerINI] = javacomponent(jScrollPaneINI,get(obj.edtINI.h.eINI,'Position'),h.wEditINI);
            delete(obj.edtINI.h.eINI);
            
            % Save the new object
            obj.edtINI.jEdit.jINI = jCodePaneINI;
            obj.edtINI.jEdit.hINI = hContainerINI;
            
            % Load INI keywords
            obj.edtINI.keywordsINI = goIniReader([obj.getSettingsDir() obj.defaultINIKeywordsFile], 0);
            obj.edtINI.keywordsINI.readFile();
            % Sections
            sections = obj.edtINI.keywordsINI.getData('INI','sections');
            obj.setGuiElStr(h.lSections, sections);
            drawnow;
            % Fields
            obj.updateFieldsINI();
            
            % Replacing the field
            jCodePaneFields = com.mathworks.widgets.SyntaxTextPane;
            jCodePaneFields.setText('');
            % Create the ScrollPanel containing the widget
            jScrollPaneFields = com.mathworks.mwswing.MJScrollPane(jCodePaneFields);
            % Substitute the eFields edit box with the Java Scroll Pane
            set(obj.edtINI.h.eFields, 'Units', 'pixels');
            [jhPanel, hContainerFields] = javacomponent(jScrollPaneFields,get(obj.edtINI.h.eFields,'Position'),h.wEditINI);
            delete(obj.edtINI.h.eFields);
            
            % Save the new object
            obj.edtINI.jEdit.jFields = jCodePaneFields;
            obj.edtINI.jEdit.hFields = hContainerFields;
            drawnow;
            obj.updateFieldsINI();
            
            % Replacing the field
            jCodePaneBrowse = com.mathworks.widgets.SyntaxTextPane;
            jCodePaneBrowse.setText('');
            % Create the ScrollPanel containing the widget
            jScrollPaneBrowse = com.mathworks.mwswing.MJScrollPane(jCodePaneBrowse);
            % Substitute the eFields edit box with the Java Scroll Pane
            set(obj.edtINI.h.sBrowse, 'Units', 'pixels');
            [jhPanel, hContainerBrowse] = javacomponent(jScrollPaneBrowse,get(obj.edtINI.h.sBrowse,'Position'),h.wEditINI);
            delete(obj.edtINI.h.sBrowse);
            
            % Save the new object
            obj.edtINI.jEdit.jBrowse = jCodePaneBrowse;
            obj.edtINI.jEdit.hBrowse = hContainerBrowse;
            drawnow;
        end
        
        % Browse INI file => select a file to be edited
        function browseINIEditInFile(obj)
            % In multi receiver mode, I read from ini file
            [filename, pathname] = uigetfile( ...
                {'*.ini;','INI configuration file (*.ini)'; ...
                '*.*',  'All Files (*.*)'}, ...
                'Choose an INI configuration file',[obj.getSettingsDir()]);
            if (filename ~= 0)
                fileName = [pathname filename];
                obj.setGuiElStr(obj.edtINI.h.sINI, fileName);
                obj.setGuiElStr(obj.edtINI.h.sINIout, fileName);
                
                % Get the name of the ini file
                if isempty(fileName)
                    fileName = [obj.getSettingsDir() obj.defaultINIFile];
                    if ~exist(fileName, 'file');
                        fileName = '';
                    end
                end
                
                % Get the content of the ini file
                if ~isempty(fileName)
                    fid = fopen(fileName,'r');
                    text = fread(fid, '*char');
                    text = text';
                    fclose(fid);
                else
                    text = '# No INI file found';
                end
                obj.edtINI.jEdit.jINI.setText(text);
                drawnow;                
            end
        end
     
        % Browse 4 INI
        function saveINI(obj)
            filename = get(obj.edtINI.h.sINIout, 'String');
            try
                fid = fopen(filename,'w');
                fwrite(fid, char(obj.edtINI.jEdit.jINI.getText()));
                fclose(fid);
                
                msgbox('The file has been saved correctly');
                % If the main goGPS interface exist
                if (ishandle(obj.goh.main_panel))
                    if (filename ~= 0)
                        obj.setElVal(obj.idUI.sINI, fullfile(filename));
                    end
                    obj.forceINIupdate();
                end
            catch e
                msgbox(['Error: ' e.message ' Please provide a valid filename(path)']);
            end
        end
        
        % Update the fields according to the selected section
        function updateFieldsINI(obj)
            curSection = get(obj.edtINI.h.lSections, 'String');
            valSection = get(obj.edtINI.h.lSections, 'Value');
            if length(curSection) > 1
                curSection = curSection{valSection};
                
                numFields = obj.edtINI.keywordsINI.getData(curSection, 'num_fields');
                strFields = obj.edtINI.keywordsINI.getData(curSection, 'str_fields');
                vectFields = obj.edtINI.keywordsINI.getData(curSection, 'vect_fields');
                
                str = sprintf('[%s]\n', curSection);
                if ~isempty(numFields)
                    if iscell(numFields)
                        for i=1:length(numFields)
                            str = sprintf('%s%s = 0\n',str, numFields{i});
                        end
                    else
                        str = sprintf('%s%s = 0\n',str, numFields);
                    end
                end
                if ~isempty(strFields)
                    if iscell(strFields)
                        for i=1:length(strFields)
                            str = sprintf('%s%s = "insert a string"\n',str, strFields{i});
                        end
                    else
                        str = sprintf('%s%s = "insert a string"\n',str, strFields);
                    end
                end
                if ~isempty(vectFields)
                    if iscell(vectFields)
                        for i=1:length(vectFields)
                            str = sprintf('%s%s = [0 0 0]\n',str, vectFields{i});
                        end
                    else
                        str = sprintf('%s%s = [0 0 0]\n',str, vectFields);
                    end
                end
                if isfield(obj.edtINI.jEdit,'jFields')
                    obj.edtINI.jEdit.jFields.setText(str);
                end
            end
        end
        
        % Browse for RINEX file
        function browse4Rin(obj)
            % In multi receiver mode, I read from ini file
            [filename, pathname] = uigetfile( ...
                    {'*.obs;*.??o;*.??O','RINEX observation files (*.obs,*.??o,*.??O)';
                    '*.obs','Observation files (*.obs)'; ...
                    '*.??o;*.??O','Observation files (*.??o,*.??O)'; ...
                    '*.*',  'All Files (*.*)'}, ...
                    'MultiSelect', 'on', ...
                    'Choose a RINEX observation file',[obj.getWorkingDir() 'data_RINEX']);
                
            if ~isempty(filename)
                obj.workingDir = [pathname];
                str = sprintf('data_path = "%s"\n', pathname);
                if iscell(filename)
                    str = sprintf('nRec = %d\n%sfile_name = [', length(filename), str);
                    for r=1:length(filename)
                        str = sprintf('%s "%s"',str, filename{r});                      
                    end
                    str = sprintf('%s ]',str);
                else
                    str = sprintf('%sfile_name = "%s"', str, filename);
                end
                obj.edtINI.jEdit.jBrowse.setText(str);
                clipboard('copy', str);
            end
        end
        
        % Browse for a navigation file
        function browse4Nav(obj)
            [filename, pathname] = uigetfile( ...
                {'*.nav;*.??n;*.??N','RINEX navigation files (*.nav,*.??n,*.??N)';
                '*.nav','Navigation files (*.nav)'; ...
                '*.??n;*.??N','Navigation files (*.??n,*.??N)'; ...
                '*.*',  'All Files (*.*)'}, ...
                'Choose a RINEX navigation file',[obj.getWorkingDir() 'data_RINEX']);
            
            if (filename ~= 0)
                obj.workingDir = [pathname];
                str = sprintf('data_path = "%s"\nfile_name = "%s"', pathname, filename);
                obj.edtINI.jEdit.jBrowse.setText(str);
                clipboard('copy', str);
            end
        end
        
        % Browse for a binary file
        function browse4Bin(obj)
            [filename, pathname] = uigetfile( ...
                {'*.bin','goGPS binary data (*.bin)'}, ...
                'Choose goGPS binary data',obj.getWorkingDir());
            
            if (filename ~= 0)
                pos = find(filename == '_');
                filename = filename(1:pos(end-1)-1);
                str = sprintf('data_path = "%s"\nfile_prefix = "%s"', pathname, filename);
                obj.edtINI.jEdit.jBrowse.setText(str);
                clipboard('copy', str);
            end
        end
        
        % Browse output folder
        function browse4Dir(obj)
            dname = uigetdir(obj.getWorkingDir(),'Choose a directory');
            if (dname ~= 0)
                obj.workingDir = [dname '../'];
                str = sprintf('data_path = "%s"', dname);
                obj.edtINI.jEdit.jBrowse.setText(str);
                clipboard('copy', str);
            end
        end
        
        % Browse for the path containing reference points for constrained solutions
        function browse4Ref(obj)
            [filename, pathname] = uigetfile('*.mat', 'Choose file containing reference path',obj.getWorkingDir());
            
            if (filename ~= 0)
                obj.workingDir = [pathname];
                str = sprintf('data_path = "%s"\nfile_name = "%s"', pathname, filename);
                obj.edtINI.jEdit.jBrowse.setText(str);
                clipboard('copy', str);
            end
            obj.updateGUI();
        end
        
        % Browse for a Generic File
        function browse4Gen(obj)
            [filename, pathname] = uigetfile( ...
                {'*.*',  'All Files (*.*)'}, ...
                'Choose a file',obj.getWorkingDir());
            
            if (filename ~= 0)
                str = sprintf('data_path = "%s"\nfile_name = "%s"', pathname, filename);
                obj.edtINI.jEdit.jBrowse.setText(str);
                clipboard('copy', str);
            end
        end
        
    end  
    
    %   GUI STATIC MODIFIERS
    % -------------------------------------------------------------------------
    % Functions to set some properties of an GUI object
    methods(Static, Access = 'private')
        
        % Enable/Disable a generic element of the interface
        function onoffGuiEl(hObject, state)
            if nargin < 2
                state = 'on';
            end
            set(hObject, 'Enable', state);
        end

        % Enable/Disable a panel element of the interface
        function onoffGuiPanel(hObject, state)
            if nargin < 2
                state = 'on';
            end
            if strcmp(state,'off')
                set(hObject, 'ForegroundColor', goGUIclass.disableCol);
            else
                set(hObject, 'ForegroundColor', goGUIclass.enableCol);
            end
        end        
     
        % Get enabled status of a generic element
        function state = isGuiElOn(hObject)
            state = strcmp(get(hObject, 'Enable'),'on');
        end
        
        % Get enabled status of a panel element
        function state = isGuiPanelOn(hObject)
            state = isequal(get(hObject, 'ForegroundColor'), goGUIclass.enableCol);
        end
        
        % Get a value from an element of the interface
        function val = getGuiElVal(hObject)
            val = get(hObject, 'Value');
        end
        
        % Get a value from an element of the interface
        function val = getGuiElColor(hObject)
            val = get(hObject, 'ForegroundColor');
        end
        
        % Get a string from an element of the interface
        function str = getGuiElStr(hObject)
            str = get(hObject, 'String');
        end

        % Get a title from an element of the interface
        function str = getGuiElTitle(hObject)
            str = get(hObject, 'Title');
        end

        % Set a value of an element of the interface
        function setGuiElVal(hObject, value)
            set(hObject, 'Value', value);
        end

        % Set a value of an element of the interface
        function setGuiElColor(hObject, color)
            set(hObject, 'ForegroundColor', color);
        end
        
        % Set a string of an element of the interface
        function setGuiElStr(hObject, str)
            set(hObject, 'String', str);
        end
        
        % Get a title from an element of the interface
        function setGuiElTitle(hObject, str)
            set(hObject, 'Title', str);
        end
        
        % Close box
        function closeGUI(src,evnt)
            global goGUI % I cannot pass the object directly
            selection = questdlg('Do you want to quit goGPS?',...
                '',...
                'Yes','No','Yes');
            switch selection,
                case 'Yes',
                    if isempty(goGUI)  % Something went wrong in the initialization phase
                        delete(gcf);
                    else
                        if isfield(goGUI.goh, 'main_panel');
                            delete(goGUI.goh.main_panel);
                        else % Something went wrong in the initialization phase
                            delete(gcf);
                        end
                    end
                case 'No'
                    return
            end
        end

    end    
        
end
