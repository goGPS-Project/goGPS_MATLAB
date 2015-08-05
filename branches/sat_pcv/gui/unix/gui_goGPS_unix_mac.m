function varargout = gui_goGPS_unix_mac(varargin)
% GUI_GOGPS_UNIX M-file for gui_goGPS_unix_mac.fig
%      GUI_GOGPS_UNIX, by itself, creates a new GUI_GOGPS_UNIX or raises the existing
%      singleton*.
%
%      H = GUI_GOGPS_UNIX returns the handle to a new GUI_GOGPS_UNIX or the handle to
%      the existing singleton*.
%
%      GUI_GOGPS_UNIX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_GOGPS_UNIX.M with the given input arguments.
%
%      GUI_GOGPS_UNIX('Property','Value',...) creates a new GUI_GOGPS_UNIX or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_goGPS_unix_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_goGPS_unix_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_goGPS_unix_mac

% Last Modified by GUIDE v2.5 19-Apr-2013 12:12:09

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Ivan Reguzzoni
%----------------------------------------------------------------------------------------------
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
%----------------------------------------------------------------------------------------------

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name', mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gui_goGPS_unix_OpeningFcn, ...
    'gui_OutputFcn',  @gui_goGPS_unix_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before gui_goGPS_unix_mac is made visible.
function gui_goGPS_unix_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_goGPS_unix_mac (see VARARGIN)
clearvars -global goGUI
global goGUI
goGUI = goGUIclass(handles, goGUIclass.isUnix);

% UIWAIT makes gui_goGPS_unix_mac wait for user response (see UIRESUME)
uiwait(hObject);


% --- Outputs from this function are returned to the command line.
function varargout = gui_goGPS_unix_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global goGUI
    % If I press the exit button
    if(~isstruct(handles))
        varargout = cell(30,1);
        return
    end
    
    varargout = goGUI.outputFun();
%close main panel
delete(gcf)

% --------------------------------------------------------------------

%   MODE
% ===============================================================

% --- Executes on selection change in mode.
function mode_Callback(hObject, eventdata, handles)
% hObject    handle to mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns mode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mode
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.lProcMode);

% --- Executes during object creation, after setting all properties.
function mode_CreateFcn(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in nav_mon.
function nav_mon_Callback(hObject, eventdata, handles)
% hObject    handle to nav_mon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns nav_mon contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nav_mon
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.lCaptMode);

% --- Executes during object creation, after setting all properties.
function nav_mon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nav_mon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in kalman_ls.
function kalman_ls_Callback(hObject, eventdata, handles)
% hObject    handle to kalman_ls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns kalman_ls contents as cell array
%        contents{get(hObject,'Value')} returns selected item from kalman_ls
%enable Kalman filters settings
global goGUI;
    goGUI.syncFromGUI(goGUI.idUI.lAlgType);
    
% --- Executes during object creation, after setting all properties.
function kalman_ls_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kalman_ls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in code_dd_sa.
function code_dd_sa_Callback(hObject, eventdata, handles)
% hObject    handle to code_dd_sa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns code_dd_sa contents as cell array
%        contents{get(hObject,'Value')} returns selected item from code_dd_sa
global goGUI;
    goGUI.syncFromGUI(goGUI.idUI.lProcType);
    
% --- Executes during object creation, after setting all properties.
function code_dd_sa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to code_dd_sa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%   INPUT FILE TYPE
% ===============================================================
    
% --- Executes when selected object is changed in file_type.
function file_type_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in file_type
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global goGUI;
    goGUI.syncFromGUI(goGUI.idGroup.gIFiles);
    
%   OPTIONS
% ===============================================================
    
% --- Executes on button press in constraint.
function constraint_Callback(hObject, eventdata, handles)
% hObject    handle to constraint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of constraint
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.cConstraint);
    
% --- Executes on button press in ref_path.
function ref_path_Callback(hObject, eventdata, handles)
% hObject    handle to ref_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ref_path
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.cRefPath);

% --- Executes on button press in plotproc.
function plotproc_Callback(hObject, eventdata, handles)
% hObject    handle to plotproc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotproc
global goGUI;
    goGUI.syncFromGUI(goGUI.idUI.cPlotProc);

% --- Executes on button press in no_skyplot_snr.
function no_skyplot_snr_Callback(hObject, eventdata, handles)
% hObject    handle to no_skyplot_snr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of no_skyplot_snr
global goGUI;
    goGUI.syncFromGUI(goGUI.idUI.cSkyPlot);

% --- Executes on button press in google_earth.
function google_earth_Callback(hObject, eventdata, handles)
% hObject    handle to google_earth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of google_earth
global goGUI;
    goGUI.syncFromGUI(goGUI.idUI.cGEarth);
    
% --- Executes on button press in err_ellipse.
function err_ellipse_Callback(hObject, eventdata, handles)
% hObject    handle to err_ellipse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of err_ellipse
global goGUI;
    goGUI.syncFromGUI(goGUI.idUI.cErrEllipse);
    
% --- Executes on button press in plot_master.
function plot_master_Callback(hObject, eventdata, handles)
% hObject    handle to plot_master (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_master
global goGUI;
    goGUI.syncFromGUI(goGUI.idUI.cPlotMaster);
    
% --- Executes on button press in plot_amb.
function plot_amb_Callback(hObject, eventdata, handles)
% hObject    handle to plot_amb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_amb
global goGUI;
    goGUI.syncFromGUI(goGUI.idUI.cPlotAmb);
    
% --- Executes on button press in use_ntrip.
function use_ntrip_Callback(hObject, eventdata, handles)
% hObject    handle to use_ntrip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_ntrip
global goGUI;
    goGUI.syncFromGUI(goGUI.idUI.cUseNTRIP);
    
% --- Executes on button press in flag_doppler.
function flag_doppler_Callback(hObject, eventdata, handles)
% hObject    handle to flag_doppler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of flag_doppler
global goGUI;
    goGUI.syncFromGUI(goGUI.idUI.cDoppler);
    
% --- Executes on button press in use_SBAS.
function use_SBAS_Callback(hObject, eventdata, handles)
% hObject    handle to use_SBAS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_SBAS
global goGUI;
    goGUI.syncFromGUI(goGUI.idUI.cUse_SBAS);
    
%   INPUT/OUTPUT FILE AND FOLDERS
% ===============================================================

% Rover/INI ---------------------------------------------------

function sINI_Callback(hObject, eventdata, handles)
% hObject    handle to sINI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sINI as text
%        str2double(get(hObject,'String')) returns contents of sINI as a double
global goGUI;
    goGUI.syncFromGUI(goGUI.idUI.sINI);

% --- Executes during object creation, after setting all properties.
function sINI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sINI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in bINI.
function bINI_Callback(hObject, eventdata, handles)
% hObject    handle to bINI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.bINI);

% --- Executes on button press in bEditINI.
function bEditINI_Callback(hObject, eventdata, handles)
% hObject    handle to bEditINI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global goGUI
    goGUI.openEditINI();

% Output ------------------------------------------------------    

function sDirGoOut_Callback(hObject, eventdata, handles)
% hObject    handle to sDirGoOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sDirGoOut as text
%        str2double(get(hObject,'String')) returns contents of sDirGoOut as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.sDirGoOut);

% --- Executes during object creation, after setting all properties.
function sDirGoOut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sDirGoOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in bDirGoOut.
function bDirGoOut_Callback(hObject, eventdata, handles)
% hObject    handle to bDirGoOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.bDirGoOut);
    
function sPrefixGoOut_Callback(hObject, eventdata, handles)
% hObject    handle to sPrefixGoOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sPrefixGoOut as text
%        str2double(get(hObject,'String')) returns contents of sPrefixGoOut as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.sPrefixGoOut);

% --- Executes during object creation, after setting all properties.
function sPrefixGoOut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sPrefixGoOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Constellation select ----------------------------------------

% --- Executes on button press in cGPS.
function cGPS_Callback(hObject, eventdata, handles)
% hObject    handle to cGPS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cGPS
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.cGPS);

% --- Executes on button press in cGLONASS.
function cGLONASS_Callback(hObject, eventdata, handles)
% hObject    handle to cGLONASS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.cGLONASS);

% --- Executes on button press in cGalileo.
function cGalileo_Callback(hObject, eventdata, handles)
% hObject    handle to cGalileo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cGalileo
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.cGalileo);

% --- Executes on button press in cBeiDou.
function cBeiDou_Callback(hObject, eventdata, handles)
% hObject    handle to cBeiDou (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cBeiDou
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.cBeiDou);

% --- Executes on button press in cQZSS.
function cQZSS_Callback(hObject, eventdata, handles)
% hObject    handle to cQZSS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cQZSS
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.cQZSS);

% --- Executes on button press in cSBAS.
function cSBAS_Callback(hObject, eventdata, handles)
% hObject    handle to cSBAS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cSBAS
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.cSBAS);

%   SETTINGS - KALMAN FILTER - STD
% ===============================================================

% East --------------------------------------------------------            

function std_X_Callback(hObject, eventdata, handles)
% hObject    handle to std_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_X as text
%        str2double(get(hObject,'String')) returns contents of std_X as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.nStdE);

% --- Executes during object creation, after setting all properties.
function std_X_CreateFcn(hObject, eventdata, handles)
% hObject    handle to std_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Nord --------------------------------------------------------            

function std_Y_Callback(hObject, eventdata, handles)
% hObject    handle to std_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_Y as text
%        str2double(get(hObject,'String')) returns contents of std_Y as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.nStdN);

% --- Executes during object creation, after setting all properties.
function std_Y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to std_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Up ----------------------------------------------------------

function std_Z_Callback(hObject, eventdata, handles)
% hObject    handle to std_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_Z as text
%        str2double(get(hObject,'String')) returns contents of std_Z as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.nStdU);

% --- Executes during object creation, after setting all properties.
function std_Z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to std_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Code --------------------------------------------------------

function std_code_Callback(hObject, eventdata, handles)
% hObject    handle to std_code (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_code as text
%        str2double(get(hObject,'String')) returns contents of std_code as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.nStdCode);

% --- Executes during object creation, after setting all properties.
function std_code_CreateFcn(hObject, eventdata, handles)
% hObject    handle to std_code (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Phase -------------------------------------------------------

% --- Executes on button press in toggle_std_phase.
function toggle_std_phase_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_std_phase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.bStdPhase);
    
function std_phase_Callback(hObject, eventdata, handles)
% hObject    handle to std_phase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_phase as text
%        str2double(get(hObject,'String')) returns contents of std_phase as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.nStdPhase);

% --- Executes during object creation, after setting all properties.
function std_phase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to std_phase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Initial State -----------------------------------------------

function std_init_Callback(hObject, eventdata, handles)
% hObject    handle to std_init (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_init as text
%        str2double(get(hObject,'String')) returns contents of std_init as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.nStdT0);

% --- Executes during object creation, after setting all properties.
function std_init_CreateFcn(hObject, eventdata, handles)
% hObject    handle to std_init (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% DTM ---------------------------------------------------------

% --- Executes on button press in toggle_std_dtm.
function toggle_std_dtm_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_std_dtm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggle_std_dtm
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.bStdDTM);
   
function std_dtm_Callback(hObject, eventdata, handles)
% hObject    handle to std_dtm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_dtm as text
%        str2double(get(hObject,'String')) returns contents of std_dtm as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.nStdDTM);

% --- Executes during object creation, after setting all properties.
function std_dtm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to std_dtm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Velocity ----------------------------------------------------

function std_vel_Callback(hObject, eventdata, handles)
% hObject    handle to std_vel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_vel as text
%        str2double(get(hObject,'String')) returns contents of std_vel as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.nStdVel);

% --- Executes during object creation, after setting all properties.
function std_vel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to std_vel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%   SETTINGS - KALMAN FILTER - WEIGHT MODEL
% ===============================================================

% --- Executes when selected object is changed in weight_select.
function weight_select_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in weight_select
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global goGUI
    goGUI.syncFromGUI(goGUI.idGroup.pW);

%   SETTINGS - KALMAN FILTER
% ===============================================================

% Cut-off thr -------------------------------------------------

function cut_off_Callback(hObject, eventdata, handles)
% hObject    handle to cut_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cut_off as text
%        str2double(get(hObject,'String')) returns contents of cut_off as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.nCutOff);

% --- Executes during object creation, after setting all properties.
function cut_off_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cut_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
    
% SNR thr -----------------------------------------------------

function snr_thres_Callback(hObject, eventdata, handles)
% hObject    handle to snr_thres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of snr_thres as text
%        str2double(get(hObject,'String')) returns contents of snr_thres as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.nSNR);

% --- Executes during object creation, after setting all properties.
function snr_thres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to snr_thres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% CS thr ------------------------------------------------------

function cs_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to cs_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cs_thresh as text
%        str2double(get(hObject,'String')) returns contents of cs_thresh as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.nCS);

% --- Executes during object creation, after setting all properties.
function cs_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cs_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Min sat number ----------------------------------------------

function min_sat_Callback(hObject, eventdata, handles)
% hObject    handle to min_sat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_sat as text
%        str2double(get(hObject,'String')) returns contents of min_sat as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.nMinNSat);


% --- Executes during object creation, after setting all properties.
function min_sat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_sat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Antenna height ----------------------------------------------

function antenna_h_Callback(hObject, eventdata, handles)
% hObject    handle to antenna_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of antenna_h as text
%        str2double(get(hObject,'String')) returns contents of antenna_h as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.nHAntenna);


% --- Executes during object creation, after setting all properties.
function antenna_h_CreateFcn(hObject, eventdata, handles)
% hObject    handle to antenna_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Stop Go Stop ------------------------------------------------

% --- Executes on button press in stopGOstop.
function stopGOstop_Callback(hObject, eventdata, handles)
% hObject    handle to stopGOstop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of stopGOstop
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.cStopGoStop);

%   SETTINGS - KALMAN FILTER - DYNAMIC MODEL
% ===============================================================

% --- Executes on selection change in dyn_mod.
function dyn_mod_Callback(hObject, eventdata, handles)
% hObject    handle to dyn_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns dyn_mod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dyn_mod
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.lDynModel);

% --- Executes during object creation, after setting all properties.
function dyn_mod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dyn_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%   SETTINGS - KALMAN FILTER - ARAA
% ===============================================================

% --- Executes on selection change in amb_select.
function amb_select_Callback(hObject, eventdata, handles)
% hObject    handle to amb_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns amb_select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from amb_select
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.lARAA);

% --- Executes during object creation, after setting all properties.
function amb_select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to amb_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%   SETTINGS - MASTER STATION
% ===============================================================

% --- Executes on button press in master_pos.
function master_pos_Callback(hObject, eventdata, handles)
% hObject    handle to master_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.cMPos);

% --- Executes on selection change in crs.
function crs_Callback(hObject, eventdata, handles)
% hObject    handle to crs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns crs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from crs
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.lCRS);

% --- Executes during object creation, after setting all properties.
function crs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to crs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function master_X_Callback(hObject, eventdata, handles)
% hObject    handle to master_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of master_X as text
%        str2double(get(hObject,'String')) returns contents of master_X as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.nMX);

% --- Executes during object creation, after setting all properties.
function master_X_CreateFcn(hObject, eventdata, handles)
% hObject    handle to master_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function master_Y_Callback(hObject, eventdata, handles)
% hObject    handle to master_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of master_Y as text
%        str2double(get(hObject,'String')) returns contents of master_Y as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.nMY);

% --- Executes during object creation, after setting all properties.
function master_Y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to master_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function master_Z_Callback(hObject, eventdata, handles)
% hObject    handle to master_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of master_Z as text
%        str2double(get(hObject,'String')) returns contents of master_Z as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.nMZ);

% --- Executes during object creation, after setting all properties.
function master_Z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to master_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function master_lat_Callback(hObject, eventdata, handles)
% hObject    handle to master_lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of master_lat as text
%        str2double(get(hObject,'String')) returns contents of master_lat as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.nMLat);

% --- Executes during object creation, after setting all properties.
function master_lat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to master_lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function master_lon_Callback(hObject, eventdata, handles)
% hObject    handle to master_lon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of master_lon as text
%        str2double(get(hObject,'String')) returns contents of master_lon as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.nMLon);

% --- Executes during object creation, after setting all properties.
function master_lon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to master_lon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function master_h_Callback(hObject, eventdata, handles)
% hObject    handle to master_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of master_h as text
%        str2double(get(hObject,'String')) returns contents of master_h as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.nMh);
    
% --- Executes during object creation, after setting all properties.
function master_h_CreateFcn(hObject, eventdata, handles)
% hObject    handle to master_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%   SETTINGS - PORTS
% ===============================================================

% --- Executes on selection change in num_receivers.
function num_receivers_Callback(hObject, eventdata, handles)
% hObject    handle to num_receivers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns num_receivers contents as cell array
%        contents{get(hObject,'Value')} returns selected item from num_receivers
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.lnPorts);
    
% --- Executes during object creation, after setting all properties.
function num_receivers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_receivers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function pumCaptureRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pumCaptureRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in com_select_0.
function com_select_0_Callback(hObject, eventdata, handles)
% hObject    handle to com_select_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns com_select_0 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from com_select_0
global goGUI
    goGUI.syncFromGUI(goGUI.idGroup.lPort0);

% --- Executes during object creation, after setting all properties.
function com_select_0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to com_select_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
    
% --- Executes on selection change in protocol_select_0.
function protocol_select_0_Callback(hObject, eventdata, handles)
% hObject    handle to protocol_select_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns protocol_select_0 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from protocol_select_0
global goGUI
    goGUI.syncFromGUI(goGUI.idGroup.lPort0);

% --- Executes during object creation, after setting all properties.
function protocol_select_0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to protocol_select_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in com_select_1.
function com_select_1_Callback(hObject, eventdata, handles)
% hObject    handle to com_select_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns com_select_1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from com_select_1
global goGUI
    goGUI.syncFromGUI(goGUI.idGroup.lPort1);
    
% --- Executes during object creation, after setting all properties.
function com_select_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to com_select_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in protocol_select_1.
function protocol_select_1_Callback(hObject, eventdata, handles)
% hObject    handle to protocol_select_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns protocol_select_1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from protocol_select_1
global goGUI
    goGUI.syncFromGUI(goGUI.idGroup.lPort1);

% --- Executes during object creation, after setting all properties.
function protocol_select_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to protocol_select_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in com_select_2.
function com_select_2_Callback(hObject, eventdata, handles)
% hObject    handle to com_select_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global goGUI
    goGUI.syncFromGUI(goGUI.idGroup.lPort2);
    
% --- Executes during object creation, after setting all properties.
function com_select_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to com_select_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in protocol_select_2.
function protocol_select_2_Callback(hObject, eventdata, handles)
% hObject    handle to protocol_select_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns protocol_select_2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from protocol_select_2
global goGUI
    goGUI.syncFromGUI(goGUI.idGroup.lPort2);

% --- Executes during object creation, after setting all properties.
function protocol_select_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to protocol_select_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in com_select_3.
function com_select_3_Callback(hObject, eventdata, handles)
% hObject    handle to com_select_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns com_select_3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from com_select_3
global goGUI
    goGUI.syncFromGUI(goGUI.idGroup.lPort3);
    
% --- Executes during object creation, after setting all properties.
function com_select_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to com_select_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in protocol_select_3.
function protocol_select_3_Callback(hObject, eventdata, handles)
% hObject    handle to protocol_select_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns protocol_select_3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from protocol_select_3
global goGUI
    goGUI.syncFromGUI(goGUI.idGroup.lPort3);

% --- Executes during object creation, after setting all properties.
function protocol_select_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to protocol_select_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%   SETTINGS - MASTER SERVER
% ===============================================================
            
function IP_address_Callback(hObject, eventdata, handles)
% hObject    handle to IP_address (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IP_address as text
%        str2double(get(hObject,'String')) returns contents of IP_address as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.sIPaddr);

% --- Executes during object creation, after setting all properties.
function IP_address_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IP_address (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function port_Callback(hObject, eventdata, handles)
% hObject    handle to port (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of port as text
%        str2double(get(hObject,'String')) returns contents of port as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.sIPport);

% --- Executes during object creation, after setting all properties.
function port_CreateFcn(hObject, eventdata, handles)
% hObject    handle to port (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mountpoint_Callback(hObject, eventdata, handles)
% hObject    handle to mountpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mountpoint as text
%        str2double(get(hObject,'String')) returns contents of mountpoint as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.sMnt);

% --- Executes during object creation, after setting all properties.
function mountpoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mountpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function username_Callback(hObject, eventdata, handles)
% hObject    handle to username (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of username as text
%        str2double(get(hObject,'String')) returns contents of username as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.sUName);

% --- Executes during object creation, after setting all properties.
function username_CreateFcn(hObject, eventdata, handles)
% hObject    handle to username (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function password_Callback(hObject, eventdata, handles)
% hObject    handle to password (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of password as text
%        str2double(get(hObject,'String')) returns contents of password as a double

% --- Executes during object creation, after setting all properties.
function password_CreateFcn(hObject, eventdata, handles)
% hObject    handle to password (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on key press with focus on password and none of its controls.
function password_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to password (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
global goGUI;
    key = eventdata.Key;
    ch = eventdata.Character;
    goGUI.modifyPassword(key, ch);
  
% --- Executes on button press in show_password.
function show_password_Callback(hObject, eventdata, handles)
% hObject    handle to show_password (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of show_password
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.bUPass);

function approx_lat_Callback(hObject, eventdata, handles)
% hObject    handle to approx_lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of approx_lat as text
%        str2double(get(hObject,'String')) returns contents of approx_lat as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.nVLat);

% --- Executes during object creation, after setting all properties.
function approx_lat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to approx_lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function approx_lon_Callback(hObject, eventdata, handles)
% hObject    handle to approx_lon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of approx_lon as text
%        str2double(get(hObject,'String')) returns contents of approx_lon as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.nVLon);

% --- Executes during object creation, after setting all properties.
function approx_lon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to approx_lon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function approx_h_Callback(hObject, eventdata, handles)
% hObject    handle to approx_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of approx_h as text
%        str2double(get(hObject,'String')) returns contents of approx_h as a double
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.nVH);

% --- Executes during object creation, after setting all properties.
function approx_h_CreateFcn(hObject, eventdata, handles)
% hObject    handle to approx_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in cLAMBDA.
function cLAMBDA_Callback(hObject, eventdata, handles)
% hObject    handle to cLAMBDA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cLAMBDA
global goGUI;
    goGUI.syncFromGUI(goGUI.idUI.cLAMBDA);

% --- Executes on button press in cMu.
function cMu_Callback(hObject, eventdata, handles)
% hObject    handle to cMu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cMu
global goGUI;
    goGUI.syncFromGUI(goGUI.idUI.cMu);
    
% --- Executes on selection change in lLAMBDAMethod.
function lLAMBDAMethod_Callback(hObject, eventdata, handles)
% hObject    handle to lLAMBDAMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lLAMBDAMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lLAMBDAMethod
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.lLAMBDAMethod);

% --- Executes during object creation, after setting all properties.
function lLAMBDAMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lLAMBDAMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nP0_Callback(hObject, eventdata, handles)
% hObject    handle to nP0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nP0 as text
%        str2double(get(hObject,'String')) returns contents of nP0 as a double
global goGUI;
    goGUI.syncFromGUI(goGUI.idUI.nP0);
    
% --- Executes during object creation, after setting all properties.
function nP0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nP0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nMu_Callback(hObject, eventdata, handles)
% hObject    handle to nMu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nMu as text
%        str2double(get(hObject,'String')) returns contents of nMu as a double
global goGUI;
    goGUI.syncFromGUI(goGUI.idUI.nMu);
    
% --- Executes during object creation, after setting all properties.
function nMu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nMu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in cP0.
function cP0_Callback(hObject, eventdata, handles)
% hObject    handle to cP0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cP0
global goGUI;
    goGUI.syncFromGUI(goGUI.idUI.cP0);




%   BUTTONS
% ===============================================================
        
% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.bExit);
    
% --- Executes on button press in load_button.
function load_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%allow the user to choose which settings to load
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.bLoad);
    
% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%allow the user to specify where to save the settings file
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.bSave);

% --- Executes on button press in go_button.
function go_button_Callback(hObject, eventdata, handles)
% hObject    handle to go_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global goGUI
    goGUI.syncFromGUI(goGUI.idUI.bGo);





    
%   MENU
% ===============================================================
% The management of the menu, is not yet in the goGUI object

% --------------------------------------------------------------------
function menu_tools_Callback(hObject, eventdata, handles)
% hObject    handle to menu_tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menu_about_Callback(hObject, eventdata, handles)
% hObject    handle to menu_about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function decode_streams_unix_Callback(hObject, eventdata, handles)
% hObject    handle to decode_streams_unix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_decode_stream_unix;

% --------------------------------------------------------------------
function merge_goGPS_bin_unix_Callback(hObject, eventdata, handles)
% hObject    handle to merge_goGPS_bin_unix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_merge_goGPSbin_unix;

% --------------------------------------------------------------------
function RINEX2goGPSbin_unix_Callback(hObject, eventdata, handles)
% hObject    handle to RINEX2goGPSbin_unix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_RINEX2goGPSbin_unix;

% --------------------------------------------------------------------
function polyline_unix_Callback(hObject, eventdata, handles)
% hObject    handle to polyline_unix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_polyline_simplification_unix;

% --------------------------------------------------------------------
function help_Callback(hObject, eventdata, handles)
% hObject    handle to help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function about_unix_Callback(hObject, eventdata, handles)
% hObject    handle to about_unix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_about_unix;


% --- Executes on key press with focus on main_panel and none of its controls.
function main_panel_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to main_panel (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
global goGUI
% Easter Egg: secret shortcut ALT+0 to test the interface
if length(eventdata.Modifier) == 1
    if strcmp(eventdata.Modifier{1},'alt') && strcmp(eventdata.Key,'0')
        goGUI.testOnOff();
    end
end
