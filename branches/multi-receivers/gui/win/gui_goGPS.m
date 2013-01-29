function varargout = gui_goGPS(varargin)
% GUI_GOGPS M-file for gui_goGPS.fig
%      GUI_GOGPS, by itself, creates a new GUI_GOGPS or raises the existing
%      singleton*.
%
%      H = GUI_GOGPS returns the handle to a new GUI_GOGPS or the handle to
%      the existing singleton*.
%
%      GUI_GOGPS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_GOGPS.M with the given input arguments.
%
%      GUI_GOGPS('Property','Value',...) creates a new GUI_GOGPS or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_goGPS_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_goGPS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_goGPS

% Last Modified by GUIDE v2.5 25-Apr-2011 20:05:52

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.1 beta
%
% Copyright (C) 2009-2012 Mirko Reguzzoni, Eugenio Realini
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
    'gui_OpeningFcn', @gui_goGPS_OpeningFcn, ...
    'gui_OutputFcn',  @gui_goGPS_OutputFcn, ...
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

% --- Executes just before gui_goGPS is made visible.
function gui_goGPS_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_goGPS (see VARARGIN)
clearvars -global goGUI
global goGUI
goGUI = goGUIclass(handles, goGUIclass.isWin);

% UIWAIT makes gui_goGPS wait for user response (see UIRESUME)
uiwait(hObject);


% --- Outputs from this function are returned to the command line.
function varargout = gui_goGPS_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global goGUI
    if(~isstruct(handles))
        varargout = cell(23,1);
        return
    end
    
    varargout = goGUI.outputFun();
%close main panel
delete(gcf)

% --------------------------------------------------------------------
function saveState(handles,filename)
state.master_pos = get(handles.master_pos, 'Value');
state.constraint = get(handles.constraint, 'Value');
state.ref_path = get(handles.ref_path, 'Value');
state.plot_master = get(handles.plot_master, 'Value');
state.google_earth = get(handles.google_earth, 'Value');
state.err_ellipse = get(handles.err_ellipse, 'Value');
state.use_ntrip = get(handles.use_ntrip, 'Value');
state.plot_amb = get(handles.plot_amb, 'Value');
state.no_skyplot_snr = get(handles.no_skyplot_snr, 'Value');
state.plotproc = get(handles.plotproc, 'Value');
state.RINEX_rover_obs = get(handles.RINEX_rover_obs,'String');
state.RINEX_master_obs = get(handles.RINEX_master_obs,'String');
state.RINEX_nav = get(handles.RINEX_nav,'String');
state.gogps_data_input = get(handles.gogps_data_input,'String');
state.gogps_data_output = get(handles.gogps_data_output,'String');
state.gogps_data_output_prefix = get(handles.gogps_data_output_prefix,'String');
state.dtm_path = get(handles.dtm_path,'String');
state.ref_path_input = get(handles.ref_path_input,'String');
state.master_X = get(handles.master_X,'String');
state.master_Y = get(handles.master_Y,'String');
state.master_Z = get(handles.master_Z,'String');
state.master_lat = get(handles.master_lat,'String');
state.master_lon = get(handles.master_lon,'String');
state.master_h = get(handles.master_h,'String');
state.crs = get(handles.crs,'Value');
state.mode = get(handles.mode,'Value');
state.nav_mon = get(handles.nav_mon,'Value');
state.kalman_ls = get(handles.kalman_ls,'Value');
state.code_dd_sa = get(handles.code_dd_sa,'Value');
state.rinex_files = get(handles.rinex_files,'Value');
state.gogps_data = get(handles.gogps_data,'Value');
state.std_X = get(handles.std_X,'String');
state.std_Y = get(handles.std_Y,'String');
state.std_Z = get(handles.std_Z,'String');
state.std_code = get(handles.std_code,'String');
state.std_phase = get(handles.std_phase,'String');
state.std_dtm = get(handles.std_dtm,'String');
state.toggle_std_phase = get(handles.toggle_std_phase,'Value');
state.toggle_std_dtm = get(handles.toggle_std_dtm,'Value');
state.std_init = get(handles.std_init,'String');
state.std_vel = get(handles.std_vel,'String');
state.cs_thresh = get(handles.cs_thresh,'String');
state.flag_doppler = get(handles.flag_doppler,'Value');
state.amb_select = get(handles.amb_select,'Value');
state.cut_off = get(handles.cut_off,'String');
state.snr_thres = get(handles.snr_thres,'String');
state.antenna_h = get(handles.antenna_h,'String');
state.min_sat = get(handles.min_sat,'String');
state.dyn_mod = get(handles.dyn_mod,'Value');
state.weight_0 = get(handles.weight_0,'Value');
state.weight_1 = get(handles.weight_1,'Value');
state.weight_2 = get(handles.weight_2,'Value');
state.weight_3 = get(handles.weight_3,'Value');
contents = cellstr(get(handles.com_select_0,'String'));
state.com_select_0 = contents{get(handles.com_select_0,'Value')};
contents = cellstr(get(handles.com_select_1,'String'));
state.com_select_1 = contents{get(handles.com_select_1,'Value')};
contents = cellstr(get(handles.com_select_2,'String'));
state.com_select_2 = contents{get(handles.com_select_2,'Value')};
contents = cellstr(get(handles.com_select_3,'String'));
state.com_select_3 = contents{get(handles.com_select_3,'Value')};
state.protocol_select_0 = get(handles.protocol_select_0,'Value');
state.protocol_select_1 = get(handles.protocol_select_1,'Value');
state.protocol_select_2 = get(handles.protocol_select_2,'Value');
state.protocol_select_3 = get(handles.protocol_select_3,'Value');
state.num_receivers = get(handles.num_receivers,'Value');
state.IP_address = get(handles.IP_address,'String');
state.port = get(handles.port,'String');
state.mountpoint = get(handles.mountpoint,'String');
state.username = get(handles.username,'String');
state.password = get(handles.password,'Userdata');
state.approx_lat = get(handles.approx_lat,'String');
state.approx_lon = get(handles.approx_lon,'String');
state.approx_h = get(handles.approx_h,'String');
state.stopGOstop = get(handles.stopGOstop,'Value');

save(filename, 'state');

% --- Executes on selection change in mode.
function mode_Callback(hObject, eventdata, handles)
% hObject    handle to mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns mode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mode
global goGUI
    goGUI.selectProcassingApproach();

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

% --- Executes on selection change in kalman_ls.
function kalman_ls_Callback(hObject, eventdata, handles)
% hObject    handle to kalman_ls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns kalman_ls contents as cell array
%        contents{get(hObject,'Value')} returns selected item from kalman_ls
%enable Kalman filters settings
global goGUI;
    goGUI.selectProcessingType();

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
    goGUI.selectProcessingSolution();

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


% --- Executes on selection change in nav_mon.
function nav_mon_Callback(hObject, eventdata, handles)
% hObject    handle to nav_mon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns nav_mon contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nav_mon
global goGUI
    goGUI.selectReceiverMode();


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

% --- Executes on button press in constraint.
function constraint_Callback(hObject, eventdata, handles)
% hObject    handle to constraint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of constraint
global goGUI
    goGUI.toggleLinearConstraint();

% --- Executes on button press in ref_path.
function ref_path_Callback(hObject, eventdata, handles)
% hObject    handle to ref_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ref_path
global goGUI
    goGUI.toggleReferencePath();

% --- Executes on button press in plot_master.
function plot_master_Callback(hObject, eventdata, handles)
% hObject    handle to plot_master (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_master


% --- Executes on button press in err_ellipse.
function err_ellipse_Callback(hObject, eventdata, handles)
% hObject    handle to err_ellipse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of err_ellipse


% --- Executes on button press in google_earth.
function google_earth_Callback(hObject, eventdata, handles)
% hObject    handle to google_earth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of google_earth

% --- Executes on button press in use_ntrip.
function use_ntrip_Callback(hObject, eventdata, handles)
% hObject    handle to use_ntrip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_ntrip
global goGUI;
    goGUI.toggleNTRIP();

% --- Executes on button press in master_pos.
function master_pos_Callback(hObject, eventdata, handles)
% hObject    handle to master_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global goGUI;
    goGUI.toggleMasterPosition();


function master_X_Callback(hObject, eventdata, handles)
% hObject    handle to master_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of master_X as text
%        str2double(get(hObject,'String')) returns contents of master_X as a double


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


% --- Executes on selection change in crs.
function crs_Callback(hObject, eventdata, handles)
% hObject    handle to crs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns crs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from crs
global goGUI
    goGUI.toggleMasterPositionType();

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

function gogps_data_input_Callback(hObject, eventdata, handles)
% hObject    handle to gogps_data_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gogps_data_input as text
%        str2double(get(hObject,'String')) returns contents of gogps_data_input as a double


% --- Executes during object creation, after setting all properties.
function gogps_data_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gogps_data_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_gogps_input.
function browse_gogps_input_Callback(hObject, eventdata, handles)
% hObject    handle to browse_gogps_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {'*.bin','goGPS binary data (*.bin)'}, ...
    'Choose goGPS binary data','../data');

if (filename ~= 0)
    pos = find(filename == '_');
    filename = filename(1:pos(end-1)-1);
    set(handles.gogps_data_input,'String',fullfile(pathname, filename));
end


function gogps_data_output_Callback(hObject, eventdata, handles)
% hObject    handle to gogps_data_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gogps_data_output as text
%        str2double(get(hObject,'String')) returns contents of gogps_data_output as a double


% --- Executes during object creation, after setting all properties.
function gogps_data_output_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gogps_data_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_gogps_data_output.
function browse_gogps_data_output_Callback(hObject, eventdata, handles)
% hObject    handle to browse_gogps_data_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dname = uigetdir('../data/','Choose a directory to store goGPS data');
if (dname ~= 0)
    set(handles.gogps_data_output,'String',dname);
end

function RINEX_rover_obs_Callback(hObject, eventdata, handles)
% hObject    handle to RINEX_rover_obs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RINEX_rover_obs as text
%        str2double(get(hObject,'String')) returns contents of RINEX_rover_obs as a double


% --- Executes during object creation, after setting all properties.
function RINEX_rover_obs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RINEX_rover_obs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_rover_obs.
function browse_rover_obs_Callback(hObject, eventdata, handles)
% hObject    handle to browse_rover_obs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {'*.obs;*.??o;*.??O','RINEX observation files (*.obs,*.??o,*.??O)';
    '*.obs','Observation files (*.obs)'; ...
    '*.??o;*.??O','Observation files (*.??o,*.??O)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose a RINEX observation file for the rover','../data/data_RINEX');

if (filename ~= 0)
    set(handles.RINEX_rover_obs,'String',fullfile(pathname, filename));
end


function RINEX_master_obs_Callback(hObject, eventdata, handles)
% hObject    handle to RINEX_master_obs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RINEX_master_obs as text
%        str2double(get(hObject,'String')) returns contents of RINEX_master_obs as a double


% --- Executes during object creation, after setting all properties.
function RINEX_master_obs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RINEX_master_obs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_master_obs.
function browse_master_obs_Callback(hObject, eventdata, handles)
% hObject    handle to browse_master_obs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {'*.obs;*.??o','RINEX observation files (*.obs,*.??o)';
    '*.obs','Observation files (*.obs)'; ...
    '*.??o','Observation files (*.??o)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose a RINEX observation file for the master','../data/data_RINEX');

if (filename ~= 0)
    set(handles.RINEX_master_obs,'String',fullfile(pathname, filename));
end


function RINEX_nav_Callback(hObject, eventdata, handles)
% hObject    handle to RINEX_nav (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RINEX_nav as text
%        str2double(get(hObject,'String')) returns contents of RINEX_nav as a double


% --- Executes during object creation, after setting all properties.
function RINEX_nav_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RINEX_nav (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_nav.
function browse_nav_Callback(hObject, eventdata, handles)
% hObject    handle to browse_nav (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {'*.nav;*.??n;*.??N','RINEX navigation files (*.nav,*.??n,*.??N)';
    '*.nav','Navigation files (*.nav)'; ...
    '*.??n;*.??N','Navigation files (*.??n,*.??N)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose a RINEX navigation file','../data/data_RINEX');

if (filename ~= 0)
    set(handles.RINEX_nav,'String',fullfile(pathname, filename));
end


% --- Executes when selected object is changed in file_type.
function file_type_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in file_type
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global goGUI;
    goGUI.selectFileType(obj.goh.gogps_data);



function gogps_data_output_prefix_Callback(hObject, eventdata, handles)
% hObject    handle to gogps_data_output_prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gogps_data_output_prefix as text
%        str2double(get(hObject,'String')) returns contents of gogps_data_output_prefix as a double


% --- Executes during object creation, after setting all properties.
function gogps_data_output_prefix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gogps_data_output_prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_amb.
function plot_amb_Callback(hObject, eventdata, handles)
% hObject    handle to plot_amb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_amb
global goGUI;
    goGUI.toggleAmbiguities();
    

function snr_thres_Callback(hObject, eventdata, handles)
% hObject    handle to snr_thres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of snr_thres as text
%        str2double(get(hObject,'String')) returns contents of snr_thres as a double


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



function cut_off_Callback(hObject, eventdata, handles)
% hObject    handle to cut_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cut_off as text
%        str2double(get(hObject,'String')) returns contents of cut_off as a double


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



function min_sat_Callback(hObject, eventdata, handles)
% hObject    handle to min_sat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_sat as text
%        str2double(get(hObject,'String')) returns contents of min_sat as a double


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



function cs_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to cs_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cs_thresh as text
%        str2double(get(hObject,'String')) returns contents of cs_thresh as a double


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


function std_vel_Callback(hObject, eventdata, handles)
% hObject    handle to std_vel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_vel as text
%        str2double(get(hObject,'String')) returns contents of std_vel as a double


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


% --- Executes on button press in toggle_std_phase.
function toggle_std_phase_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_std_phase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global goGUI;
    goGUI.togglePhaseStd();


function std_phase_Callback(hObject, eventdata, handles)
% hObject    handle to std_phase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_phase as text
%        str2double(get(hObject,'String')) returns contents of std_phase as a double


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



function std_code_Callback(hObject, eventdata, handles)
% hObject    handle to std_code (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_code as text
%        str2double(get(hObject,'String')) returns contents of std_code as a double


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


% --- Executes on button press in toggle_std_dtm.
function toggle_std_dtm_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_std_dtm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggle_std_dtm
global goGUI;
    goGUI.toggleDtmStd();


function std_dtm_Callback(hObject, eventdata, handles)
% hObject    handle to std_dtm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_dtm as text
%        str2double(get(hObject,'String')) returns contents of std_dtm as a double


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



function std_init_Callback(hObject, eventdata, handles)
% hObject    handle to std_init (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_init as text
%        str2double(get(hObject,'String')) returns contents of std_init as a double


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



function std_Z_Callback(hObject, eventdata, handles)
% hObject    handle to std_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_Z as text
%        str2double(get(hObject,'String')) returns contents of std_Z as a double


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



function std_Y_Callback(hObject, eventdata, handles)
% hObject    handle to std_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_Y as text
%        str2double(get(hObject,'String')) returns contents of std_Y as a double


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



function std_X_Callback(hObject, eventdata, handles)
% hObject    handle to std_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_X as text
%        str2double(get(hObject,'String')) returns contents of std_X as a double


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



function dtm_path_Callback(hObject, eventdata, handles)
% hObject    handle to dtm_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dtm_path as text
%        str2double(get(hObject,'String')) returns contents of dtm_path as a double


% --- Executes during object creation, after setting all properties.
function dtm_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dtm_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_dtm_path.
function browse_dtm_path_Callback(hObject, eventdata, handles)
% hObject    handle to browse_dtm_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dname = uigetdir('../data/dtm','Choose a directory containing DTM data');
if (dname ~= 0)
    set(handles.dtm_path,'String',dname);
end

function ref_path_input_Callback(hObject, eventdata, handles)
% hObject    handle to ref_path_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ref_path_input as text
%        str2double(get(hObject,'String')) returns contents of ref_path_input as a double


% --- Executes during object creation, after setting all properties.
function ref_path_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ref_path_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_ref_path_input.
function browse_ref_path_input_Callback(hObject, eventdata, handles)
% hObject    handle to browse_ref_path_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile('*.mat', 'Choose file containing reference path','../data');

if (filename ~= 0)
    set(handles.ref_path_input,'String',fullfile(pathname, filename));
end


function IP_address_Callback(hObject, eventdata, handles)
% hObject    handle to IP_address (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IP_address as text
%        str2double(get(hObject,'String')) returns contents of IP_address as a double


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



function username_Callback(hObject, eventdata, handles)
% hObject    handle to username (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of username as text
%        str2double(get(hObject,'String')) returns contents of username as a double


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



function mountpoint_Callback(hObject, eventdata, handles)
% hObject    handle to mountpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mountpoint as text
%        str2double(get(hObject,'String')) returns contents of mountpoint as a double


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


% --- Executes on selection change in dyn_mod.
function dyn_mod_Callback(hObject, eventdata, handles)
% hObject    handle to dyn_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns dyn_mod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dyn_mod
global goGUI
    goGUI.selectDynMode();


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


function approx_lat_Callback(hObject, eventdata, handles)
% hObject    handle to approx_lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of approx_lat as text
%        str2double(get(hObject,'String')) returns contents of approx_lat as a double


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



function antenna_h_Callback(hObject, eventdata, handles)
% hObject    handle to antenna_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of antenna_h as text
%        str2double(get(hObject,'String')) returns contents of antenna_h as a double


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


% --- Executes on button press in go_button.
function go_button_Callback(hObject, eventdata, handles)
% hObject    handle to go_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global goObj goGUI

%master station coordinates
crs_contents = cellstr(get(handles.crs,'String'));
master_X = str2double(get(handles.master_X,'String'));
master_Y = str2double(get(handles.master_Y,'String'));
master_Z = str2double(get(handles.master_Z,'String'));
master_lat = str2double(get(handles.master_lat,'String'));
master_lon = str2double(get(handles.master_lon,'String'));
master_h = str2double(get(handles.master_h,'String'));
%KF parameters
std_init = str2double(get(handles.std_init,'String'));
std_X = str2double(get(handles.std_X,'String'));
std_Y = str2double(get(handles.std_Y,'String'));
std_Z = str2double(get(handles.std_Z,'String'));
std_vel = str2double(get(handles.std_vel,'String'));
std_code = str2double(get(handles.std_code,'String'));
if (get(handles.toggle_std_phase,'Value'))
    std_phase = str2double(get(handles.std_phase,'String'));
else
    std_phase = 1e30;
end
if (get(handles.toggle_std_dtm,'Value'))
    std_dtm = str2double(get(handles.std_dtm,'String'));
else
    std_dtm = 1e30;
end
min_nsat = str2double(get(handles.min_sat,'String'));
cutoff = str2double(get(handles.cut_off,'String'));
snr_threshold = str2double(get(handles.snr_thres,'String'));
cs_threshold = str2double(get(handles.cs_thresh,'String'));
antenna_h = str2double(get(handles.antenna_h,'String'));
contents_dyn_mod = cellstr(get(handles.dyn_mod,'String'));
flag_stopGOstop = get(handles.stopGOstop,'Value');
%input files
filerootIN = get(handles.gogps_data_input,'String');
filerootOUT = [get(handles.gogps_data_output,'String') '\' get(handles.gogps_data_output_prefix,'String')];
filerootIN(filerootIN == '\') = '/';
filerootOUT(filerootOUT == '\') = '/';
filename_R_obs = get(handles.RINEX_rover_obs,'String');
filename_M_obs = get(handles.RINEX_master_obs,'String');
filename_nav = get(handles.RINEX_nav,'String');
ref_path = get(handles.ref_path, 'Value');
filename_ref = get(handles.ref_path_input,'String');
dtm_dir = get(handles.dtm_path,'String');
%serial communication
% global COMportR
contents = cellstr(get(handles.com_select_0,'String'));
COMportR0 = contents{get(handles.com_select_0,'Value')};
% contents = cellstr(get(handles.com_select_1,'String'));
% COMportR1 = contents{get(handles.com_select_1,'Value')};
% contents = cellstr(get(handles.com_select_2,'String'));
% COMportR2 = contents{get(handles.com_select_2,'Value')};
% contents = cellstr(get(handles.com_select_3,'String'));
% COMportR3 = contents{get(handles.com_select_3,'Value')};
%TCPIP / NTRIP
flag_NTRIP = get(handles.use_ntrip,'Value');
master_ip = get(handles.IP_address,'String');
master_port = str2double(get(handles.port,'String'));
ntrip_mountpoint = get(handles.mountpoint,'String');
%functioning mode
mode = goGUI.getMode();

% If I'm in a mode that uses objects instead of regular code, set goObj flag to 1
if (mode == 15)
	goObj = true;    
    % init the objects:    
else
    goObj = false;	% set to 1 when goGPS objects are used instead of the regular code
end

ready = 1;

%check if everything is OK before starting
if isempty(dir(get(handles.gogps_data_output,'String')))
    msgbox('Output folder does not exist. Please browse to an existing folder.'); ready = 0;
end

if ~isempty(dir(get(handles.gogps_data_output,'String')))
    i = 1;
    j = length(filerootOUT);
    while (~isempty(dir([filerootOUT '_????_rover.bin'])) | ...
            ~isempty(dir([filerootOUT '_master*.bin'])) | ...
            ~isempty(dir([filerootOUT '_????_obs*.bin'])) | ...
            ~isempty(dir([filerootOUT '_????_eph*.bin'])) | ...
            ~isempty(dir([filerootOUT '_????_dyn*.bin'])) | ...
            ~isempty(dir([filerootOUT '_sat*.bin'])) | ...
            ~isempty(dir([filerootOUT '_kal*.bin'])) | ...
            ~isempty(dir([filerootOUT '_dt*.bin'])) | ...
            ~isempty(dir([filerootOUT '_conf*.bin'])) | ...
            ~isempty(dir([filerootOUT '_dop*.bin'])) | ...
            ~isempty(dir([filerootOUT '_ECEF*.txt'])) | ...
            ~isempty(dir([filerootOUT '_geod*.txt'])) | ...
            ~isempty(dir([filerootOUT '_plan*.txt'])) | ...
            ~isempty(dir([filerootOUT '_????_NMEA*.txt'])) | ...
            ~isempty(dir([filerootOUT '.kml'])) )
        
        if (i == 100)
            msgbox('Automatic numbering of output file has reached the maximum (99). Please provide a different output prefix, select a different output folder or remove some files from the current output folder.');
            ready = 0;
            break
        end
        
        filerootOUT(j+1:j+3) = ['_' num2str(i,'%02d')];
        i = i + 1;
    end
end

if (goObj) && (get(handles.file_type, 'SelectedObject') == handles.rinex_files & isempty(dir(filename_R_obs)))
        msgbox('Ini file not found.'); ready = 0;
end
    
if (mode <= 20) && (~goObj) %if post-processing
    if (get(handles.file_type, 'SelectedObject') == handles.gogps_data & (isempty(dir([filerootIN '_obs*.bin'])) | isempty(dir([filerootIN '_eph*.bin']))))
        msgbox('Input goGPS binary files not found (both *_obs_* and *_eph_* files are needed).'); ready = 0;
    elseif (get(handles.file_type, 'SelectedObject') == handles.rinex_files & isempty(dir(filename_R_obs)))
        msgbox('Input rover observation RINEX file not found.'); ready = 0;
    elseif (mode > 10 & get(handles.file_type, 'SelectedObject') == handles.rinex_files & isempty(dir(filename_M_obs)))
        msgbox('Input master observation RINEX file not found.'); ready = 0;
    elseif (get(handles.file_type, 'SelectedObject') == handles.rinex_files & isempty(dir(filename_nav)))
        msgbox('Input navigation RINEX file not found.'); ready = 0;
    end
end

if (mode <= 20 | mode == 24) && (~goObj) %if not rover and/or master monitor
    if (ref_path & isempty(dir(filename_ref)))
        msgbox('Input reference path file not found.'); ready = 0;
    elseif ((mode == 14 | mode == 4) & get(handles.toggle_std_dtm,'Value'))
        try
            load([dtm_dir '/tiles/tile_header'], 'tile_header');
            load([dtm_dir '/tiles/tile_georef'], 'tile_georef');
        catch
            msgbox('DTM directory does not contain valid data.'); ready = 0;
        end
    elseif (mode > 10 & strcmp(get(handles.crs,'Enable'),'on') & ...
            strcmp(crs_contents{get(handles.crs,'Value')},'ECEF (X,Y,Z)') & ...
            ((isnan(master_X) | isnan(master_Y) | isnan(master_Z)) | ...
            (master_X == 0 & master_Y == 0 & master_Z == 0)))
        msgbox('Please provide valid values for the master station ECEF coordinates.'); ready = 0;
    elseif (mode > 10 & strcmp(get(handles.crs,'Enable'),'on') & ...
            strcmp(crs_contents{get(handles.crs,'Value')},'Geodetic coord. (lat,lon,h)') & ...
            ((isnan(master_lat) | isnan(master_lon) | isnan(master_h)) | ...
            (abs(master_lat) > 90 | abs(master_lon) > 180)))
        msgbox('Please provide valid values for the master station geodetic coordinates.'); ready = 0;
    elseif (isnan(std_init) | std_init < 0)
        msgbox('Please provide a valid value for the initial state error standard deviation.'); ready = 0;
    elseif (isnan(std_X) | std_X < 0)
        msgbox('Please provide a valid value for the East coordinate error standard deviation.'); ready = 0;
    elseif (isnan(std_Y) | std_Y < 0)
        msgbox('Please provide a valid value for the North coordinate error standard deviation.'); ready = 0;
    elseif (isnan(std_Z) | std_Z < 0)
        msgbox('Please provide a valid value for the Up coordinate error standard deviation.'); ready = 0;
    elseif (isnan(std_vel) | std_vel < 0)
        msgbox('Please provide a valid value for the velocity coordinate error standard deviation.'); ready = 0;
    elseif (isnan(std_code) | std_code < 0)
        msgbox('Please provide a valid value for the code error standard deviation.'); ready = 0;
    elseif (isnan(std_phase) | std_phase < 0)
        msgbox('Please provide a valid value for the phase error standard deviation.'); ready = 0;
    elseif (isnan(std_dtm) | std_dtm < 0)
        msgbox('Please provide a valid value for the DTM error standard deviation.'); ready = 0;
    elseif (isnan(min_nsat) | min_nsat < 0)
        msgbox('Please provide a valid value for the minimum number of satellites.'); ready = 0;
    elseif (isnan(cutoff) | cutoff < 0 | cutoff > 90)
        msgbox('Please provide a valid value for the cutoff (between 0 and 90 degrees).'); ready = 0;
    elseif (isnan(snr_threshold) | snr_threshold < 0 | snr_threshold > 60)
        msgbox('Please provide a valid value for the SNR threshold (between 0 and 60 dB).'); ready = 0;
    elseif (isnan(cs_threshold) | cs_threshold < 0)
        msgbox('Please provide a valid value for the cycle slip threshold.'); ready = 0;
    elseif (isnan(antenna_h) | antenna_h < 0)
        msgbox('Please provide a valid value for the antenna height.'); ready = 0;
    end
end

%check if the dataset was surveyed with a variable dynamic model
d = dir([filerootIN '_dyn_00.bin']);
if (mode <= 20 & (flag_stopGOstop | strcmp(contents_dyn_mod{get(handles.dyn_mod,'Value')},'Variable')) & isempty(d))
    msgbox('The selected dataset was not surveyed with a variable dynamic model: please select another dynamic model.'); ready = 0;
end

if (mode == 21 | mode == 23 | mode == 24) %if a COM connection to the rover is required
    if(strcmp(COMportR0, 'NA'))
        msgbox('Please select an existing COM port.'); ready = 0;
    end
end

if (mode == 22 | mode == 23 | mode == 24) %if a TCP/IP connection to the master is required
    if (isempty(master_ip))
        msgbox('Please provide an IP address for the connection to the master.'); ready = 0;
    elseif (isnan(master_port) | master_port < 0 | master_port > 65535)
        msgbox('Please provide a valid port number for the connection to the master (between 0 and 65535).'); ready = 0;
    end
    if (flag_NTRIP) %if a NTRIP connection is required
        if (isempty(ntrip_mountpoint))
            msgbox('Please provide a mountpoint for the NTRIP connection.'); ready = 0;
        end
    end
end

if (ready)
    saveState(handles,'../data/settings/last_settings.mat');
    uiresume(handles.main_panel);
end

% --- Executes on button press in load_button.
function load_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%allow the user to choose which settings to load
global goGUI
[filename, pathname] = uigetfile('*.mat', 'Choose file with saved settings','../data/settings');

if pathname == 0 %if the user pressed cancelled, then we exit this callback
    return
end

%construct the path name of the file to be loaded
loadDataName = fullfile(pathname,filename);

%load the settings, which creates a new gui
goGUI.loadState(loadDataName);

% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%allow the user to specify where to save the settings file
[filename,pathname] = uiputfile('*.mat','Save your GUI settings','../data/settings');

if pathname == 0 %if the user pressed cancelled, then we exit this callback
    return
end
%construct the path name of the save location
saveDataName = fullfile(pathname,filename);

%saves the gui data
saveState(handles, saveDataName);

% --- Executes when selected object is changed in weight_select.
function weight_select_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in weight_select
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in com_select_0.
function com_select_0_Callback(hObject, eventdata, handles)
% hObject    handle to com_select_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns com_select_0 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from com_select_0
% contents1 = cellstr(get(handles.com_select_1,'String'));
% contents2 = cellstr(get(handles.com_select_2,'String'));
% contents3 = cellstr(get(handles.com_select_3,'String'));
% serialInfo = instrhwinfo('serial');
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents1{get(handles.com_select_1,'Value')})) = [];
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents2{get(handles.com_select_2,'Value')})) = [];
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents3{get(handles.com_select_3,'Value')})) = [];
% set(hObject, 'String', serialInfo.AvailableSerialPorts);

% contents  = cellstr(get(hObject,'String'));
% COMportR0 = contents{get(hObject,'Value')};



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
try
    serialInfo = instrhwinfo('serial');
catch
    serialInfo.AvailableSerialPorts = [];
end
if (~isempty(serialInfo.AvailableSerialPorts))
    set(hObject, 'String', serialInfo.AvailableSerialPorts);
else
    notAvailable = cell(1);
    notAvailable{1} = 'NA';
    set(hObject, 'String', notAvailable);
end

% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)


% --- Executes on button press in no_skyplot_snr.
function no_skyplot_snr_Callback(hObject, eventdata, handles)
% hObject    handle to no_skyplot_snr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of no_skyplot_snr


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
    goGUI.showPassword();

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
function decode_streams_Callback(hObject, eventdata, handles)
% hObject    handle to decode_streams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_decode_stream;


% --------------------------------------------------------------------
function merge_goGPS_bin_Callback(hObject, eventdata, handles)
% hObject    handle to merge_goGPS_bin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_merge_goGPSbin;


% --------------------------------------------------------------------
function RINEX2goGPSbin_Callback(hObject, eventdata, handles)
% hObject    handle to RINEX2goGPSbin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_RINEX2goGPSbin;


% --- Executes on button press in plotproc.
function plotproc_Callback(hObject, eventdata, handles)
% hObject    handle to plotproc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotproc
global goGUI;
    goGUI.togglePlotOnProcessing();


% --- Executes on button press in stopGOstop.
function stopGOstop_Callback(hObject, eventdata, handles)
% hObject    handle to stopGOstop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of stopGOstop
global goGUI;
    goGUI.toggleStopGoStop();


% --------------------------------------------------------------------
function polyline_Callback(hObject, eventdata, handles)
% hObject    handle to polyline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_polyline_simplification;


% --------------------------------------------------------------------
function help_Callback(hObject, eventdata, handles)
% hObject    handle to help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function about_Callback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_about;


% --- Executes on selection change in protocol_select_0.
function protocol_select_0_Callback(hObject, eventdata, handles)
% hObject    handle to protocol_select_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns protocol_select_0 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from protocol_select_0


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

% --- Executes on button press in flag_doppler.
function flag_doppler_Callback(hObject, eventdata, handles)
% hObject    handle to flag_doppler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of flag_doppler


% --- Executes on selection change in amb_select.
function amb_select_Callback(hObject, eventdata, handles)
% hObject    handle to amb_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns amb_select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from amb_select
global goGUI
    goGUI = selectAmbiguityRestartMethod();

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


% --- Executes on selection change in com_select_1.
function com_select_1_Callback(hObject, eventdata, handles)
% hObject    handle to com_select_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns com_select_1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from com_select_1

% contents0 = cellstr(get(handles.com_select_0,'String'));
% contents2 = cellstr(get(handles.com_select_2,'String'));
% contents3 = cellstr(get(handles.com_select_3,'String'));
% serialInfo = instrhwinfo('serial');
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents0{get(handles.com_select_0,'Value')})) = [];
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents2{get(handles.com_select_2,'Value')})) = [];
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents3{get(handles.com_select_3,'Value')})) = [];
% set(hObject, 'String', serialInfo.AvailableSerialPorts);

% contents1  = cellstr(get(hObject,'String'));
% COMportR1 = contents1{get(hObject,'Value')};

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
try
    serialInfo = instrhwinfo('serial');
catch
    serialInfo.AvailableSerialPorts = [];
end
if (~isempty(serialInfo.AvailableSerialPorts))
%     serialInfo.AvailableSerialPorts(1) = [];
%     serialInfo.AvailableSerialPorts(3) = [];
%     serialInfo.AvailableSerialPorts(4) = [];
    set(hObject, 'String', serialInfo.AvailableSerialPorts);
end
if (isempty(serialInfo.AvailableSerialPorts))
    notAvailable = cell(1);
    notAvailable{1} = 'NA';
    set(hObject, 'String', notAvailable);
end


% --- Executes on selection change in protocol_select_1.
function protocol_select_1_Callback(hObject, eventdata, handles)
% hObject    handle to protocol_select_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns protocol_select_1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from protocol_select_1


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

% Hints: contents = cellstr(get(hObject,'String')) returns com_select_2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from com_select_2
% contents0 = cellstr(get(handles.com_select_0,'String'));
% contents1 = cellstr(get(handles.com_select_1,'String'));
% contents3 = cellstr(get(handles.com_select_3,'String'));
% serialInfo = instrhwinfo('serial');
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents0{get(handles.com_select_0,'Value')})) = [];
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents1{get(handles.com_select_1,'Value')})) = [];
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents3{get(handles.com_select_3,'Value')})) = [];
% set(hObject, 'String', serialInfo.AvailableSerialPorts);

% contents  = cellstr(get(hObject,'String'));
% COMportR2 = contents{get(hObject,'Value')};

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
try
    serialInfo = instrhwinfo('serial');
catch
    serialInfo.AvailableSerialPorts = [];
end
if (~isempty(serialInfo.AvailableSerialPorts))
%     serialInfo.AvailableSerialPorts(1) = [];
%     serialInfo.AvailableSerialPorts(2) = [];
%     serialInfo.AvailableSerialPorts(4) = [];
    set(hObject, 'String', serialInfo.AvailableSerialPorts);
end
if (isempty(serialInfo.AvailableSerialPorts))
    notAvailable = cell(1);
    notAvailable{1} = 'NA';
    set(hObject, 'String', notAvailable);
end


% --- Executes on selection change in protocol_select_2.
function protocol_select_2_Callback(hObject, eventdata, handles)
% hObject    handle to protocol_select_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns protocol_select_2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from protocol_select_2


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
% contents0 = cellstr(get(handles.com_select_0,'String'));
% contents1 = cellstr(get(handles.com_select_1,'String'));
% contents2 = cellstr(get(handles.com_select_2,'String'));
% serialInfo = instrhwinfo('serial');
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents0{get(handles.com_select_0,'Value')})) = [];
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents1{get(handles.com_select_1,'Value')})) = [];
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents2{get(handles.com_select_2,'Value')})) = [];
% set(hObject, 'String', serialInfo.AvailableSerialPorts);

% contents  = cellstr(get(hObject,'String'));
% COMportR3 = contents{get(hObject,'Value')};


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
try
    serialInfo = instrhwinfo('serial');
catch
    serialInfo.AvailableSerialPorts = [];
end
if (~isempty(serialInfo.AvailableSerialPorts))
%     serialInfo.AvailableSerialPorts(1) = [];
%     serialInfo.AvailableSerialPorts(2) = [];
%     serialInfo.AvailableSerialPorts(3) = [];
    set(hObject, 'String', serialInfo.AvailableSerialPorts);
end
if (isempty(serialInfo.AvailableSerialPorts))
    notAvailable = cell(1);
    notAvailable{1} = 'NA';
    set(hObject, 'String', notAvailable);
end


% --- Executes on selection change in protocol_select_3.
function protocol_select_3_Callback(hObject, eventdata, handles)
% hObject    handle to protocol_select_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns protocol_select_3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from protocol_select_3


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


% --- Executes on selection change in num_receivers.
function num_receivers_Callback(hObject, eventdata, handles)
% hObject    handle to num_receivers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns num_receivers contents as cell array
%        contents{get(hObject,'Value')} returns selected item from num_receivers
global goGUI
    goGUI.setNumReceiverIn();
    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------


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

contents = {'1';'2';'3';'4'};
try
    serialInfo = instrhwinfo('serial');
    num_ports = size(serialInfo.AvailableSerialPorts,1);
catch
    num_ports = 0;
end
if num_ports == 0
    set(hObject,'String','1');
elseif num_ports <= size(contents,1);
    set(hObject,'String',contents(1:num_ports));
else
    set(hObject,'String',contents);
end
