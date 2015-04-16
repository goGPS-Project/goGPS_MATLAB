function varargout = gui_RINEX2goGPSbin_unix(varargin)
% GUI_RINEX2GOGPSBIN_UNIX M-file for gui_RINEX2goGPSbin_unix.fig
%      GUI_RINEX2GOGPSBIN_UNIX, by itself, creates a new GUI_RINEX2GOGPSBIN_UNIX or raises the existing
%      singleton*.
%
%      H = GUI_RINEX2GOGPSBIN_UNIX returns the handle to a new GUI_RINEX2GOGPSBIN_UNIX or the handle to
%      the existing singleton*.
%
%      GUI_RINEX2GOGPSBIN_UNIX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_RINEX2GOGPSBIN_UNIX.M with the given input arguments.
%
%      GUI_RINEX2GOGPSBIN_UNIX('Property','Value',...) creates a new GUI_RINEX2GOGPSBIN_UNIX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_RINEX2goGPSbin_unix_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_RINEX2goGPSbin_unix_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_RINEX2goGPSbin_unix

% Last Modified by GUIDE v2.5 15-Jun-2010 15:53:51

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
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
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_RINEX2goGPSbin_unix_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_RINEX2goGPSbin_unix_OutputFcn, ...
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


% --- Executes just before gui_RINEX2goGPSbin_unix is made visible.
function gui_RINEX2goGPSbin_unix_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_RINEX2goGPSbin_unix (see VARARGIN)

% Choose default command line output for gui_RINEX2goGPSbin_unix
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%pixels
set(hObject, 'Units', 'pixels' );

%get display size
screenSize = get(0, 'ScreenSize');

%calculate the center of the display
position = get(hObject, 'Position');
position(1) = (screenSize(3)-position(3))/2;
position(2) = (screenSize(4)-position(4))/2;

%center the window
set(hObject, 'Position', position);

% --- Outputs from this function are returned to the command line.
function varargout = gui_RINEX2goGPSbin_unix_OutputFcn(hObject, eventdata, handles)  %#ok<*STOUT,*INUSD>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure



function data_out_folder_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to data_out_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of data_out_folder as text
%        str2double(get(hObject,'String')) returns contents of data_out_folder as a double


% --- Executes during object creation, after setting all properties.
function data_out_folder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_out_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_data_out_folder.
function browse_data_out_folder_Callback(hObject, eventdata, handles)
% hObject    handle to browse_data_out_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dname = uigetdir('../data/','Choose a directory to store output data');
if (dname ~= 0)
    set(handles.data_out_folder,'String',dname);
end

function data_out_Callback(hObject, eventdata, handles)
% hObject    handle to data_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of data_out as text
%        str2double(get(hObject,'String')) returns contents of data_out as a double


% --- Executes during object creation, after setting all properties.
function data_out_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in convert_button.
function convert_button_Callback(hObject, eventdata, handles)
% hObject    handle to convert_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename_R_obs = get(handles.rover_obs_file,'String');
filename_R_nav = get(handles.rover_nav_file,'String');
filename_M_obs = get(handles.master_obs_file,'String');
filename_M_nav = get(handles.master_nav_file,'String');
filerootOUT = [get(handles.data_out_folder,'String') '\' get(handles.data_out,'String')];
filename_R_obs(filename_R_obs == '\') = '/';
filename_R_nav(filename_R_nav == '\') = '/';
filename_M_obs(filename_M_obs == '\') = '/';
filename_M_nav(filename_M_nav == '\') = '/';
filerootOUT(filerootOUT == '\') = '/';

wait_dlg = waitbar(0,'Please wait...');

RINEX2goGPSbin(filename_R_obs, filename_R_nav, filename_M_obs, filename_M_nav, filerootOUT, wait_dlg);

close(wait_dlg)

close(handles.RINEX2goGPS_panel)

% --- Executes on button press in exit_button.
function exit_button_Callback(hObject, eventdata, handles)
% hObject    handle to exit_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.RINEX2goGPS_panel)


function master_obs_file_Callback(hObject, eventdata, handles)
% hObject    handle to master_obs_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of master_obs_file as text
%        str2double(get(hObject,'String')) returns contents of master_obs_file as a double


% --- Executes during object creation, after setting all properties.
function master_obs_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to master_obs_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_master_obs_file.
function browse_master_obs_file_Callback(hObject, eventdata, handles)
% hObject    handle to browse_master_obs_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {'*.obs;*.??o','RINEX observation files (*.obs,*.??o)';
    '*.obs','Observation files (*.obs)'; ...
    '*.??o','Observation files (*.??o)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose a RINEX observation file for the master','../data/data_RINEX');

if (filename ~= 0)
    set(handles.master_obs_file,'String',fullfile(pathname, filename));
end


function master_nav_file_Callback(hObject, eventdata, handles)
% hObject    handle to master_nav_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of master_nav_file as text
%        str2double(get(hObject,'String')) returns contents of master_nav_file as a double


% --- Executes during object creation, after setting all properties.
function master_nav_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to master_nav_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_master_nav_file.
function browse_master_nav_file_Callback(hObject, eventdata, handles)
% hObject    handle to browse_master_nav_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {'*.nav;*.??n','RINEX navigation files (*.nav,*.??n)';
    '*.nav','Navigation files (*.obs)'; ...
    '*.??n','Navigation files (*.??o)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose a RINEX navigation file for the master','../data/data_RINEX');

if (filename ~= 0)
    set(handles.master_nav_file,'String',fullfile(pathname, filename));
    if isempty(get(handles.rover_nav_file,'String'))
        set(handles.rover_nav_file,'String',fullfile(pathname, filename));
    end
end


function rover_obs_file_Callback(hObject, eventdata, handles)
% hObject    handle to rover_obs_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rover_obs_file as text
%        str2double(get(hObject,'String')) returns contents of rover_obs_file as a double


% --- Executes during object creation, after setting all properties.
function rover_obs_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rover_obs_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_rover_obs_file.
function browse_rover_obs_file_Callback(hObject, eventdata, handles)
% hObject    handle to browse_rover_obs_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {'*.obs;*.??o','RINEX observation files (*.obs,*.??o)';
    '*.obs','Observation files (*.obs)'; ...
    '*.??o','Observation files (*.??o)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose a RINEX observation file for the rover','../data/data_RINEX');

if (filename ~= 0)
    set(handles.rover_obs_file,'String',fullfile(pathname, filename));
end


function rover_nav_file_Callback(hObject, eventdata, handles)
% hObject    handle to rover_nav_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rover_nav_file as text
%        str2double(get(hObject,'String')) returns contents of rover_nav_file as a double


% --- Executes during object creation, after setting all properties.
function rover_nav_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rover_nav_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_rover_nav_file.
function browse_rover_nav_file_Callback(hObject, eventdata, handles)
% hObject    handle to browse_rover_nav_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {'*.nav;*.??n','RINEX navigation files (*.nav,*.??n)';
    '*.nav','Navigation files (*.obs)'; ...
    '*.??n','Navigation files (*.??o)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose a RINEX navigation file for the rover','../data/data_RINEX');

if (filename ~= 0)
    set(handles.rover_nav_file,'String',fullfile(pathname, filename));
    if isempty(get(handles.master_nav_file,'String'))
        set(handles.master_nav_file,'String',fullfile(pathname, filename));
    end
end
