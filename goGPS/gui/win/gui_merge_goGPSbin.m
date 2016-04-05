function varargout = gui_merge_goGPSbin(varargin)
% GUI_MERGE_GOGPSBIN M-file for gui_merge_goGPSbin.fig
%      GUI_MERGE_GOGPSBIN, by itself, creates a new GUI_MERGE_GOGPSBIN or raises the existing
%      singleton*.
%
%      H = GUI_MERGE_GOGPSBIN returns the handle to a new GUI_MERGE_GOGPSBIN or the handle to
%      the existing singleton*.
%
%      GUI_MERGE_GOGPSBIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_MERGE_GOGPSBIN.M with the given input arguments.
%
%      GUI_MERGE_GOGPSBIN('Property','Value',...) creates a new GUI_MERGE_GOGPSBIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_merge_goGPSbin_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_merge_goGPSbin_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_merge_goGPSbin

% Last Modified by GUIDE v2.5 28-Jun-2010 15:16:10

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
gui_State = struct('gui_Name', mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gui_merge_goGPSbin_OpeningFcn, ...
    'gui_OutputFcn',  @gui_merge_goGPSbin_OutputFcn, ...
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


% --- Executes just before gui_merge_goGPSbin is made visible.
function gui_merge_goGPSbin_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_merge_goGPSbin (see VARARGIN)

% Choose default command line output for gui_merge_goGPSbin
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
function varargout = gui_merge_goGPSbin_OutputFcn(hObject, eventdata, handles)  %#ok<*STOUT,*INUSD>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure


function rover_dataset_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to rover_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rover_dataset as text
%        str2double(get(hObject,'String')) returns contents of rover_dataset as a double


% --- Executes during object creation, after setting all properties.
function rover_dataset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rover_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_rover_dataset.
function browse_rover_dataset_Callback(hObject, eventdata, handles)
% hObject    handle to browse_rover_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {'*.bin','Rover dataset (*.bin)'}, ...
    'Choose rover observations dataset','../data');

if (filename ~= 0)
    pos = find(filename == '_');
    filename = filename(1:pos(end-1)-1);
    set(handles.rover_dataset,'String',fullfile(pathname, filename));
end


% --- Executes on button press in browse_master_dataset.
function browse_master_dataset_Callback(hObject, eventdata, handles)
% hObject    handle to browse_master_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {'*.bin','Rover dataset (*.bin)'}, ...
    'Choose master observations dataset','../data');

if (filename ~= 0)
    pos = find(filename == '_');
    filename = filename(1:pos(end-1)-1);
    set(handles.master_dataset,'String',fullfile(pathname, filename));
end

function data_out_folder_Callback(hObject, eventdata, handles)
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


function data_out_name_Callback(hObject, eventdata, handles)
% hObject    handle to data_out_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of data_out_name as text
%        str2double(get(hObject,'String')) returns contents of data_out_name as a double


% --- Executes during object creation, after setting all properties.
function data_out_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_out_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in merge_button.
function merge_button_Callback(hObject, eventdata, handles)
% hObject    handle to merge_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filerootR = get(handles.rover_dataset,'String');
filerootM = get(handles.master_dataset,'String');
filerootOUT = [get(handles.data_out_folder,'String') '\' get(handles.data_out_name,'String')];
filerootR(filerootR == '\') = '/';
filerootM(filerootM == '\') = '/';
filerootOUT(filerootOUT == '\') = '/';

wait_dlg = waitbar(0,'Please wait...');

goGPSbinMerge(filerootR, filerootM, filerootOUT, wait_dlg);

close(wait_dlg)

close(handles.converter_panel)

% --- Executes on button press in exit_button.
function exit_button_Callback(hObject, eventdata, handles)
% hObject    handle to exit_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.converter_panel)


function master_dataset_Callback(hObject, eventdata, handles)
% hObject    handle to master_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of master_dataset as text
%        str2double(get(hObject,'String')) returns contents of master_dataset as a double


% --- Executes during object creation, after setting all properties.
function master_dataset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to master_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
