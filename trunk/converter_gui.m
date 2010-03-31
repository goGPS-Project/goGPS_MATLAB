function varargout = converter_gui(varargin)
% CONVERTER_GUI M-file for converter_gui.fig
%      CONVERTER_GUI, by itself, creates a new CONVERTER_GUI or raises the existing
%      singleton*.
%
%      H = CONVERTER_GUI returns the handle to a new CONVERTER_GUI or the handle to
%      the existing singleton*.
%
%      CONVERTER_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONVERTER_GUI.M with the given input arguments.
%
%      CONVERTER_GUI('Property','Value',...) creates a new CONVERTER_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before converter_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to converter_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help converter_gui

% Last Modified by GUIDE v2.5 16-Mar-2010 18:57:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @converter_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @converter_gui_OutputFcn, ...
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


% --- Executes just before converter_gui is made visible.
function converter_gui_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to converter_gui (see VARARGIN)

% Choose default command line output for converter_gui
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
function varargout = converter_gui_OutputFcn(hObject, eventdata, handles)  %#ok<*STOUT,*INUSD>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure


function data_stream_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to data_stream (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of data_stream as text
%        str2double(get(hObject,'String')) returns contents of data_stream as a double


% --- Executes during object creation, after setting all properties.
function data_stream_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_stream (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_data_stream.
function browse_data_stream_Callback(hObject, eventdata, handles)
% hObject    handle to browse_data_stream (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {'*.bin','Stream data (*.bin)'}, ...
    'Choose stream data','../data');

if (filename ~= 0)
    pos = find(filename == '_');
    filename = filename(1:pos(end-1)-1);
    set(handles.data_stream,'String',fullfile(pathname, filename));
    set(handles.data_out_name,'String',filename);
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


% --- Executes on button press in flag_rover_stream.
function flag_rover_stream_Callback(hObject, eventdata, handles)
% hObject    handle to flag_rover_stream (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of flag_rover_stream


% --- Executes on button press in flag_master_stream.
function flag_master_stream_Callback(hObject, eventdata, handles)
% hObject    handle to flag_master_stream (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of flag_master_stream


% --- Executes on button press in convert_button.
function convert_button_Callback(hObject, eventdata, handles)
% hObject    handle to convert_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filerootIN = get(handles.data_stream,'String');
filerootOUT = [get(handles.data_out_folder,'String') '\' get(handles.data_out_name,'String')];
filerootIN(filerootIN == '\') = '/';
filerootOUT(filerootOUT == '\') = '/';

wait_dlg = waitbar(0,'Please wait...');

%check if RINEX or goGPS data is requested
if (get(handles.output_type, 'SelectedObject') == handles.out_gogps_binary)
    streams2goGPSbin(filerootIN, filerootOUT, wait_dlg);
elseif (get(handles.output_type, 'SelectedObject') == handles.out_rinex)
    week = 0;
    if (get(handles.flag_rover_stream,'Value'))
        week = streamR2RINEX(filerootIN, [filerootOUT '_rover'], wait_dlg);
    end
    if (~week)
        week = GPS_week_gui;
    end
    if (week) & (get(handles.flag_master_stream,'Value'))
        streamM2RINEX(filerootIN, [filerootOUT '_master'], week, wait_dlg);
    end
end

close(wait_dlg)

close(handles.converter_panel)

% --- Executes on button press in exit_button.
function exit_button_Callback(hObject, eventdata, handles)
% hObject    handle to exit_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.converter_panel)


% --- Executes when selected object is changed in output_type.
function output_type_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in output_type 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if (hObject == handles.out_rinex)
    set(handles.flag_rover_stream, 'Enable', 'on');
    set(handles.flag_master_stream, 'Enable', 'on');
else
    set(handles.flag_rover_stream, 'Enable', 'off');
    set(handles.flag_master_stream, 'Enable', 'off');
end
