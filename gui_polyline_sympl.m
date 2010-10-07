function varargout = gui_polyline_sympl(varargin)
% GUI_POLYLINE_SYMPL M-file for gui_polyline_sympl.fig
%      GUI_POLYLINE_SYMPL, by itself, creates a new GUI_POLYLINE_SYMPL or raises the existing
%      singleton*.
%
%      H = GUI_POLYLINE_SYMPL returns the handle to a new GUI_POLYLINE_SYMPL or the handle to
%      the existing singleton*.
%
%      GUI_POLYLINE_SYMPL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_POLYLINE_SYMPL.M with the given input arguments.
%
%      GUI_POLYLINE_SYMPL('Property','Value',...) creates a new GUI_POLYLINE_SYMPL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_polyline_sympl_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_polyline_sympl_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_polyline_sympl

% Last Modified by GUIDE v2.5 07-Oct-2010 17:31:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name', mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gui_polyline_sympl_OpeningFcn, ...
    'gui_OutputFcn',  @gui_polyline_sympl_OutputFcn, ...
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


% --- Executes just before gui_polyline_sympl is made visible.
function gui_polyline_sympl_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_polyline_sympl (see VARARGIN)

% Choose default command line output for gui_polyline_sympl
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
function varargout = gui_polyline_sympl_OutputFcn(hObject, eventdata, handles)  %#ok<*STOUT,*INUSD>
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
    {'*_position.txt;*_cov.txt','Input dataset (*_position.txt and *_cov.txt)'}, ...
    'Choose input dataset','../data');

if (filename ~= 0)
    pos = find(filename == '_');
    filename = filename(1:pos(end)-1);
    set(handles.data_stream,'String',fullfile(pathname, filename));
end


% --- Executes on button press in execute_button.
function execute_button_Callback(hObject, eventdata, handles)
% hObject    handle to execute_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filerootIN = get(handles.data_stream,'String');
filerootIN(filerootIN == '\') = '/';
angle = str2double(get(handles.angle_thres,'String'));
dN1 = str2double(get(handles.disreg_begin,'String'));
dN2 = str2double(get(handles.disreg_end,'String'));

% wait_dlg = waitbar(0,'Please wait...');

%check if input data are available
d1 = dir([filerootIN '_position.txt']);
d2 = dir([filerootIN '_cov.txt']);
if ~isempty(d1) & ~isempty(d2)
    polyline(filerootIN, angle, dN1, dN2);
else
    msgbox('Both *_position.txt and *_cov.txt files are needed to run the polyline simplification algorithm.');
end

% close(wait_dlg)

close(handles.poyline_sympl_panel)

% --- Executes on button press in exit_button.
function exit_button_Callback(hObject, eventdata, handles)
% hObject    handle to exit_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.poyline_sympl_panel)


function angle_thres_Callback(hObject, eventdata, handles)
% hObject    handle to angle_thres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of angle_thres as text
%        str2double(get(hObject,'String')) returns contents of angle_thres as a double


% --- Executes during object creation, after setting all properties.
function angle_thres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle_thres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function disreg_begin_Callback(hObject, eventdata, handles)
% hObject    handle to disreg_begin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of disreg_begin as text
%        str2double(get(hObject,'String')) returns contents of disreg_begin as a double


% --- Executes during object creation, after setting all properties.
function disreg_begin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to disreg_begin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function disreg_end_Callback(hObject, eventdata, handles)
% hObject    handle to disreg_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of disreg_end as text
%        str2double(get(hObject,'String')) returns contents of disreg_end as a double


% --- Executes during object creation, after setting all properties.
function disreg_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to disreg_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
