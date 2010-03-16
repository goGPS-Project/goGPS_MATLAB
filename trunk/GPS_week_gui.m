function varargout = GPS_week_gui(varargin)
% GPS_WEEK_GUI M-file for GPS_week_gui.fig
%      GPS_WEEK_GUI, by itself, creates a new GPS_WEEK_GUI or raises the existing
%      singleton*.
%
%      H = GPS_WEEK_GUI returns the handle to a new GPS_WEEK_GUI or the handle to
%      the existing singleton*.
%
%      GPS_WEEK_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GPS_WEEK_GUI.M with the given input arguments.
%
%      GPS_WEEK_GUI('Property','Value',...) creates a new GPS_WEEK_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GPS_week_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GPS_week_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GPS_week_gui

% Last Modified by GUIDE v2.5 16-Mar-2010 18:57:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GPS_week_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @GPS_week_gui_OutputFcn, ...
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


% --- Executes just before GPS_week_gui is made visible.
function GPS_week_gui_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GPS_week_gui (see VARARGIN)

% Choose default command line output for GPS_week_gui
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

% UIWAIT makes GPS_week_gui wait for user response (see UIRESUME)
uiwait(handles.GPS_week_panel);


% --- Outputs from this function are returned to the command line.
function varargout = GPS_week_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
contents_year  = cellstr(get(handles.year,'String'));
contents_month = cellstr(get(handles.month,'String'));
contents_day   = cellstr(get(handles.day,'String'));

year  = str2num(contents_year{get(handles.year,'Value')});
month = str2num(contents_month{get(handles.month,'Value')});
day   = str2num(contents_day{get(handles.day,'Value')});

date0 = datenum(1980,1,6,0,0,0);
date1 = datenum(year,month,day,0,0,0);

GPS_week = floor((date1 - date0)/7);

varargout{1} = GPS_week;

close(handles.GPS_week_panel)


% --- Executes on selection change in year.
function year_Callback(hObject, eventdata, handles) %#ok<*DEFNU,*INUSD>
% hObject    handle to year (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns year contents as cell array
%        contents{get(hObject,'Value')} returns selected item from year


% --- Executes during object creation, after setting all properties.
function year_CreateFcn(hObject, eventdata, handles)
% hObject    handle to year (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in month.
function month_Callback(hObject, eventdata, handles)
% hObject    handle to month (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns month contents as cell array
%        contents{get(hObject,'Value')} returns selected item from month


% --- Executes during object creation, after setting all properties.
function month_CreateFcn(hObject, eventdata, handles)
% hObject    handle to month (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in day.
function day_Callback(hObject, eventdata, handles)
% hObject    handle to day (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns day contents as cell array
%        contents{get(hObject,'Value')} returns selected item from day


% --- Executes during object creation, after setting all properties.
function day_CreateFcn(hObject, eventdata, handles)
% hObject    handle to day (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ok.
function ok_Callback(hObject, eventdata, handles)
% hObject    handle to ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.GPS_week_panel);
