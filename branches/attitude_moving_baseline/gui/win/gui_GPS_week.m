function varargout = gui_GPS_week(varargin)
% GUI_GPS_WEEK M-file for gui_GPS_week.fig
%      GUI_GPS_WEEK, by itself, creates a new GUI_GPS_WEEK or raises the existing
%      singleton*.
%
%      H = GUI_GPS_WEEK returns the handle to a new GUI_GPS_WEEK or the handle to
%      the existing singleton*.
%
%      GUI_GPS_WEEK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_GPS_WEEK.M with the given input arguments.
%
%      GUI_GPS_WEEK('Property','Value',...) creates a new GUI_GPS_WEEK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_GPS_week_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_GPS_week_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_GPS_week

% Last Modified by GUIDE v2.5 28-Jun-2010 16:09:09

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
                   'gui_OpeningFcn', @gui_GPS_week_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_GPS_week_OutputFcn, ...
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


% --- Executes just before gui_GPS_week is made visible.
function gui_GPS_week_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_GPS_week (see VARARGIN)

% Choose default command line output for gui_GPS_week
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

% UIWAIT makes gui_GPS_week wait for user response (see UIRESUME)
uiwait(handles.GPS_week_panel);


% --- Outputs from this function are returned to the command line.
function varargout = gui_GPS_week_OutputFcn(hObject, eventdata, handles) 
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
