function varargout = gui_decode_stream(varargin)
% GUI_DECODE_STREAM M-file for gui_decode_stream.fig
%      GUI_DECODE_STREAM, by itself, creates a new GUI_DECODE_STREAM or raises the existing
%      singleton*.
%
%      H = GUI_DECODE_STREAM returns the handle to a new GUI_DECODE_STREAM or the handle to
%      the existing singleton*.
%
%      GUI_DECODE_STREAM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_DECODE_STREAM.M with the given input arguments.
%
%      GUI_DECODE_STREAM('Property','Value',...) creates a new GUI_DECODE_STREAM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_decode_stream_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_decode_stream_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_decode_stream

% Last Modified by GUIDE v2.5 28-Jun-2010 15:18:49

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.1 beta
%
% Copyright (C) 2009-2013 Mirko Reguzzoni, Eugenio Realini
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
    'gui_OpeningFcn', @gui_decode_stream_OpeningFcn, ...
    'gui_OutputFcn',  @gui_decode_stream_OutputFcn, ...
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


% --- Executes just before gui_decode_stream is made visible.
function gui_decode_stream_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_decode_stream (see VARARGIN)

% Choose default command line output for gui_decode_stream
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
if (~isempty(varargin))
    set(handles.data_stream,'String',varargin{1});
    temp = varargin{1};
    temp = temp( max(strfind(temp, '/'))+1:end );
    set(handles.data_out_name,'String',temp);
    clear temp
end

% --- Outputs from this function are returned to the command line.
function varargout = gui_decode_stream_OutputFcn(hObject, eventdata, handles)  %#ok<*STOUT,*INUSD>
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
    {'*.ubx;*.UBX;*_rover_*.bin','UBX stream data (*.ubx,*_rover_*.bin)'; ...
     '*.stq;*.STQ;*_rover_*.bin','SkyTraq stream data (*.stq,*_rover_*.bin)'; ...
     '*.rtcm;*.RTCM;*_master_*.bin','RTCM stream data (*.rtcm,*_master_*.bin)'}, ...
    'Choose stream data','../data');

if (filename ~= 0)
    if (~strcmp(filename(end-3:end),'.ubx') & ~strcmp(filename(end-3:end),'.UBX') & ...
        ~strcmp(filename(end-3:end),'.stq') & ~strcmp(filename(end-3:end),'.STQ') & ...
        ~strcmp(filename(end-4:end),'.rtcm') & ~strcmp(filename(end-4:end),'.RTCM'))
        pos = find(filename == '_');
        filename = filename(1:pos(end-1)-1);
        set(handles.data_stream,'String',fullfile(pathname, filename));
        set(handles.data_out_name,'String',filename);
    else
        if strcmp(filename(end-3:end),'.ubx')
            pos = strfind(filename,'.ubx');
        elseif strcmp(filename(end-3:end),'.UBX')
            pos = strfind(filename,'.UBX');
        elseif strcmp(filename(end-3:end),'.stq')
            pos = strfind(filename,'.stq');
        elseif strcmp(filename(end-3:end),'.STQ')
            pos = strfind(filename,'.STQ');
        elseif strcmp(filename(end-4:end),'.rtcm')
            pos = strfind(filename,'.rtcm');
        elseif strcmp(filename(end-4:end),'.RTCM')
            pos = strfind(filename,'.RTCM');
        end
        filename_out = filename(1:pos(end)-1);
        set(handles.data_stream,'String',fullfile(pathname, filename));
        set(handles.data_out_name,'String',filename_out);
    end
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

flag_no_master = 0;
if strcmp(filerootIN(end-3:end),'.ubx') | strcmp(filerootIN(end-3:end),'.UBX') | ...
   strcmp(filerootIN(end-3:end),'.stq') | strcmp(filerootIN(end-3:end),'.STQ')
    flag_no_master = 1;
end

wait_dlg = waitbar(0,'Please wait...');

%check if RINEX or goGPS data is requested
if (get(handles.output_type, 'SelectedObject') == handles.out_gogps_binary)
    streams2goGPSbin(filerootIN, filerootOUT, wait_dlg);
elseif (get(handles.output_type, 'SelectedObject') == handles.out_rinex)
    week = 0;
    if (get(handles.flag_rover_stream,'Value'))
        week = streamR2RINEX(filerootIN, [filerootOUT '_rover'], wait_dlg);
    end
    
    if (get(handles.flag_master_stream,'Value')) & (~flag_no_master)
        if (~week)
            week = gui_GPS_week;
        end
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
    %if the output folder is the default one
    if (strcmp(get(handles.data_out_folder, 'String'), '../data/data_goGPS'))
        set(handles.data_out_folder, 'String', '../data/data_RINEX');
    end
else
    set(handles.flag_rover_stream, 'Enable', 'off');
    set(handles.flag_master_stream, 'Enable', 'off');
    %if the output folder is the default one
    if (strcmp(get(handles.data_out_folder, 'String'), '../data/data_RINEX'))
        set(handles.data_out_folder, 'String', '../data/data_goGPS');
    end
end
