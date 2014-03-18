function varargout = gui_seid(varargin)
% GUI_SEID MATLAB code for gui_seid.fig
%      GUI_SEID, by itself, creates a new GUI_SEID or raises the existing
%      singleton*.
%
%      H = GUI_SEID returns the handle to a new GUI_SEID or the handle to
%      the existing singleton*.
%
%      GUI_SEID('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SEID.M with the given input arguments.
%
%      GUI_SEID('Property','Value',...) creates a new GUI_SEID or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_seid_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_seid_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_seid

% Last Modified by GUIDE v2.5 18-Mar-2014 17:15:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_seid_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_seid_OutputFcn, ...
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


% --- Executes just before gui_seid is made visible.
function gui_seid_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_seid (see VARARGIN)

% Choose default command line output for gui_seid
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_seid wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_seid_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function eOutput_folder_Callback(hObject, eventdata, handles)
% hObject    handle to eOutput_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eOutput_folder as text
%        str2double(get(hObject,'String')) returns contents of eOutput_folder as a double


% --- Executes during object creation, after setting all properties.
function eOutput_folder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eOutput_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pOutput_folder.
function pOutput_folder_Callback(hObject, eventdata, handles)
% hObject    handle to pOutput_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dname = uigetdir('../data/','Choose a directory to store output data');
if (dname ~= 0)
    set(handles.eOutput_folder,'String',dname);
end

% --- Executes on button press in pBrowse_INI.
function pBrowse_INI_Callback(hObject, eventdata, handles)
% hObject    handle to pBrowse_INI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {'*.ini;','INI configuration file (*.ini)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose an INI configuration file','./settings');
if (filename ~= 0)
    set(handles.eBrowse_INI,'String',fullfile(pathname, filename));
end

% --- Executes on button press in pEdit_INI.
function pEdit_INI_Callback(hObject, eventdata, handles)
% hObject    handle to pEdit_INI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global goGUI
global goIni;
goIni = goIniReader;
goIni.setFileName(get(handles.eBrowse_INI,'String'));
    goGUI.openEditINI();


function eBrowse_INI_Callback(hObject, eventdata, handles)
% hObject    handle to eBrowse_INI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eBrowse_INI as text
%        str2double(get(hObject,'String')) returns contents of eBrowse_INI as a double


% --- Executes during object creation, after setting all properties.
function eBrowse_INI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eBrowse_INI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pExecuteSEID.
function pExecuteSEID_Callback(hObject, eventdata, handles)
% hObject    handle to pExecuteSEID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
