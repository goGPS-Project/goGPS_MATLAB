function varargout = guiEditINI(varargin)
% GUIEDITINI MATLAB code for guiEditINI.fig
%      GUIEDITINI, by itself, creates a new GUIEDITINI or raises the existing
%      singleton*.
%
%      H = GUIEDITINI returns the handle to a new GUIEDITINI or the handle to
%      the existing singleton*.
%
%      GUIEDITINI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIEDITINI.M with the given input arguments.
%
%      GUIEDITINI('Property','Value',...) creates a new GUIEDITINI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guiEditINI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guiEditINI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guiEditINI

% Last Modified by GUIDE v2.5 04-Mar-2013 18:10:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guiEditINI_OpeningFcn, ...
                   'gui_OutputFcn',  @guiEditINI_OutputFcn, ...
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


% --- Executes just before guiEditINI is made visible.
function guiEditINI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guiEditINI (see VARARGIN)

% Choose default command line output for guiEditINI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes guiEditINI wait for user response (see UIRESUME)
% uiwait(handles.wEditINI);
global goGUI
    goGUI.initEditINI(handles);

% --- Outputs from this function are returned to the command line.
function varargout = guiEditINI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function eINI_Callback(hObject, eventdata, handles)
% hObject    handle to eINI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eINI as text
%        str2double(get(hObject,'String')) returns contents of eINI as a double


% --- Executes during object creation, after setting all properties.
function eINI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eINI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sINI_Callback(hObject, eventdata, handles)
% hObject    handle to sINI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sINI as text
%        str2double(get(hObject,'String')) returns contents of sINI as a double


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
    goGUI.browseINIEditInFile();

% --- Executes on selection change in lSections.
function lSections_Callback(hObject, eventdata, handles)
% hObject    handle to lSections (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lSections contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lSections
global goGUI
goGUI.updateFieldsINI();

% --- Executes during object creation, after setting all properties.
function lSections_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lSections (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eFields_Callback(hObject, eventdata, handles)
% hObject    handle to eFields (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eFields as text
%        str2double(get(hObject,'String')) returns contents of eFields as a double


% --- Executes during object creation, after setting all properties.
function eFields_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eFields (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bSave.
function bSave_Callback(hObject, eventdata, handles)
% hObject    handle to bSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global goGUI
    goGUI.saveINI();



function sINIout_Callback(hObject, eventdata, handles)
% hObject    handle to sINIout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sINIout as text
%        str2double(get(hObject,'String')) returns contents of sINIout as a double


% --- Executes during object creation, after setting all properties.
function sINIout_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sINIout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bBrowse4Rin.
function bBrowse4Rin_Callback(hObject, eventdata, handles)
% hObject    handle to bBrowse4Rin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global goGUI
    goGUI.browse4Rin();


function sBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to sBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sBrowse as text
%        str2double(get(hObject,'String')) returns contents of sBrowse as a double


% --- Executes during object creation, after setting all properties.
function sBrowse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bBrowse4Bin.
function bBrowse4Bin_Callback(hObject, eventdata, handles)
% hObject    handle to bBrowse4Bin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global goGUI
    goGUI.browse4Bin();


% --- Executes on button press in bBrowse4Nav.
function bBrowse4Nav_Callback(hObject, eventdata, handles)
% hObject    handle to bBrowse4Nav (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global goGUI
    goGUI.browse4Nav();


% --- Executes on button press in bBrowse4Gen.
function bBrowse4Gen_Callback(hObject, eventdata, handles)
% hObject    handle to bBrowse4Gen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global goGUI
    goGUI.browse4Gen();


% --- Executes on button press in bBrowse4Path.
function bBrowse4Path_Callback(hObject, eventdata, handles)
% hObject    handle to bBrowse4Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global goGUI
    goGUI.browse4Dir();


% --- Executes on button press in bBrowse4Ref.
function bBrowse4Ref_Callback(hObject, eventdata, handles)
% hObject    handle to bBrowse4Ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global goGUI
    goGUI.browse4Ref();
    
% --- Executes on key press with focus on wEditINI and none of its controls.
function wEditINI_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to wEditINI (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
global goGUI
if length(eventdata.Modifier) == 1
    if (strcmp(eventdata.Modifier{1},'control') || strcmp(eventdata.Modifier{1},'command')) && strcmp(eventdata.Key,'s')
        goGUI.saveINI();
    end
end
