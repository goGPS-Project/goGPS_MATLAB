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

% Last Modified by GUIDE v2.5 19-Mar-2014 14:46:29

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
global goGUI goIni oldIniFile
oldIniFile = goIni.getFileName;
goIni.setFileName(get(handles.eBrowse_INI,'String'));
goGUI.setElVal(goGUI.idUI.sINI, get(handles.eBrowse_INI,'String'));
goGUI.forceINIupdate();

% Choose default command line output for gui_seid
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

pCheckFiles_Callback(handles.pCheckFiles, [], handles);

% UIWAIT makes gui_seid wait for user response (see UIRESUME)
uiwait(handles.seid_panel);

% --- Outputs from this function are returned to the command line.
function varargout = gui_seid_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
global goIni oldIniFile goGUI

goIni.setFileName(oldIniFile);
goGUI.setElVal(goGUI.idUI.sINI, oldIniFile);
goGUI.forceINIupdate();


function eOutput_folder_Callback(hObject, eventdata, handles)
% hObject    handle to eOutput_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eOutput_folder as text
%        str2double(get(hObject,'String')) returns contents of eOutput_folder as a double
pCheckFiles_Callback(handles.pCheckFiles, [], handles);

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
resetLEDs(handles);
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
resetLEDs(handles);
goGUI.openEditINI();


function eBrowse_INI_Callback(hObject, eventdata, handles)
% hObject    handle to eBrowse_INI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eBrowse_INI as text
%        str2double(get(hObject,'String')) returns contents of eBrowse_INI as a double
pCheckFiles_Callback(handles.pCheckFiles, [], handles);

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
global goIni oldIniFile goGUI
data_path_src = goIni.getData('Receivers','data_path');
file_name_src = goIni.getData('Receivers','file_name');

data_path_tgt = goIni.getData('Master','data_path');
file_name_tgt = goIni.getData('Master','file_name');

data_path_nav = goIni.getData('Navigational','data_path');
file_name_nav = goIni.getData('Navigational','file_name');

data_path_pcv = goIni.getData('Antenna PCV','data_path');
file_name_pcv = goIni.getData('Antenna PCV','file_name');

goIni.setFileName(oldIniFile);
goGUI.setElVal(goGUI.idUI.sINI, oldIniFile);
goGUI.forceINIupdate();

close(handles.seid_panel)


% --- Executes on button press in pCheckFiles.
function pCheckFiles_Callback(hObject, eventdata, handles)
% hObject    handle to pCheckFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global goGUI goIni

ready = 0;

% 1. INI file -------------------------------------
ini_file = get(handles.eBrowse_INI,'String');
if (isempty(ini_file))
    set(handles.fINI,'ForegroundColor',goGUI.red)
else
    % Check the presence of all the files
    if exist(ini_file,'file')
        set(handles.fINI,'ForegroundColor',goGUI.green);
        ready = ready + 1;
    else
        set(handles.fINI,'ForegroundColor',goGUI.yellow);
    end
end

% 2. Source files --------------------------------------
nR = goIni.getData('Receivers','nRec');
data_path = goIni.getData('Receivers','data_path');
file_name = goIni.getData('Receivers','file_name');

if (isempty(nR))
    if iscell(file_name)
        nR = length(file_name);
    else
        nR = 1;
    end
    goIni.addKey('Receivers','nRec',nR);
end
set(handles.tNumSources,'String',['x ' num2str(nR)]);

if (isempty(data_path))
    data_path = '';
end
if (isempty(file_name))
    set(handles.fSrc,'ForegroundColor',goGUI.red);
else
    % If I have more than one receiver
    if iscell(file_name)
        % The number of receiver is = to the number of files?
        if nR ~= length((file_name))   % Declared number of file ~= number of files
            set(handles.fSrc,'ForegroundColor',goGUI.yellow);
        else
            % Check the presence of all the files
            fileOk = true;
            for r = 1:nR
                if ~exist([data_path file_name{r}],'file')
                    fileOk = false;
                end
            end
            if fileOk
                set(handles.fSrc,'ForegroundColor',goGUI.green);
                ready = ready + 1;
            else
                set(handles.fSrc,'ForegroundColor',goGUI.yellow);
            end
        end
    else
        % The number of receiver is = to the number of files?
        if nR > 1   % Declared number of file ~= number of files
            set(handles.fSrc,'ForegroundColor',goGUI.yellow);
        else
            % Check the presence of all the files
            if exist([data_path file_name],'file')
                set(handles.fSrc,'ForegroundColor',goGUI.green);
                ready = ready + 1;
            else
                set(handles.fSrc,'ForegroundColor',goGUI.yellow);
            end
        end
    end
end

% 3. Target file -----------------------------------------
data_path = goIni.getData('Master','data_path');
file_name = goIni.getData('Master','file_name');
if (isempty(data_path))
    data_path = '';
end
if (isempty(file_name))
    set(handles.fTgt,'ForegroundColor',goGUI.red);
else
    % Check the presence of all the files
    if exist([data_path file_name],'file')
        set(handles.fTgt,'ForegroundColor',goGUI.green);
        ready = ready + 1;
    else
        set(handles.fTgt,'ForegroundColor',goGUI.yellow);
    end
end

% 4. Navigation file -------------------------------------
data_path = goIni.getData('Navigational','data_path');
file_name = goIni.getData('Navigational','file_name');
if (isempty(data_path))
    data_path = '';
end
if (isempty(file_name))
    set(handles.fNav,'ForegroundColor',goGUI.red);
else
    % Check the presence of all the files
    if exist([data_path file_name],'file')
        set(handles.fNav,'ForegroundColor',goGUI.green);
        ready = ready + 1;
    else
        set(handles.fNav,'ForegroundColor',goGUI.yellow);
    end
end

% 5. Antenna PCO/PCV file -------------------------------------
data_path = goIni.getData('Antenna PCV','data_path');
file_name = goIni.getData('Antenna PCV','file_name');
if (isempty(data_path))
    data_path = '';
end
if (isempty(file_name))
    set(handles.fPcv,'ForegroundColor',goGUI.red)
else
    % Check the presence of all the files
    if exist([data_path file_name],'file')
        set(handles.fPcv,'ForegroundColor',goGUI.green);
        ready = ready + 1;
    else
        set(handles.fNav,'ForegroundColor',goGUI.yellow);
    end
end

% 6. Output directory -------------------------------------
out_path = get(handles.eOutput_folder,'String');
if (isempty(out_path))
    set(handles.fOut,'ForegroundColor',goGUI.red)
else
    % Check the presence of all the files
    if exist(out_path,'dir')
        set(handles.fOut,'ForegroundColor',goGUI.green);
        ready = ready + 1;
    else
        set(handles.fOut,'ForegroundColor',goGUI.yellow);
    end
end

if (ready == 6)
    set(handles.pExecuteSEID,'Enable','on')
else
    set(handles.pExecuteSEID,'Enable','off')
end

function resetLEDs(handles)
set(handles.pExecuteSEID,'Enable','off')
set(handles.fINI,'ForegroundColor','black')
set(handles.fSrc,'ForegroundColor','black')
set(handles.fTgt,'ForegroundColor','black')
set(handles.fNav,'ForegroundColor','black')
set(handles.fPcv,'ForegroundColor','black')
set(handles.fOut,'ForegroundColor','black')
