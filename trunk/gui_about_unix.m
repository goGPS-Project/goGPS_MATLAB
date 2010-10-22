function varargout = gui_about_unix(varargin)
% GUI_ABOUT_UNIX M-file for gui_about_unix.fig
%      GUI_ABOUT_UNIX, by itself, creates a new GUI_ABOUT_UNIX or raises the existing
%      singleton*.
%
%      H = GUI_ABOUT_UNIX returns the handle to a new GUI_ABOUT_UNIX or the handle to
%      the existing singleton*.
%
%      GUI_ABOUT_UNIX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_ABOUT_UNIX.M with the given input arguments.
%
%      GUI_ABOUT_UNIX('Property','Value',...) creates a new GUI_ABOUT_UNIX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_about_unix_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_about_unix_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_about_unix

% Last Modified by GUIDE v2.5 22-Oct-2010 10:43:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_about_unix_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_about_unix_OutputFcn, ...
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


% --- Executes just before gui_about_unix is made visible.
function gui_about_unix_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_about_unix (see VARARGIN)

% Choose default command line output for gui_about_unix
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_about_unix wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_about_unix_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)
