function varargout = gui_goGPS_unix(varargin)
% GUI_GOGPS_UNIX M-file for gui_goGPS_unix.fig
%      GUI_GOGPS_UNIX, by itself, creates a new GUI_GOGPS_UNIX or raises the existing
%      singleton*.
%
%      H = GUI_GOGPS_UNIX returns the handle to a new GUI_GOGPS_UNIX or the handle to
%      the existing singleton*.
%
%      GUI_GOGPS_UNIX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_GOGPS_UNIX.M with the given input arguments.
%
%      GUI_GOGPS_UNIX('Property','Value',...) creates a new GUI_GOGPS_UNIX or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_goGPS_unix_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_goGPS_unix_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_goGPS_unix

% Last Modified by GUIDE v2.5 14-May-2011 12:06:31

%----------------------------------------------------------------------------------------------
%                           goGPS v0.2.0 beta
%
% Copyright (C) 2009-2011 Mirko Reguzzoni, Eugenio Realini
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
    'gui_OpeningFcn', @gui_goGPS_unix_OpeningFcn, ...
    'gui_OutputFcn',  @gui_goGPS_unix_OutputFcn, ...
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

% --- Executes just before gui_goGPS_unix is made visible.
function gui_goGPS_unix_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_goGPS_unix (see VARARGIN)

% Choose default command line output for gui_goGPS_unix
handles.output = hObject;

set(hObject,'CloseRequestFcn',@closeGUI);

% Update handles structure
guidata(hObject, handles);

if exist('../data/settings/last_settings.mat','file')
    loadState(handles, '../data/settings/last_settings.mat');
else
    loadState(handles, '../data/settings/default_settings.mat');
end

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

% UIWAIT makes gui_goGPS_unix wait for user response (see UIRESUME)
uiwait(handles.main_panel);


% --- Outputs from this function are returned to the command line.
function varargout = gui_goGPS_unix_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(~isstruct(handles))
    varargout = cell(22,1);
    return
end
mode = select_mode(handles);
mode_vinc = get(handles.constraint,'Value');
if (get(handles.file_type, 'SelectedObject') == handles.rinex_files)
    mode_data = 0;
else %goGPS data
    mode_data = 1;
end
contents_dyn_mod = cellstr(get(handles.dyn_mod,'String'));
if (strcmp(contents_dyn_mod{get(handles.dyn_mod,'Value')},'Variable') | get(handles.stopGOstop,'Value'))
    flag_var_dyn_model = 1;
else
    flag_var_dyn_model = 0;
end
mode_ref = get(handles.ref_path,'Value');
flag_ms_pos = get(handles.master_pos,'Value');
flag_ms = get(handles.plot_master,'Value');
flag_ge = get(handles.google_earth,'Value');
flag_cov = get(handles.err_ellipse,'Value');
flag_NTRIP = get(handles.use_ntrip,'Value');
flag_amb = get(handles.plot_amb,'Value');
flag_skyplot = get(handles.no_skyplot_snr,'Value');
flag_plotproc = get(handles.plotproc,'Value');
flag_stopGOstop = get(handles.stopGOstop,'Value');
filerootIN = get(handles.gogps_data_input,'String');
filerootOUT = [get(handles.gogps_data_output,'String') '\' get(handles.gogps_data_output_prefix,'String')];
filerootIN(filerootIN == '\') = '/';
filerootOUT(filerootOUT == '\') = '/';
i = 1;
j = length(filerootOUT);
while (~isempty(dir([filerootOUT '_????_rover.bin'])) | ...
       ~isempty(dir([filerootOUT '_master*.bin'])) | ...
       ~isempty(dir([filerootOUT '_????_obs*.bin'])) | ...
       ~isempty(dir([filerootOUT '_????_eph*.bin'])) | ...
       ~isempty(dir([filerootOUT '_????_dyn*.bin'])) | ...
       ~isempty(dir([filerootOUT '_sat*.bin'])) | ...
       ~isempty(dir([filerootOUT '_kal*.bin'])) | ...
       ~isempty(dir([filerootOUT '_dt*.bin'])) | ...
       ~isempty(dir([filerootOUT '_conf*.bin'])) | ...
       ~isempty(dir([filerootOUT '_dop*.bin'])) | ...
       ~isempty(dir([filerootOUT '_ECEF*.txt'])) | ...
       ~isempty(dir([filerootOUT '_geod*.txt'])) | ...
       ~isempty(dir([filerootOUT '_plan*.txt'])) | ...
       ~isempty(dir([filerootOUT '_????_NMEA*.txt'])) | ...
       ~isempty(dir([filerootOUT '.kml'])) )

   filerootOUT(j+1:j+3) = ['_' num2str(i,'%02d')];
   i = i + 1;
end
filename_R_obs = get(handles.RINEX_rover_obs,'String');
filename_M_obs = get(handles.RINEX_master_obs,'String');
filename_nav = get(handles.RINEX_nav,'String');
filename_ref = get(handles.ref_path_input,'String');

contents = cellstr(get(handles.crs,'String'));
if (strcmp(contents{get(handles.crs,'Value')},'ECEF (X,Y,Z)'))
    XM = str2double(get(handles.master_X,'String'));
    YM = str2double(get(handles.master_Y,'String'));
    ZM = str2double(get(handles.master_Z,'String'));
else
    latM = str2double(get(handles.master_lat,'String'));
    lonM = str2double(get(handles.master_lon,'String'));
    hM = str2double(get(handles.master_h,'String'));
    [XM, YM, ZM] = geod2cart (latM*pi/180, lonM*pi/180, hM, 6378137, 1/298.257222101);
end
pos_M_man = [XM; YM; ZM];

contents = cellstr(get(handles.num_receivers,'String'));
num_rec = str2double(contents{get(handles.num_receivers,'Value')});

if num_rec >= 1
    contentsProt = cellstr(get(handles.protocol_select_0,'String'));
    if (strcmp(contentsProt{get(handles.protocol_select_0,'Value')},'UBX (u-blox)'))
        protocol_idx(1) = 0;
    elseif (strcmp(contentsProt{get(handles.protocol_select_0,'Value')},'iTalk (Fastrax)'))
        protocol_idx(1) = 1;
    elseif (strcmp(contentsProt{get(handles.protocol_select_0,'Value')},'SkyTraq'))
        protocol_idx(1) = 2;
    end

    if num_rec >= 2
        contentsProt = cellstr(get(handles.protocol_select_1,'String'));
        if (strcmp(contentsProt{get(handles.protocol_select_1,'Value')},'UBX (u-blox)'))
            protocol_idx(2) = 0;
        elseif (strcmp(contentsProt{get(handles.protocol_select_1,'Value')},'iTalk (Fastrax)'))
            protocol_idx(2) = 1;
        elseif (strcmp(contentsProt{get(handles.protocol_select_1,'Value')},'SkyTraq'))
            protocol_idx(2) = 2;
        end

        if num_rec >= 3
            contentsProt = cellstr(get(handles.protocol_select_2,'String'));
            if (strcmp(contentsProt{get(handles.protocol_select_2,'Value')},'UBX (u-blox)'))
                protocol_idx(3) = 0;
            elseif (strcmp(contentsProt{get(handles.protocol_select_2,'Value')},'iTalk (Fastrax)'))
                protocol_idx(3) = 1;
            elseif (strcmp(contentsProt{get(handles.protocol_select_2,'Value')},'SkyTraq'))
                protocol_idx(3) = 2;
            end
            
            if num_rec >= 4
                contentsProt = cellstr(get(handles.protocol_select_3,'String'));
                if (strcmp(contentsProt{get(handles.protocol_select_3,'Value')},'UBX (u-blox)'))
                    protocol_idx(4) = 0;
                elseif (strcmp(contentsProt{get(handles.protocol_select_3,'Value')},'iTalk (Fastrax)'))
                    protocol_idx(4) = 1;
                elseif (strcmp(contentsProt{get(handles.protocol_select_3,'Value')},'SkyTraq'))
                    protocol_idx(4) = 2;
                end
            end
        end
    end
end

varargout{1} = mode;
varargout{2} = mode_vinc;
varargout{3} = mode_data;
varargout{4} = mode_ref;
varargout{5} = flag_ms_pos;
varargout{6} = flag_ms;
varargout{7} = flag_ge;
varargout{8} = flag_cov;
varargout{9} = flag_NTRIP;
varargout{10} = flag_amb;
varargout{11} = flag_skyplot;
varargout{12} = flag_plotproc;
varargout{13} = flag_var_dyn_model;
varargout{14} = flag_stopGOstop;
varargout{15} = filerootIN;
varargout{16} = filerootOUT;
varargout{17} = filename_R_obs;
varargout{18} = filename_M_obs;
varargout{19} = filename_nav;
varargout{20} = filename_ref;
varargout{21} = pos_M_man;
varargout{22} = protocol_idx;

global sigmaq0 sigmaq_vE sigmaq_vN sigmaq_vU sigmaq_vel
global sigmaq_cod1 sigmaq_cod2 sigmaq_ph sigmaq0_N sigmaq_dtm
global min_nsat cutoff snr_threshold cs_threshold weights snr_a snr_0 snr_1 snr_A order o1 o2 o3
global h_antenna
global tile_header tile_georef dtm_dir
global master_ip master_port ntrip_user ntrip_pw ntrip_mountpoint
global nmea_init
global flag_doppler_cs
global COMportR

contents = cellstr(get(handles.com_select_0,'String'));
COMportR0 = contents{get(handles.com_select_0,'Value')};
contents = cellstr(get(handles.com_select_1,'String'));
COMportR1 = contents{get(handles.com_select_1,'Value')};
contents = cellstr(get(handles.com_select_2,'String'));
COMportR2 = contents{get(handles.com_select_2,'Value')};
contents = cellstr(get(handles.com_select_3,'String'));
COMportR3 = contents{get(handles.com_select_3,'Value')};

if num_rec >= 1
    COMportR{1,1} = COMportR0;
    if num_rec >= 2
        COMportR{2,1} = COMportR1;
        if num_rec >= 3
            COMportR{3,1} = COMportR2;
            if num_rec >= 4
                COMportR{4,1} = COMportR3;
            end
        end
    end
end

flag_doppler_cs = get(handles.flag_doppler,'Value');
sigmaq0 = str2double(get(handles.std_init,'String'))^2;
sigmaq_vE = str2double(get(handles.std_X,'String'))^2;
sigmaq_vN = str2double(get(handles.std_Y,'String'))^2;
sigmaq_vU = str2double(get(handles.std_Z,'String'))^2;
sigmaq_vel = str2double(get(handles.std_vel,'String'))^2;
sigmaq_cod1 = str2double(get(handles.std_code,'String'))^2;
sigmaq_cod2 = 0.16;
if (get(handles.toggle_std_phase,'Value'))
    sigmaq_ph = str2double(get(handles.std_phase,'String'))^2;
else
    sigmaq_ph = 1e30;
end
sigmaq0_N = 100;
if (get(handles.toggle_std_dtm,'Value'))
    sigmaq_dtm = str2double(get(handles.std_dtm,'String'))^2;
else
    sigmaq_dtm = 1e30;
end
min_nsat = str2double(get(handles.min_sat,'String'));
if (mode == 2)
    disp('Minimum number of satellites is forced to 4 (for stand-alone positioning)');
    min_nsat = 4;
end
cutoff = str2double(get(handles.cut_off,'String'));
snr_threshold = str2double(get(handles.snr_thres,'String'));
cs_threshold = str2double(get(handles.cs_thresh,'String'));
if (get(handles.weight_select, 'SelectedObject') == handles.weight_0)
    weights = 0;
elseif (get(handles.weight_select, 'SelectedObject') == handles.weight_1)
    weights = 1;
elseif (get(handles.weight_select, 'SelectedObject') == handles.weight_2)
    weights = 2;
elseif (get(handles.weight_select, 'SelectedObject') == handles.weight_3)
    weights = 3;
end
snr_a = 30;
snr_0 = 10;
snr_1 = 50;
snr_A = 30;
amb_select_Callback(handles.amb_select, eventdata, handles);
dyn_mod_Callback(handles.dyn_mod, eventdata, handles);
o1 = order;
o2 = order*2;
o3 = order*3;
h_antenna = str2double(get(handles.antenna_h,'String'));
dtm_dir = get(handles.dtm_path,'String');
try
    load([dtm_dir '/tiles/tile_header'], 'tile_header');
    load([dtm_dir '/tiles/tile_georef'], 'tile_georef');
catch
    tile_header.nrows = 0;
    tile_header.ncols = 0;
    tile_header.cellsize = 0;
    tile_header.nodata = 0;
    tile_georef = zeros(1,1,4);
end
master_ip = get(handles.IP_address,'String');
master_port = str2double(get(handles.port,'String'));
ntrip_user = get(handles.username,'String');
ntrip_pw = get(handles.password,'Userdata');
ntrip_mountpoint = get(handles.mountpoint,'String');
phiApp = str2double(get(handles.approx_lat,'String'));
lamApp = str2double(get(handles.approx_lon,'String'));
hApp = str2double(get(handles.approx_h,'String'));
[XApp,YApp,ZApp] = geod2cart (phiApp*pi/180, lamApp*pi/180, hApp, 6378137, 1/298.257222101);
if ~isnan(XApp) & ~isnan(YApp) & ~isnan(ZApp)
    nmea_init = NMEA_GGA_gen([XApp YApp ZApp],10);
else
    nmea_init = '';
end
    
%close main panel
delete(gcf)

% --------------------------------------------------------------------
function closeGUI(src,evnt) %#ok<*INUSD>

selection = questdlg('Do you want to quit goGPS?',...
    'Close Request Function',...
    'Yes','No','Yes');
switch selection,
    case 'Yes',
        delete(gcf)
    case 'No'
        return
end

% --------------------------------------------------------------------
function saveState(handles,filename)
state.master_pos = get(handles.master_pos, 'Value');
state.constraint = get(handles.constraint, 'Value');
state.ref_path = get(handles.ref_path, 'Value');
state.plot_master = get(handles.plot_master, 'Value');
state.google_earth = get(handles.google_earth, 'Value');
state.err_ellipse = get(handles.err_ellipse, 'Value');
state.use_ntrip = get(handles.use_ntrip, 'Value');
state.plot_amb = get(handles.plot_amb, 'Value');
state.no_skyplot_snr = get(handles.no_skyplot_snr, 'Value');
state.plotproc = get(handles.plotproc, 'Value');
state.RINEX_rover_obs = get(handles.RINEX_rover_obs,'String');
state.RINEX_master_obs = get(handles.RINEX_master_obs,'String');
state.RINEX_nav = get(handles.RINEX_nav,'String');
state.gogps_data_input = get(handles.gogps_data_input,'String');
state.gogps_data_output = get(handles.gogps_data_output,'String');
state.gogps_data_output_prefix = get(handles.gogps_data_output_prefix,'String');
state.dtm_path = get(handles.dtm_path,'String');
state.ref_path_input = get(handles.ref_path_input,'String');
state.master_X = get(handles.master_X,'String');
state.master_Y = get(handles.master_Y,'String');
state.master_Z = get(handles.master_Z,'String');
state.master_lat = get(handles.master_lat,'String');
state.master_lon = get(handles.master_lon,'String');
state.master_h = get(handles.master_h,'String');
state.crs = get(handles.crs,'Value');
state.mode = get(handles.mode,'Value');
state.nav_mon = get(handles.nav_mon,'Value');
state.kalman_ls = get(handles.kalman_ls,'Value');
state.code_dd_sa = get(handles.code_dd_sa,'Value');
state.rinex_files = get(handles.rinex_files,'Value');
state.gogps_data = get(handles.gogps_data,'Value');
state.std_X = get(handles.std_X,'String');
state.std_Y = get(handles.std_Y,'String');
state.std_Z = get(handles.std_Z,'String');
state.std_code = get(handles.std_code,'String');
state.std_phase = get(handles.std_phase,'String');
state.std_dtm = get(handles.std_dtm,'String');
state.toggle_std_phase = get(handles.toggle_std_phase,'Value');
state.toggle_std_dtm = get(handles.toggle_std_dtm,'Value');
state.std_init = get(handles.std_init,'String');
state.std_vel = get(handles.std_vel,'String');
state.cs_thresh = get(handles.cs_thresh,'String');
state.flag_doppler = get(handles.flag_doppler,'Value');
state.amb_select = get(handles.amb_select,'Value');
state.cut_off = get(handles.cut_off,'String');
state.snr_thres = get(handles.snr_thres,'String');
state.antenna_h = get(handles.antenna_h,'String');
state.min_sat = get(handles.min_sat,'String');
state.dyn_mod = get(handles.dyn_mod,'Value');
state.weight_0 = get(handles.weight_0,'Value');
state.weight_1 = get(handles.weight_1,'Value');
state.weight_2 = get(handles.weight_2,'Value');
state.weight_3 = get(handles.weight_3,'Value');
contents = cellstr(get(handles.com_select_0,'String'));
state.com_select_0 = contents{get(handles.com_select_0,'Value')};
contents = cellstr(get(handles.com_select_1,'String'));
state.com_select_1 = contents{get(handles.com_select_1,'Value')};
contents = cellstr(get(handles.com_select_2,'String'));
state.com_select_2 = contents{get(handles.com_select_2,'Value')};
contents = cellstr(get(handles.com_select_3,'String'));
state.com_select_3 = contents{get(handles.com_select_3,'Value')};
state.protocol_select_0 = get(handles.protocol_select_0,'Value');
state.protocol_select_1 = get(handles.protocol_select_1,'Value');
state.protocol_select_2 = get(handles.protocol_select_2,'Value');
state.protocol_select_3 = get(handles.protocol_select_3,'Value');
state.num_receivers = get(handles.num_receivers,'Value');
state.IP_address = get(handles.IP_address,'String');
state.port = get(handles.port,'String');
state.mountpoint = get(handles.mountpoint,'String');
state.username = get(handles.username,'String');
state.password = get(handles.password,'Userdata');
state.approx_lat = get(handles.approx_lat,'String');
state.approx_lon = get(handles.approx_lon,'String');
state.approx_h = get(handles.approx_h,'String');
state.stopGOstop = get(handles.stopGOstop,'Value');

save(filename, 'state');

% --------------------------------------------------------------------
function loadState(handles,filename)

load(filename);

set(handles.master_pos, 'Value', state.master_pos);
set(handles.constraint, 'Value', state.constraint);
set(handles.ref_path, 'Value', state.ref_path);
set(handles.plot_master, 'Value', state.plot_master);
set(handles.google_earth, 'Value', state.google_earth);
set(handles.err_ellipse, 'Value', state.err_ellipse);
set(handles.use_ntrip, 'Value', state.use_ntrip);
set(handles.plot_amb, 'Value', state.plot_amb);
set(handles.no_skyplot_snr, 'Value', state.no_skyplot_snr);
set(handles.plotproc, 'Value', state.plotproc);
set(handles.RINEX_rover_obs,'String', state.RINEX_rover_obs);
set(handles.RINEX_master_obs,'String', state.RINEX_master_obs);
set(handles.RINEX_nav,'String', state.RINEX_nav);
set(handles.gogps_data_input,'String', state.gogps_data_input);
set(handles.gogps_data_output,'String', state.gogps_data_output);
set(handles.gogps_data_output_prefix,'String', state.gogps_data_output_prefix);
set(handles.ref_path_input,'String', state.ref_path_input);
set(handles.dtm_path,'String', state.dtm_path);
set(handles.master_X,'String', state.master_X);
set(handles.master_Y,'String', state.master_Y);
set(handles.master_Z,'String', state.master_Z);
set(handles.master_lat,'String', state.master_lat);
set(handles.master_lon,'String', state.master_lon);
set(handles.master_h,'String', state.master_h);
set(handles.crs,'Value', state.crs);
set(handles.mode,'Value', state.mode);
set(handles.nav_mon,'Value', state.nav_mon);
set(handles.kalman_ls,'Value', state.kalman_ls);
set(handles.code_dd_sa,'Value', state.code_dd_sa);
set(handles.rinex_files,'Value', state.rinex_files);
set(handles.gogps_data,'Value', state.gogps_data);
set(handles.std_X,'String', state.std_X);
set(handles.std_Y,'String', state.std_Y);
set(handles.std_Z,'String', state.std_Z);
set(handles.std_code,'String', state.std_code);
set(handles.std_phase,'String', state.std_phase);
set(handles.std_dtm,'String', state.std_dtm);
set(handles.toggle_std_phase,'Value', state.toggle_std_phase);
set(handles.toggle_std_dtm,'Value', state.toggle_std_dtm);
set(handles.std_init,'String', state.std_init);
set(handles.std_vel,'String', state.std_vel);
set(handles.cs_thresh,'String', state.cs_thresh);
set(handles.flag_doppler,'Value', state.flag_doppler);
set(handles.amb_select,'Value', state.amb_select);
set(handles.cut_off,'String', state.cut_off);
set(handles.snr_thres,'String', state.snr_thres);
set(handles.antenna_h,'String', state.antenna_h);
set(handles.min_sat,'String', state.min_sat);
set(handles.dyn_mod,'Value', state.dyn_mod);
set(handles.weight_0,'Value', state.weight_0);
set(handles.weight_1,'Value', state.weight_1);
set(handles.weight_2,'Value', state.weight_2);
set(handles.weight_3,'Value', state.weight_3);
contents = get(handles.com_select_0,'String');
select_0 = 1; select_1 = 1; select_2 = 1; select_3 = 1;
for i = 1 : numel(contents)
    if (strcmp(contents{i},state.com_select_0))
        select_0 = i;
    elseif (strcmp(contents{i},state.com_select_1))
        select_1 = i;
    elseif (strcmp(contents{i},state.com_select_2))
        select_2 = i;
    elseif (strcmp(contents{i},state.com_select_3))
        select_3 = i;
    end
end
set(handles.com_select_0,'Value', select_0);
set(handles.com_select_1,'Value', select_1);
set(handles.com_select_2,'Value', select_2);
set(handles.com_select_3,'Value', select_3);
set(handles.protocol_select_0,'Value', state.protocol_select_0);
set(handles.protocol_select_1,'Value', state.protocol_select_1);
set(handles.protocol_select_2,'Value', state.protocol_select_2);
set(handles.protocol_select_3,'Value', state.protocol_select_3);
set(handles.num_receivers,'Value', state.num_receivers);
set(handles.IP_address,'String', state.IP_address);
set(handles.port,'String', state.port);
set(handles.mountpoint,'String', state.mountpoint);
set(handles.username,'String', state.username);
set(handles.password,'Userdata', state.password);
show_password_Callback(handles.show_password, [], handles);
set(handles.approx_lat,'String', state.approx_lat);
set(handles.approx_lon,'String', state.approx_lon);
set(handles.approx_h,'String', state.approx_h);
set(handles.stopGOstop,'Value', state.stopGOstop);

plot_amb_Callback(handles.plot_amb, [], handles);
constraint_Callback(handles.constraint, [], handles);
plotproc_Callback(handles.constraint, [], handles);
master_pos_Callback(handles.master_pos, [], handles);
kalman_ls_Callback(handles.kalman_ls, [], handles);
dyn_mod_Callback(handles.dyn_mod, [], handles);
mode_Callback(handles.mode, [], handles);
num_receivers_Callback(handles.num_receivers, [], handles);
stopGOstop_Callback(handles.stopGOstop, [], handles);
if(get(handles.file_type, 'SelectedObject') == handles.rinex_files);
    file_type_SelectionChangeFcn(handles.rinex_files, [], handles);
else
    file_type_SelectionChangeFcn(handles.gogps_data, [], handles);
end


% --- Executes on selection change in mode.
function mode_Callback(hObject, eventdata, handles)
% hObject    handle to mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns mode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mode
contents = cellstr(get(hObject,'String'));
if (strcmp(contents{get(hObject,'Value')},'Real-time'))
    try
        instrhwinfo;
    catch
        warndlg('Instrument Control Toolbox is needed to run goGPS in real-time mode.', 'Warning');
        set(handles.mode, 'Value', 2);
        mode_Callback(hObject, eventdata, handles);
        return
    end
    set(handles.nav_mon, 'Enable', 'on');
    set(handles.kalman_ls, 'Enable', 'off');
    set(handles.kalman_ls, 'Value', 1);
    set(handles.code_dd_sa, 'Enable', 'off');
    set(handles.code_dd_sa, 'Value', 4);
    set(handles.rinex_files, 'Enable', 'off');
    set(handles.gogps_data, 'Enable', 'off');
    
    set(handles.plot_amb, 'Enable', 'off');
    set(handles.no_skyplot_snr, 'Enable', 'on');
    set(handles.plotproc, 'Enable', 'on');
    set(handles.flag_doppler, 'Enable', 'off');
    plotproc_Callback(handles.plotproc, eventdata, handles);

    nav_mon_Callback(handles.nav_mon, eventdata, handles);

    %disable file input fields
    set(handles.RINEX_rover_obs, 'Enable', 'off');
    set(handles.RINEX_master_obs, 'Enable', 'off');
    set(handles.RINEX_nav, 'Enable', 'off');
    set(handles.browse_rover_obs, 'Enable', 'off');
    set(handles.browse_master_obs, 'Enable', 'off');
    set(handles.browse_nav, 'Enable', 'off');
    set(handles.text_RINEX_rover_obs, 'Enable', 'off');
    set(handles.text_RINEX_master_obs, 'Enable', 'off');
    set(handles.text_RINEX_nav, 'Enable', 'off');
    set(handles.gogps_data_input, 'Enable', 'off');
    set(handles.browse_gogps_input, 'Enable', 'off');
    set(handles.text_gogps_input, 'Enable', 'off');

else
    set(handles.nav_mon, 'Enable', 'off');
    set(handles.nav_mon, 'Value', 1);

    nav_mon_Callback(handles.nav_mon, eventdata, handles);

    set(handles.kalman_ls, 'Enable', 'on');
    set(handles.code_dd_sa, 'Enable', 'on');
    set(handles.rinex_files, 'Enable', 'on');
    set(handles.gogps_data, 'Enable', 'on');
    set(handles.text_num_receivers, 'Enable', 'off');
    set(handles.num_receivers, 'Enable', 'off');
    set(handles.com_select_0, 'Enable', 'off');
    set(handles.com_select_1, 'Enable', 'off');
    set(handles.com_select_2, 'Enable', 'off');
    set(handles.com_select_3, 'Enable', 'off');
    set(handles.protocol_select_0, 'Enable', 'off');
    set(handles.protocol_select_1, 'Enable', 'off');
    set(handles.protocol_select_2, 'Enable', 'off');
    set(handles.protocol_select_3, 'Enable', 'off');    
    set(handles.use_ntrip, 'Enable', 'off');
    contents = cellstr(get(handles.code_dd_sa,'String'));
    if (get(handles.plotproc,'Value') & (strcmp(contents{get(handles.code_dd_sa,'Value')}, ...
        'Code and phase double difference') | strcmp(contents{get(handles.code_dd_sa,'Value')},'Code and phase stand-alone')))
        set(handles.plot_amb, 'Enable', 'on');
        plot_amb_Callback(handles.plot_amb, [], handles);
    end

    %enable/disable file input fields
    if(get(handles.file_type, 'SelectedObject') == handles.rinex_files);
        file_type_SelectionChangeFcn(handles.rinex_files, eventdata, handles);
    else
        file_type_SelectionChangeFcn(handles.gogps_data, eventdata, handles);
    end

    kalman_ls_Callback(handles.kalman_ls, eventdata, handles);

    %disable approximate position
    set(handles.text_approx_pos, 'Enable', 'off');
    set(handles.approx_lat, 'Enable', 'off');
    set(handles.approx_lon, 'Enable', 'off');
    set(handles.approx_h, 'Enable', 'off');
    set(handles.text_approx_lat, 'Enable', 'off');
    set(handles.text_approx_lon, 'Enable', 'off');
    set(handles.text_approx_h, 'Enable', 'off');
    set(handles.text_approx_lat_unit, 'Enable', 'off');
    set(handles.text_approx_lon_unit, 'Enable', 'off');
    set(handles.text_approx_h_unit, 'Enable', 'off');

    %disable NTRIP parameters
    set(handles.IP_address, 'Enable', 'off');
    set(handles.port, 'Enable', 'off');
    set(handles.mountpoint, 'Enable', 'off');
    set(handles.username, 'Enable', 'off');
    set(handles.password, 'Enable', 'off');
    set(handles.show_password, 'Enable', 'off');
    set(handles.text_IP_address, 'Enable', 'off');
    set(handles.text_port, 'Enable', 'off');
    set(handles.text_mountpoint, 'Enable', 'off');
    set(handles.text_username, 'Enable', 'off');
    set(handles.text_password, 'Enable', 'off');
    
    stopGOstop_Callback(handles.stopGOstop, [], handles);

end

% --- Executes during object creation, after setting all properties.
function mode_CreateFcn(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in kalman_ls.
function kalman_ls_Callback(hObject, eventdata, handles)
% hObject    handle to kalman_ls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns kalman_ls contents as cell array
%        contents{get(hObject,'Value')} returns selected item from kalman_ls
%enable Kalman filters settings
contents = cellstr(get(hObject,'String'));
if (strcmp(contents{get(hObject,'Value')},'Kalman filter'))
    cell_contents = cell(4,1);
    cell_contents{1} = 'Code stand-alone';
    cell_contents{2} = 'Code double difference';
    cell_contents{3} = 'Code and phase stand-alone';
    cell_contents{4} = 'Code and phase double difference';
    set(handles.code_dd_sa, 'String', cell_contents);

    set(handles.std_code, 'Enable', 'on');
    set(handles.text_std_code, 'Enable', 'on');
    set(handles.text_std_code_unit, 'Enable', 'on');
    set(handles.std_init, 'Enable', 'on');
    set(handles.text_std_init, 'Enable', 'on');
    set(handles.text_std_init_unit, 'Enable', 'on');

    set(handles.toggle_std_phase, 'Enable', 'on');
    set(handles.toggle_std_dtm, 'Enable', 'on');
    toggle_std_phase_Callback(handles.toggle_std_phase, eventdata, handles);
    toggle_std_dtm_Callback(handles.toggle_std_dtm, eventdata, handles);

    constraint_Callback(handles.constraint, eventdata, handles);
    code_dd_sa_Callback(handles.code_dd_sa, eventdata, handles);

    set(handles.cut_off, 'Enable', 'on');
    set(handles.text_cut_off, 'Enable', 'on');
    set(handles.text_cut_off_unit, 'Enable', 'on');
    set(handles.snr_thres, 'Enable', 'on');
    set(handles.text_snr_thres, 'Enable', 'on');
    set(handles.text_snr_thres_unit, 'Enable', 'on');
    set(handles.min_sat, 'Enable', 'on');
    set(handles.text_min_sat, 'Enable', 'on');
    set(handles.dyn_mod, 'Enable', 'on');
    %set(handles.text_dyn_mod, 'Enable', 'on');
    
    dyn_mod_Callback(handles.dyn_mod, eventdata, handles);
    
    stopGOstop_Callback(handles.stopGOstop, [], handles);
else
    cell_contents = cell(2,1);
    cell_contents{1} = 'Code stand-alone';
    cell_contents{2} = 'Code double difference';
    old_value = get(handles.code_dd_sa, 'Value');
    if (old_value == 3), set(handles.code_dd_sa, 'Value', 1); end
    if (old_value == 4), set(handles.code_dd_sa, 'Value', 2); end
    set(handles.code_dd_sa, 'String', cell_contents);

    code_dd_sa_Callback(handles.code_dd_sa, eventdata, handles);

    %disable Kalman filters settings
    set(handles.std_X, 'Enable', 'off');
    set(handles.std_Y, 'Enable', 'off');
    set(handles.std_Z, 'Enable', 'off');
    set(handles.text_std_X, 'Enable', 'off');
    set(handles.text_std_Y, 'Enable', 'off');
    set(handles.text_std_Z, 'Enable', 'off');
    set(handles.text_std_X_unit, 'Enable', 'off');
    set(handles.text_std_Y_unit, 'Enable', 'off');
    set(handles.text_std_Z_unit, 'Enable', 'off');
    set(handles.std_code, 'Enable', 'off');
    set(handles.std_phase, 'Enable', 'off');
    set(handles.text_std_code, 'Enable', 'off');
    set(handles.toggle_std_phase, 'Enable', 'off');
    set(handles.text_std_code_unit, 'Enable', 'off');
    set(handles.text_std_phase_unit, 'Enable', 'off');
    set(handles.std_init, 'Enable', 'off');
    set(handles.std_dtm, 'Enable', 'off');
    set(handles.std_vel, 'Enable', 'off');
    set(handles.text_std_init, 'Enable', 'off');
    set(handles.toggle_std_dtm, 'Enable', 'off');
    set(handles.text_std_vel, 'Enable', 'off');
    set(handles.text_std_init_unit, 'Enable', 'off');
    set(handles.text_std_dtm_unit, 'Enable', 'off');
    set(handles.text_std_vel_unit, 'Enable', 'off');
    set(handles.snr_thres, 'Enable', 'off');
    set(handles.text_snr_thres, 'Enable', 'off');
    set(handles.text_snr_thres_unit, 'Enable', 'off');
    set(handles.cs_thresh, 'Enable', 'off');
    set(handles.text_cs_thresh, 'Enable', 'off');
    set(handles.text_cs_thresh_unit, 'Enable', 'off');
    set(handles.min_sat, 'Enable', 'off');
    set(handles.text_min_sat, 'Enable', 'off');
    set(handles.antenna_h, 'Enable', 'off');
    set(handles.text_antenna_h, 'Enable', 'off');
    set(handles.text_antenna_h_unit, 'Enable', 'off');
    set(handles.dtm_path, 'Enable', 'off');
    set(handles.text_dtm_path, 'Enable', 'off');
    set(handles.browse_dtm_path, 'Enable', 'off');
    set(handles.stopGOstop, 'Enable', 'off');
    set(handles.text_stopGOstop, 'Enable', 'off');
    set(handles.dyn_mod, 'Enable', 'off');
    %set(handles.text_dyn_mod, 'Enable', 'off');
    set(handles.flag_doppler, 'Enable', 'off');
    set(handles.amb_select, 'Enable', 'off');
end

% --- Executes during object creation, after setting all properties.
function kalman_ls_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kalman_ls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in code_dd_sa.
function code_dd_sa_Callback(hObject, eventdata, handles)
% hObject    handle to code_dd_sa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns code_dd_sa contents as cell array
%        contents{get(hObject,'Value')} returns selected item from code_dd_sa
contents = cellstr(get(hObject,'String'));
if strcmp(contents{get(hObject,'Value')},'Code and phase double difference')
    check_mode = cellstr(get(handles.mode,'String'));
    if (~strcmp(check_mode{get(handles.mode,'Value')},'Real-time')) & ...
            (get(handles.plotproc,'Value'))
        set(handles.plot_amb, 'Enable', 'on');
        plot_amb_Callback(handles.plot_amb, [], handles);
    end
    set(handles.cs_thresh, 'Enable', 'on');
    set(handles.text_cs_thresh, 'Enable', 'on');
    set(handles.text_cs_thresh_unit, 'Enable', 'on');
    set(handles.toggle_std_phase, 'Enable', 'on');
    toggle_std_phase_Callback(handles.toggle_std_phase, eventdata, handles);
    set(handles.snr_thres, 'Enable', 'on');
    set(handles.text_snr_thres, 'Enable', 'on');
    set(handles.text_snr_thres_unit, 'Enable', 'on');
    set(handles.constraint, 'Enable', 'on');
    set(handles.ref_path, 'Enable', 'on');
    ref_path_Callback(handles.ref_path, eventdata, handles);
    set(handles.stopGOstop, 'Enable', 'on');
    set(handles.text_stopGOstop, 'Enable', 'on');
    stopGOstop_Callback(handles.stopGOstop, [], handles);
    cell_contents = cell(4,1);
    cell_contents{1} = 'Const. velocity';
    cell_contents{2} = 'Const. acceleration';
    cell_contents{3} = 'Static';
    cell_contents{4} = 'Variable';
    set(handles.dyn_mod, 'String', cell_contents);
    if (~strcmp(check_mode{get(handles.mode,'Value')},'Real-time'))
        set(handles.flag_doppler, 'Enable', 'on');
    end
    set(handles.amb_select, 'Enable', 'on');
elseif strcmp(contents{get(hObject,'Value')},'Code and phase stand-alone')
    check_mode = cellstr(get(handles.mode,'String'));
    if (~strcmp(check_mode{get(handles.mode,'Value')},'Real-time')) & ...
            (get(handles.plotproc,'Value'))
        set(handles.plot_amb, 'Enable', 'on');
        plot_amb_Callback(handles.plot_amb, [], handles);
    end
    set(handles.cs_thresh, 'Enable', 'on');
    set(handles.text_cs_thresh, 'Enable', 'on');
    set(handles.text_cs_thresh_unit, 'Enable', 'on');
    set(handles.toggle_std_phase, 'Enable', 'on');
    toggle_std_phase_Callback(handles.toggle_std_phase, eventdata, handles);
    set(handles.snr_thres, 'Enable', 'on');
    set(handles.text_snr_thres, 'Enable', 'on');
    set(handles.text_snr_thres_unit, 'Enable', 'on');
    set(handles.constraint, 'Value', 0);
    constraint_Callback(handles.constraint, eventdata, handles);
    set(handles.constraint, 'Enable', 'off');
    set(handles.stopGOstop, 'Enable', 'off');
    set(handles.text_stopGOstop, 'Enable', 'off');
    set(handles.dyn_mod, 'Enable', 'on');
    cell_contents = cell(3,1);
    cell_contents{1} = 'Const. velocity';
    cell_contents{2} = 'Const. acceleration';
    cell_contents{3} = 'Static';
    old_value = get(handles.dyn_mod, 'Value');
    if (old_value == 4), set(handles.dyn_mod, 'Value', 1); end
    set(handles.dyn_mod, 'String', cell_contents);
    set(handles.flag_doppler, 'Enable', 'on');
    set(handles.amb_select, 'Enable', 'off');
else
    set(handles.plot_amb, 'Enable', 'off');
    set(handles.no_skyplot_snr, 'Enable', 'on');
    set(handles.plotproc, 'Enable', 'on');
    plotproc_Callback(handles.plotproc, eventdata, handles);
    set(handles.cs_thresh, 'Enable', 'off');
    set(handles.text_cs_thresh, 'Enable', 'off');
    set(handles.text_cs_thresh_unit, 'Enable', 'off');
    set(handles.toggle_std_phase, 'Enable', 'off');
    set(handles.std_phase, 'Enable', 'off');
    set(handles.text_std_phase_unit, 'Enable', 'off');
    set(handles.constraint, 'Value', 0);
    constraint_Callback(handles.constraint, eventdata, handles);
    set(handles.constraint, 'Enable', 'off');
    set(handles.stopGOstop, 'Enable', 'off');
    set(handles.text_stopGOstop, 'Enable', 'off');
    check_KF = cellstr(get(handles.kalman_ls,'String'));
    if (strcmp(check_KF{get(handles.kalman_ls,'Value')},'Kalman filter'))
        set(handles.dyn_mod, 'Enable', 'on');
    end
    cell_contents = cell(3,1);
    cell_contents{1} = 'Const. velocity';
    cell_contents{2} = 'Const. acceleration';
    cell_contents{3} = 'Static';
    old_value = get(handles.dyn_mod, 'Value');
    if (old_value == 4), set(handles.dyn_mod, 'Value', 1); end
    set(handles.dyn_mod, 'String', cell_contents);
    set(handles.flag_doppler, 'Enable', 'off');
    set(handles.amb_select, 'Enable', 'off');
end

if strcmp(contents{get(hObject,'Value')},'Code and phase stand-alone') | ...
        strcmp(contents{get(hObject,'Value')},'Code stand-alone')
    
    set(handles.RINEX_master_obs, 'Enable', 'off');
    set(handles.text_RINEX_master_obs, 'Enable', 'off');
    set(handles.browse_master_obs, 'Enable', 'off');
    set(handles.plot_master, 'Enable', 'off');
    set(handles.master_pos, 'Enable', 'off');
    set(handles.crs, 'Enable', 'off');
    set(handles.master_X, 'Enable', 'off');
    set(handles.master_Y, 'Enable', 'off');
    set(handles.master_Z, 'Enable', 'off');
    set(handles.text_master_X, 'Enable', 'off');
    set(handles.text_master_Y, 'Enable', 'off');
    set(handles.text_master_Z, 'Enable', 'off');
    set(handles.text_master_X_unit, 'Enable', 'off');
    set(handles.text_master_Y_unit, 'Enable', 'off');
    set(handles.text_master_Z_unit, 'Enable', 'off');
    set(handles.master_lat, 'Enable', 'off');
    set(handles.master_lon, 'Enable', 'off');
    set(handles.master_h, 'Enable', 'off');
    set(handles.text_master_lat, 'Enable', 'off');
    set(handles.text_master_lon, 'Enable', 'off');
    set(handles.text_master_h, 'Enable', 'off');
    set(handles.text_master_lat_unit, 'Enable', 'off');
    set(handles.text_master_lon_unit, 'Enable', 'off');
    set(handles.text_master_h_unit, 'Enable', 'off');
else
    contents = cellstr(get(handles.mode,'String'));
    if(get(handles.file_type, 'SelectedObject') == handles.rinex_files & ~strcmp(contents{get(handles.mode,'Value')},'Real-time'));
        set(handles.RINEX_master_obs, 'Enable', 'on');
        set(handles.text_RINEX_master_obs, 'Enable', 'on');
        set(handles.browse_master_obs, 'Enable', 'on');
    end
    set(handles.plot_master, 'Enable', 'on');
    set(handles.master_pos, 'Enable', 'on');
    master_pos_Callback(handles.master_pos, [], handles);
end

% --- Executes during object creation, after setting all properties.
function code_dd_sa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to code_dd_sa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nav_mon.
function nav_mon_Callback(hObject, eventdata, handles)
% hObject    handle to nav_mon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns nav_mon contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nav_mon

contents = cellstr(get(hObject,'String'));
if (strcmp(contents{get(hObject,'Value')},'Navigation'))

    %enable options
    set(handles.constraint, 'Enable', 'on');
    set(handles.ref_path, 'Enable', 'on');
    ref_path_Callback(handles.ref_path, eventdata, handles);
    set(handles.plot_master, 'Enable', 'on');
    set(handles.err_ellipse, 'Enable', 'on');
    set(handles.google_earth, 'Enable', 'on');
    set(handles.text_num_receivers, 'Enable', 'off');
    set(handles.num_receivers, 'Value', 1);
    set(handles.num_receivers, 'Enable', 'off');
    set(handles.com_select_0, 'Enable', 'on');
    set(handles.com_select_1, 'Enable', 'off');
    set(handles.com_select_2, 'Enable', 'off');
    set(handles.com_select_3, 'Enable', 'off');
    set(handles.protocol_select_0, 'Enable', 'on');
    set(handles.protocol_select_1, 'Enable', 'off');
    set(handles.protocol_select_2, 'Enable', 'off');
    set(handles.protocol_select_3, 'Enable', 'off');
    set(handles.use_ntrip, 'Enable', 'on');
    set(handles.no_skyplot_snr, 'Enable', 'on');
    
    set(handles.gogps_data_output, 'Enable', 'on');
    set(handles.text_gogps_data_output, 'Enable', 'on');
    set(handles.browse_gogps_data_output, 'Enable', 'on');
    set(handles.gogps_data_output_prefix, 'Enable', 'on');
    set(handles.text_gogps_data_output_prefix, 'Enable', 'on');

    kalman_ls_Callback(handles.kalman_ls, eventdata, handles);
    
    set(handles.plotproc, 'Enable', 'on');
    plotproc_Callback(handles.plotproc, eventdata, handles);
    
    set(handles.stopGOstop, 'Enable', 'on');
    set(handles.text_stopGOstop, 'Enable', 'on');
    stopGOstop_Callback(handles.stopGOstop, [], handles);
    
    cell_contents = cell(4,1);
    cell_contents{1} = 'Const. velocity';
    cell_contents{2} = 'Const. acceleration';
    cell_contents{3} = 'Static';
    cell_contents{4} = 'Variable';
    set(handles.dyn_mod, 'String', cell_contents);
    
    %enable weights
    set(handles.weight_0, 'Enable', 'on');
    set(handles.weight_1, 'Enable', 'on');
    set(handles.weight_2, 'Enable', 'on');
    set(handles.weight_3, 'Enable', 'on');

    %enable master connection
    set(handles.IP_address, 'Enable', 'on');
    set(handles.port, 'Enable', 'on');
    set(handles.text_IP_address, 'Enable', 'on');
    set(handles.text_port, 'Enable', 'on');
    use_ntrip_Callback(handles.use_ntrip, eventdata, handles);

    %disable approximate position
    set(handles.text_approx_pos, 'Enable', 'off');
    set(handles.approx_lat, 'Enable', 'off');
    set(handles.approx_lon, 'Enable', 'off');
    set(handles.approx_h, 'Enable', 'off');
    set(handles.text_approx_lat, 'Enable', 'off');
    set(handles.text_approx_lon, 'Enable', 'off');
    set(handles.text_approx_h, 'Enable', 'off');
    set(handles.text_approx_lat_unit, 'Enable', 'off');
    set(handles.text_approx_lon_unit, 'Enable', 'off');
    set(handles.text_approx_h_unit, 'Enable', 'off');

else

    %disable some options
    set(handles.constraint, 'Enable', 'off');
    set(handles.ref_path, 'Enable', 'off');
    set(handles.plot_master, 'Enable', 'off');
    set(handles.err_ellipse, 'Enable', 'off');
    set(handles.google_earth, 'Enable', 'off');
    set(handles.no_skyplot_snr, 'Enable', 'off');
    set(handles.plotproc, 'Enable', 'off');
    set(handles.plot_amb, 'Enable', 'off');
    %set(handles.gogps_data_output, 'Enable', 'off');
    %set(handles.text_gogps_data_output, 'Enable', 'off');
    %set(handles.browse_gogps_data_output, 'Enable', 'off');
    %set(handles.gogps_data_output_prefix, 'Enable', 'off');
    %set(handles.text_gogps_data_output_prefix, 'Enable', 'off');
    set(handles.ref_path_input, 'Enable', 'off');
    set(handles.text_ref_path_input, 'Enable', 'off');
    set(handles.browse_ref_path_input, 'Enable', 'off');

    %disable Kalman filters settings
    set(handles.std_X, 'Enable', 'off');
    set(handles.std_Y, 'Enable', 'off');
    set(handles.std_Z, 'Enable', 'off');
    set(handles.text_std_X, 'Enable', 'off');
    set(handles.text_std_Y, 'Enable', 'off');
    set(handles.text_std_Z, 'Enable', 'off');
    set(handles.text_std_X_unit, 'Enable', 'off');
    set(handles.text_std_Y_unit, 'Enable', 'off');
    set(handles.text_std_Z_unit, 'Enable', 'off');
    set(handles.std_code, 'Enable', 'off');
    set(handles.std_phase, 'Enable', 'off');
    set(handles.text_std_code, 'Enable', 'off');
    set(handles.toggle_std_phase, 'Enable', 'off');
    set(handles.text_std_code_unit, 'Enable', 'off');
    set(handles.text_std_phase_unit, 'Enable', 'off');
    set(handles.std_init, 'Enable', 'off');
    set(handles.std_dtm, 'Enable', 'off');
    set(handles.std_vel, 'Enable', 'off');
    set(handles.text_std_init, 'Enable', 'off');
    set(handles.toggle_std_dtm, 'Enable', 'off');
    set(handles.text_std_vel, 'Enable', 'off');
    set(handles.text_std_init_unit, 'Enable', 'off');
    set(handles.text_std_dtm_unit, 'Enable', 'off');
    set(handles.text_std_vel_unit, 'Enable', 'off');
    set(handles.cs_thresh, 'Enable', 'off');
    set(handles.text_cs_thresh, 'Enable', 'off');
    set(handles.text_cs_thresh_unit, 'Enable', 'off');
    set(handles.cut_off, 'Enable', 'off');
    set(handles.text_cut_off, 'Enable', 'off');
    set(handles.text_cut_off_unit, 'Enable', 'off');
    set(handles.snr_thres, 'Enable', 'off');
    set(handles.text_snr_thres, 'Enable', 'off');
    set(handles.text_snr_thres_unit, 'Enable', 'off');
    set(handles.min_sat, 'Enable', 'off');
    set(handles.text_min_sat, 'Enable', 'off');
    set(handles.antenna_h, 'Enable', 'off');
    set(handles.text_antenna_h, 'Enable', 'off');
    set(handles.text_antenna_h_unit, 'Enable', 'off');
    set(handles.dtm_path, 'Enable', 'off');
    set(handles.text_dtm_path, 'Enable', 'off');
    set(handles.browse_dtm_path, 'Enable', 'off');
    set(handles.weight_0, 'Enable', 'off');
    set(handles.weight_1, 'Enable', 'off');
    set(handles.weight_2, 'Enable', 'off');
    set(handles.weight_3, 'Enable', 'off');
    set(handles.master_pos, 'Enable', 'off');
    set(handles.crs, 'Enable', 'off');
    set(handles.master_X, 'Enable', 'off');
    set(handles.master_Y, 'Enable', 'off');
    set(handles.master_Z, 'Enable', 'off');
    set(handles.text_master_X, 'Enable', 'off');
    set(handles.text_master_Y, 'Enable', 'off');
    set(handles.text_master_Z, 'Enable', 'off');
    set(handles.text_master_X_unit, 'Enable', 'off');
    set(handles.text_master_Y_unit, 'Enable', 'off');
    set(handles.text_master_Z_unit, 'Enable', 'off');
    set(handles.master_lat, 'Enable', 'off');
    set(handles.master_lon, 'Enable', 'off');
    set(handles.master_h, 'Enable', 'off');
    set(handles.text_master_lat, 'Enable', 'off');
    set(handles.text_master_lon, 'Enable', 'off');
    set(handles.text_master_h, 'Enable', 'off');
    set(handles.text_master_lat_unit, 'Enable', 'off');
    set(handles.text_master_lon_unit, 'Enable', 'off');
    set(handles.text_master_h_unit, 'Enable', 'off');
    set(handles.flag_doppler, 'Enable', 'off');
    set(handles.amb_select, 'Enable', 'off');
    
    cell_contents = cell(2,1);
    cell_contents{1} = 'Constant';
    cell_contents{2} = 'Variable';
    old_value = get(handles.dyn_mod, 'Value');
    if (old_value == 3), set(handles.dyn_mod, 'Value', 1); end
    if (old_value == 4), set(handles.dyn_mod, 'Value', 2); end
    set(handles.dyn_mod, 'String', cell_contents);

    if (strcmp(contents{get(hObject,'Value')},'Rover monitor'))

        num_receivers_Callback(handles.num_receivers, [], handles);
        set(handles.use_ntrip, 'Enable', 'off');

        %disable approximate position
        set(handles.text_approx_pos, 'Enable', 'off');
        set(handles.approx_lat, 'Enable', 'off');
        set(handles.approx_lon, 'Enable', 'off');
        set(handles.approx_h, 'Enable', 'off');
        set(handles.text_approx_lat, 'Enable', 'off');
        set(handles.text_approx_lon, 'Enable', 'off');
        set(handles.text_approx_h, 'Enable', 'off');
        set(handles.text_approx_lat_unit, 'Enable', 'off');
        set(handles.text_approx_lon_unit, 'Enable', 'off');
        set(handles.text_approx_h_unit, 'Enable', 'off');

        %disable NTRIP parameters
        set(handles.IP_address, 'Enable', 'off');
        set(handles.port, 'Enable', 'off');
        set(handles.mountpoint, 'Enable', 'off');
        set(handles.username, 'Enable', 'off');
        set(handles.password, 'Enable', 'off');
        set(handles.show_password, 'Enable', 'off');
        set(handles.text_IP_address, 'Enable', 'off');
        set(handles.text_port, 'Enable', 'off');
        set(handles.text_mountpoint, 'Enable', 'off');
        set(handles.text_username, 'Enable', 'off');
        set(handles.text_password, 'Enable', 'off');
        
        set(handles.stopGOstop, 'Enable', 'on');
        set(handles.text_stopGOstop, 'Enable', 'on');
        stopGOstop_Callback(handles.stopGOstop, [], handles);

    elseif (strcmp(contents{get(hObject,'Value')},'Master monitor'))

        set(handles.text_num_receivers, 'Enable', 'off');
        set(handles.num_receivers, 'Enable', 'off');        
        set(handles.com_select_0, 'Enable', 'off');
        set(handles.com_select_1, 'Enable', 'off');
        set(handles.com_select_2, 'Enable', 'off');
        set(handles.com_select_3, 'Enable', 'off');
        set(handles.protocol_select_0, 'Enable', 'off');
        set(handles.protocol_select_1, 'Enable', 'off');
        set(handles.protocol_select_2, 'Enable', 'off');
        set(handles.protocol_select_3, 'Enable', 'off');        
        set(handles.use_ntrip, 'Enable', 'on');
        set(handles.dyn_mod, 'Enable', 'off');
        %set(handles.text_dyn_mod, 'Enable', 'off');

        %enable master connection
        set(handles.IP_address, 'Enable', 'on');
        set(handles.port, 'Enable', 'on');
        set(handles.text_IP_address, 'Enable', 'on');
        set(handles.text_port, 'Enable', 'on');
        use_ntrip_Callback(handles.use_ntrip, eventdata, handles);

        %enable approximate position
        set(handles.text_approx_pos, 'Enable', 'on');
        set(handles.approx_lat, 'Enable', 'on');
        set(handles.approx_lon, 'Enable', 'on');
        set(handles.approx_h, 'Enable', 'on');
        set(handles.text_approx_lat, 'Enable', 'on');
        set(handles.text_approx_lon, 'Enable', 'on');
        set(handles.text_approx_h, 'Enable', 'on');
        set(handles.text_approx_lat_unit, 'Enable', 'on');
        set(handles.text_approx_lon_unit, 'Enable', 'on');
        set(handles.text_approx_h_unit, 'Enable', 'on');
        
        set(handles.stopGOstop, 'Enable', 'off');
        set(handles.text_stopGOstop, 'Enable', 'off');

    elseif (strcmp(contents{get(hObject,'Value')},'Rover and Master monitor'))

        set(handles.text_num_receivers, 'Enable', 'off');
        set(handles.num_receivers, 'Value', 1);
        set(handles.num_receivers, 'Enable', 'off');
        set(handles.com_select_0, 'Enable', 'on');
        set(handles.com_select_1, 'Enable', 'off');
        set(handles.com_select_2, 'Enable', 'off');
        set(handles.com_select_3, 'Enable', 'off');
        set(handles.protocol_select_0, 'Enable', 'on');
        set(handles.protocol_select_1, 'Enable', 'off');
        set(handles.protocol_select_2, 'Enable', 'off');
        set(handles.protocol_select_3, 'Enable', 'off');
        set(handles.use_ntrip, 'Enable', 'on');
        
        %enable master connection
        set(handles.IP_address, 'Enable', 'on');
        set(handles.port, 'Enable', 'on');
        set(handles.text_IP_address, 'Enable', 'on');
        set(handles.text_port, 'Enable', 'on');
        use_ntrip_Callback(handles.use_ntrip, eventdata, handles);

        %disable approximate position
        set(handles.text_approx_pos, 'Enable', 'off');
        set(handles.approx_lat, 'Enable', 'off');
        set(handles.approx_lon, 'Enable', 'off');
        set(handles.approx_h, 'Enable', 'off');
        set(handles.text_approx_lat, 'Enable', 'off');
        set(handles.text_approx_lon, 'Enable', 'off');
        set(handles.text_approx_h, 'Enable', 'off');
        set(handles.text_approx_lat_unit, 'Enable', 'off');
        set(handles.text_approx_lon_unit, 'Enable', 'off');
        set(handles.text_approx_h_unit, 'Enable', 'off');
        
        set(handles.stopGOstop, 'Enable', 'on');
        set(handles.text_stopGOstop, 'Enable', 'on');
        stopGOstop_Callback(handles.stopGOstop, [], handles);
    end
end


% --- Executes during object creation, after setting all properties.
function nav_mon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nav_mon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in constraint.
function constraint_Callback(hObject, eventdata, handles)
% hObject    handle to constraint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of constraint
if (get(hObject,'Value'))
    set(handles.err_ellipse, 'Enable', 'off');
    set(handles.std_vel, 'Enable', 'on');
    set(handles.text_std_vel, 'Enable', 'on');
    set(handles.text_std_vel_unit, 'Enable', 'on');
else
    if (get(handles.plotproc,'Value'))
        set(handles.err_ellipse, 'Enable', 'on');
    end
    set(handles.std_vel, 'Enable', 'off');
    set(handles.text_std_vel, 'Enable', 'off');
    set(handles.text_std_vel_unit, 'Enable', 'off');
end

% --- Executes on button press in ref_path.
function ref_path_Callback(hObject, eventdata, handles)
% hObject    handle to ref_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ref_path
contents = cellstr(get(handles.code_dd_sa,'String'));
if (get(hObject,'Value')) & (strcmp(contents{get(handles.code_dd_sa,'Value')},'Code and phase double difference'))
    set(handles.ref_path_input, 'Enable', 'on');
    set(handles.text_ref_path_input, 'Enable', 'on');
    set(handles.browse_ref_path_input, 'Enable', 'on');
    set(handles.constraint, 'Enable', 'on');
else
    set(handles.ref_path_input, 'Enable', 'off');
    set(handles.text_ref_path_input, 'Enable', 'off');
    set(handles.browse_ref_path_input, 'Enable', 'off');
    set(handles.constraint, 'Value', 0);
    constraint_Callback(handles.constraint, eventdata, handles);
    set(handles.constraint, 'Enable', 'off');
end

% --- Executes on button press in plot_master.
function plot_master_Callback(hObject, eventdata, handles)
% hObject    handle to plot_master (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_master


% --- Executes on button press in err_ellipse.
function err_ellipse_Callback(hObject, eventdata, handles)
% hObject    handle to err_ellipse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of err_ellipse


% --- Executes on button press in google_earth.
function google_earth_Callback(hObject, eventdata, handles)
% hObject    handle to google_earth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of google_earth

% --- Executes on button press in use_ntrip.
function use_ntrip_Callback(hObject, eventdata, handles)
% hObject    handle to use_ntrip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_ntrip

if (get(hObject,'Value'))
    %enable NTRIP parameters
    set(handles.mountpoint, 'Enable', 'on');
    set(handles.username, 'Enable', 'on');
    set(handles.password, 'Enable', 'on');
    set(handles.show_password, 'Enable', 'on');
    set(handles.text_mountpoint, 'Enable', 'on');
    set(handles.text_username, 'Enable', 'on');
    set(handles.text_password, 'Enable', 'on');
else
    %disable NTRIP parameters
    set(handles.mountpoint, 'Enable', 'off');
    set(handles.username, 'Enable', 'off');
    set(handles.password, 'Enable', 'off');
    set(handles.show_password, 'Enable', 'off');
    set(handles.text_mountpoint, 'Enable', 'off');
    set(handles.text_username, 'Enable', 'off');
    set(handles.text_password, 'Enable', 'off');

    %disable approximate position
    set(handles.text_approx_pos, 'Enable', 'off');
    set(handles.approx_lat, 'Enable', 'off');
    set(handles.approx_lon, 'Enable', 'off');
    set(handles.approx_h, 'Enable', 'off');
    set(handles.text_approx_lat, 'Enable', 'off');
    set(handles.text_approx_lon, 'Enable', 'off');
    set(handles.text_approx_h, 'Enable', 'off');
    set(handles.text_approx_lat_unit, 'Enable', 'off');
    set(handles.text_approx_lon_unit, 'Enable', 'off');
    set(handles.text_approx_h_unit, 'Enable', 'off');
end

% --- Executes on button press in master_pos.
function master_pos_Callback(hObject, eventdata, handles)
% hObject    handle to master_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of master_pos
if (get(hObject,'Value'))
    set(handles.crs, 'Enable', 'off');
    set(handles.master_X, 'Enable', 'off');
    set(handles.master_Y, 'Enable', 'off');
    set(handles.master_Z, 'Enable', 'off');
    set(handles.text_master_X, 'Enable', 'off');
    set(handles.text_master_Y, 'Enable', 'off');
    set(handles.text_master_Z, 'Enable', 'off');
    set(handles.text_master_X_unit, 'Enable', 'off');
    set(handles.text_master_Y_unit, 'Enable', 'off');
    set(handles.text_master_Z_unit, 'Enable', 'off');
    set(handles.master_lat, 'Enable', 'off');
    set(handles.master_lon, 'Enable', 'off');
    set(handles.master_h, 'Enable', 'off');
    set(handles.text_master_lat, 'Enable', 'off');
    set(handles.text_master_lon, 'Enable', 'off');
    set(handles.text_master_h, 'Enable', 'off');
    set(handles.text_master_lat_unit, 'Enable', 'off');
    set(handles.text_master_lon_unit, 'Enable', 'off');
    set(handles.text_master_h_unit, 'Enable', 'off');
else
    set(handles.crs, 'Enable', 'on');
    crs_Callback(handles.crs, eventdata, handles);
end


function master_X_Callback(hObject, eventdata, handles)
% hObject    handle to master_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of master_X as text
%        str2double(get(hObject,'String')) returns contents of master_X as a double


% --- Executes during object creation, after setting all properties.
function master_X_CreateFcn(hObject, eventdata, handles)
% hObject    handle to master_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function master_Y_Callback(hObject, eventdata, handles)
% hObject    handle to master_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of master_Y as text
%        str2double(get(hObject,'String')) returns contents of master_Y as a double


% --- Executes during object creation, after setting all properties.
function master_Y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to master_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function master_Z_Callback(hObject, eventdata, handles)
% hObject    handle to master_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of master_Z as text
%        str2double(get(hObject,'String')) returns contents of master_Z as a double


% --- Executes during object creation, after setting all properties.
function master_Z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to master_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function master_lat_Callback(hObject, eventdata, handles)
% hObject    handle to master_lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of master_lat as text
%        str2double(get(hObject,'String')) returns contents of master_lat as a double


% --- Executes during object creation, after setting all properties.
function master_lat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to master_lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function master_lon_Callback(hObject, eventdata, handles)
% hObject    handle to master_lon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of master_lon as text
%        str2double(get(hObject,'String')) returns contents of master_lon as a double


% --- Executes during object creation, after setting all properties.
function master_lon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to master_lon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function master_h_Callback(hObject, eventdata, handles)
% hObject    handle to master_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of master_h as text
%        str2double(get(hObject,'String')) returns contents of master_h as a double


% --- Executes during object creation, after setting all properties.
function master_h_CreateFcn(hObject, eventdata, handles)
% hObject    handle to master_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in crs.
function crs_Callback(hObject, eventdata, handles)
% hObject    handle to crs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns crs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from crs
contents = cellstr(get(hObject,'String'));
if (strcmp(contents{get(hObject,'Value')},'ECEF (X,Y,Z)'))
    set(handles.master_X, 'Enable', 'on');
    set(handles.master_Y, 'Enable', 'on');
    set(handles.master_Z, 'Enable', 'on');
    set(handles.text_master_X, 'Enable', 'on');
    set(handles.text_master_Y, 'Enable', 'on');
    set(handles.text_master_Z, 'Enable', 'on');
    set(handles.text_master_X_unit, 'Enable', 'on');
    set(handles.text_master_Y_unit, 'Enable', 'on');
    set(handles.text_master_Z_unit, 'Enable', 'on');
    set(handles.master_lat, 'Enable', 'off');
    set(handles.master_lon, 'Enable', 'off');
    set(handles.master_h, 'Enable', 'off');
    set(handles.text_master_lat, 'Enable', 'off');
    set(handles.text_master_lon, 'Enable', 'off');
    set(handles.text_master_h, 'Enable', 'off');
    set(handles.text_master_lat_unit, 'Enable', 'off');
    set(handles.text_master_lon_unit, 'Enable', 'off');
    set(handles.text_master_h_unit, 'Enable', 'off');
else
    set(handles.master_X, 'Enable', 'off');
    set(handles.master_Y, 'Enable', 'off');
    set(handles.master_Z, 'Enable', 'off');
    set(handles.text_master_X, 'Enable', 'off');
    set(handles.text_master_Y, 'Enable', 'off');
    set(handles.text_master_Z, 'Enable', 'off');
    set(handles.text_master_X_unit, 'Enable', 'off');
    set(handles.text_master_Y_unit, 'Enable', 'off');
    set(handles.text_master_Z_unit, 'Enable', 'off');
    set(handles.master_lat, 'Enable', 'on');
    set(handles.master_lon, 'Enable', 'on');
    set(handles.master_h, 'Enable', 'on');
    set(handles.text_master_lat, 'Enable', 'on');
    set(handles.text_master_lon, 'Enable', 'on');
    set(handles.text_master_h, 'Enable', 'on');
    set(handles.text_master_lat_unit, 'Enable', 'on');
    set(handles.text_master_lon_unit, 'Enable', 'on');
    set(handles.text_master_h_unit, 'Enable', 'on');
end

% --- Executes during object creation, after setting all properties.
function crs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to crs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function gogps_data_input_Callback(hObject, eventdata, handles)
% hObject    handle to gogps_data_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gogps_data_input as text
%        str2double(get(hObject,'String')) returns contents of gogps_data_input as a double


% --- Executes during object creation, after setting all properties.
function gogps_data_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gogps_data_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_gogps_input.
function browse_gogps_input_Callback(hObject, eventdata, handles)
% hObject    handle to browse_gogps_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {'*.bin','goGPS binary data (*.bin)'}, ...
    'Choose goGPS binary data','../data');

if (filename ~= 0)
    pos = find(filename == '_');
    filename = filename(1:pos(end-1)-1);
    set(handles.gogps_data_input,'String',fullfile(pathname, filename));
end


function gogps_data_output_Callback(hObject, eventdata, handles)
% hObject    handle to gogps_data_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gogps_data_output as text
%        str2double(get(hObject,'String')) returns contents of gogps_data_output as a double


% --- Executes during object creation, after setting all properties.
function gogps_data_output_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gogps_data_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_gogps_data_output.
function browse_gogps_data_output_Callback(hObject, eventdata, handles)
% hObject    handle to browse_gogps_data_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dname = uigetdir('../data/','Choose a directory to store goGPS data');
if (dname ~= 0)
    set(handles.gogps_data_output,'String',dname);
end

function RINEX_rover_obs_Callback(hObject, eventdata, handles)
% hObject    handle to RINEX_rover_obs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RINEX_rover_obs as text
%        str2double(get(hObject,'String')) returns contents of RINEX_rover_obs as a double


% --- Executes during object creation, after setting all properties.
function RINEX_rover_obs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RINEX_rover_obs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_rover_obs.
function browse_rover_obs_Callback(hObject, eventdata, handles)
% hObject    handle to browse_rover_obs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {'*.obs;*.??o','RINEX observation files (*.obs,*.??o)';
    '*.obs','Observation files (*.obs)'; ...
    '*.??o','Observation files (*.??o)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose a RINEX observation file for the rover','../data/data_RINEX');

if (filename ~= 0)
    set(handles.RINEX_rover_obs,'String',fullfile(pathname, filename));
end


function RINEX_master_obs_Callback(hObject, eventdata, handles)
% hObject    handle to RINEX_master_obs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RINEX_master_obs as text
%        str2double(get(hObject,'String')) returns contents of RINEX_master_obs as a double


% --- Executes during object creation, after setting all properties.
function RINEX_master_obs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RINEX_master_obs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_master_obs.
function browse_master_obs_Callback(hObject, eventdata, handles)
% hObject    handle to browse_master_obs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {'*.obs;*.??o','RINEX observation files (*.obs,*.??o)';
    '*.obs','Observation files (*.obs)'; ...
    '*.??o','Observation files (*.??o)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose a RINEX observation file for the master','../data/data_RINEX');

if (filename ~= 0)
    set(handles.RINEX_master_obs,'String',fullfile(pathname, filename));
end


function RINEX_nav_Callback(hObject, eventdata, handles)
% hObject    handle to RINEX_nav (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RINEX_nav as text
%        str2double(get(hObject,'String')) returns contents of RINEX_nav as a double


% --- Executes during object creation, after setting all properties.
function RINEX_nav_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RINEX_nav (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_nav.
function browse_nav_Callback(hObject, eventdata, handles)
% hObject    handle to browse_nav (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {'*.nav;*.??n','RINEX navigation files (*.nav,*.??n)';
    '*.nav','Navigation files (*.obs)'; ...
    '*.??n','Navigation files (*.??o)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose a RINEX navigation file for the master','../data/data_RINEX');

if (filename ~= 0)
    set(handles.RINEX_nav,'String',fullfile(pathname, filename));
end


% --- Executes when selected object is changed in file_type.
function file_type_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in file_type
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if (strcmp(get(handles.rinex_files, 'Enable'),'on'))
    if (hObject == handles.rinex_files)
        
        set(handles.RINEX_rover_obs, 'Enable', 'on');
        set(handles.RINEX_nav, 'Enable', 'on');
        set(handles.browse_rover_obs, 'Enable', 'on');
        set(handles.browse_nav, 'Enable', 'on');
        set(handles.text_RINEX_rover_obs, 'Enable', 'on');
        set(handles.text_RINEX_nav, 'Enable', 'on');
        
        code_dd_sa_Callback(handles.code_dd_sa, [], handles);
        
        set(handles.gogps_data_input, 'Enable', 'off');
        set(handles.browse_gogps_input, 'Enable', 'off');
        set(handles.text_gogps_input, 'Enable', 'off');
        
        set(handles.stopGOstop, 'Enable', 'off');
        set(handles.stopGOstop, 'Value', 0);
        set(handles.text_stopGOstop, 'Enable', 'off');
        stopGOstop_Callback(handles.stopGOstop, [], handles);
        
        cell_contents = cell(3,1);
        cell_contents{1} = 'Const. velocity';
        cell_contents{2} = 'Const. acceleration';
        cell_contents{3} = 'Static';
        old_value = get(handles.dyn_mod, 'Value');
        if (old_value == 4), set(handles.dyn_mod, 'Value', 1); end
        set(handles.dyn_mod, 'String', cell_contents);
        
    else
        set(handles.RINEX_rover_obs, 'Enable', 'off');
        set(handles.RINEX_master_obs, 'Enable', 'off');
        set(handles.RINEX_nav, 'Enable', 'off');
        set(handles.browse_rover_obs, 'Enable', 'off');
        set(handles.browse_master_obs, 'Enable', 'off');
        set(handles.browse_nav, 'Enable', 'off');
        set(handles.text_RINEX_rover_obs, 'Enable', 'off');
        set(handles.text_RINEX_master_obs, 'Enable', 'off');
        set(handles.text_RINEX_nav, 'Enable', 'off');
        
        set(handles.gogps_data_input, 'Enable', 'on');
        set(handles.browse_gogps_input, 'Enable', 'on');
        set(handles.text_gogps_input, 'Enable', 'on');
        
        set(handles.stopGOstop, 'Enable', 'on');
        set(handles.text_stopGOstop, 'Enable', 'on');
        stopGOstop_Callback(handles.stopGOstop, [], handles);
        
        code_dd_sa_Callback(handles.code_dd_sa, [], handles);
    end
end


function gogps_data_output_prefix_Callback(hObject, eventdata, handles)
% hObject    handle to gogps_data_output_prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gogps_data_output_prefix as text
%        str2double(get(hObject,'String')) returns contents of gogps_data_output_prefix as a double


% --- Executes during object creation, after setting all properties.
function gogps_data_output_prefix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gogps_data_output_prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_amb.
function plot_amb_Callback(hObject, eventdata, handles)
% hObject    handle to plot_amb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_amb
if (get(hObject,'Value'))
    set(handles.no_skyplot_snr, 'Enable', 'off');
else
    set(handles.no_skyplot_snr, 'Enable', 'on');
end


function snr_thres_Callback(hObject, eventdata, handles)
% hObject    handle to snr_thres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of snr_thres as text
%        str2double(get(hObject,'String')) returns contents of snr_thres as a double


% --- Executes during object creation, after setting all properties.
function snr_thres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to snr_thres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cut_off_Callback(hObject, eventdata, handles)
% hObject    handle to cut_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cut_off as text
%        str2double(get(hObject,'String')) returns contents of cut_off as a double


% --- Executes during object creation, after setting all properties.
function cut_off_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cut_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function min_sat_Callback(hObject, eventdata, handles)
% hObject    handle to min_sat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_sat as text
%        str2double(get(hObject,'String')) returns contents of min_sat as a double


% --- Executes during object creation, after setting all properties.
function min_sat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_sat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cs_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to cs_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cs_thresh as text
%        str2double(get(hObject,'String')) returns contents of cs_thresh as a double


% --- Executes during object creation, after setting all properties.
function cs_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cs_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function std_vel_Callback(hObject, eventdata, handles)
% hObject    handle to std_vel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_vel as text
%        str2double(get(hObject,'String')) returns contents of std_vel as a double


% --- Executes during object creation, after setting all properties.
function std_vel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to std_vel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in toggle_std_phase.
function toggle_std_phase_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_std_phase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggle_std_phase
if (get(hObject,'Value'))
    set(handles.std_phase, 'Enable', 'on');
    set(handles.text_std_phase_unit, 'Enable', 'on');
else
    set(handles.std_phase, 'Enable', 'off');
    set(handles.text_std_phase_unit, 'Enable', 'off');
end


function std_phase_Callback(hObject, eventdata, handles)
% hObject    handle to std_phase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_phase as text
%        str2double(get(hObject,'String')) returns contents of std_phase as a double


% --- Executes during object creation, after setting all properties.
function std_phase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to std_phase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function std_code_Callback(hObject, eventdata, handles)
% hObject    handle to std_code (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_code as text
%        str2double(get(hObject,'String')) returns contents of std_code as a double


% --- Executes during object creation, after setting all properties.
function std_code_CreateFcn(hObject, eventdata, handles)
% hObject    handle to std_code (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in toggle_std_dtm.
function toggle_std_dtm_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_std_dtm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggle_std_dtm
if (get(hObject,'Value'))
    set(handles.std_dtm, 'Enable', 'on');
    set(handles.text_std_dtm_unit, 'Enable', 'on');
    set(handles.antenna_h, 'Enable', 'on');
    set(handles.text_antenna_h, 'Enable', 'on');
    set(handles.text_antenna_h_unit, 'Enable', 'on');
    set(handles.dtm_path, 'Enable', 'on');
    set(handles.text_dtm_path, 'Enable', 'on');
    set(handles.browse_dtm_path, 'Enable', 'on');
else
    set(handles.std_dtm, 'Enable', 'off');
    set(handles.text_std_dtm_unit, 'Enable', 'off');
    set(handles.antenna_h, 'Enable', 'off');
    set(handles.text_antenna_h, 'Enable', 'off');
    set(handles.text_antenna_h_unit, 'Enable', 'off');
    set(handles.dtm_path, 'Enable', 'off');
    set(handles.text_dtm_path, 'Enable', 'off');
    set(handles.browse_dtm_path, 'Enable', 'off');
end


function std_dtm_Callback(hObject, eventdata, handles)
% hObject    handle to std_dtm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_dtm as text
%        str2double(get(hObject,'String')) returns contents of std_dtm as a double


% --- Executes during object creation, after setting all properties.
function std_dtm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to std_dtm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function std_init_Callback(hObject, eventdata, handles)
% hObject    handle to std_init (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_init as text
%        str2double(get(hObject,'String')) returns contents of std_init as a double


% --- Executes during object creation, after setting all properties.
function std_init_CreateFcn(hObject, eventdata, handles)
% hObject    handle to std_init (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function std_Z_Callback(hObject, eventdata, handles)
% hObject    handle to std_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_Z as text
%        str2double(get(hObject,'String')) returns contents of std_Z as a double


% --- Executes during object creation, after setting all properties.
function std_Z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to std_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function std_Y_Callback(hObject, eventdata, handles)
% hObject    handle to std_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_Y as text
%        str2double(get(hObject,'String')) returns contents of std_Y as a double


% --- Executes during object creation, after setting all properties.
function std_Y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to std_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function std_X_Callback(hObject, eventdata, handles)
% hObject    handle to std_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_X as text
%        str2double(get(hObject,'String')) returns contents of std_X as a double


% --- Executes during object creation, after setting all properties.
function std_X_CreateFcn(hObject, eventdata, handles)
% hObject    handle to std_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dtm_path_Callback(hObject, eventdata, handles)
% hObject    handle to dtm_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dtm_path as text
%        str2double(get(hObject,'String')) returns contents of dtm_path as a double


% --- Executes during object creation, after setting all properties.
function dtm_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dtm_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_dtm_path.
function browse_dtm_path_Callback(hObject, eventdata, handles)
% hObject    handle to browse_dtm_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dname = uigetdir('../data/dtm','Choose a directory containing DTM data');
if (dname ~= 0)
    set(handles.dtm_path,'String',dname);
end

function ref_path_input_Callback(hObject, eventdata, handles)
% hObject    handle to ref_path_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ref_path_input as text
%        str2double(get(hObject,'String')) returns contents of ref_path_input as a double


% --- Executes during object creation, after setting all properties.
function ref_path_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ref_path_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_ref_path_input.
function browse_ref_path_input_Callback(hObject, eventdata, handles)
% hObject    handle to browse_ref_path_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile('*.mat', 'Choose file containing reference path','../data');

if (filename ~= 0)
    set(handles.ref_path_input,'String',fullfile(pathname, filename));
end


function IP_address_Callback(hObject, eventdata, handles)
% hObject    handle to IP_address (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IP_address as text
%        str2double(get(hObject,'String')) returns contents of IP_address as a double


% --- Executes during object creation, after setting all properties.
function IP_address_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IP_address (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function port_Callback(hObject, eventdata, handles)
% hObject    handle to port (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of port as text
%        str2double(get(hObject,'String')) returns contents of port as a double


% --- Executes during object creation, after setting all properties.
function port_CreateFcn(hObject, eventdata, handles)
% hObject    handle to port (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function username_Callback(hObject, eventdata, handles)
% hObject    handle to username (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of username as text
%        str2double(get(hObject,'String')) returns contents of username as a double


% --- Executes during object creation, after setting all properties.
function username_CreateFcn(hObject, eventdata, handles)
% hObject    handle to username (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function password_Callback(hObject, eventdata, handles)
% hObject    handle to password (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of password as text
%        str2double(get(hObject,'String')) returns contents of password as a double


% --- Executes during object creation, after setting all properties.
function password_CreateFcn(hObject, eventdata, handles)
% hObject    handle to password (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mountpoint_Callback(hObject, eventdata, handles)
% hObject    handle to mountpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mountpoint as text
%        str2double(get(hObject,'String')) returns contents of mountpoint as a double


% --- Executes during object creation, after setting all properties.
function mountpoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mountpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in dyn_mod.
function dyn_mod_Callback(hObject, eventdata, handles)
% hObject    handle to dyn_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns dyn_mod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dyn_mod
global order

contents = cellstr(get(hObject,'String'));
if (strcmp(contents{get(hObject,'Value')},'Static'))
    order = 1;
    set(handles.std_X, 'Enable', 'off');
    set(handles.std_Y, 'Enable', 'off');
    set(handles.std_Z, 'Enable', 'off');
    set(handles.text_std_X, 'Enable', 'off');
    set(handles.text_std_Y, 'Enable', 'off');
    set(handles.text_std_Z, 'Enable', 'off');
    set(handles.text_std_X_unit, 'Enable', 'off');
    set(handles.text_std_Y_unit, 'Enable', 'off');
    set(handles.text_std_Z_unit, 'Enable', 'off');
else
    mode = cellstr(get(handles.nav_mon,'String'));
    if (strcmp(mode{get(handles.nav_mon,'Value')},'Navigation'))
        set(handles.std_X, 'Enable', 'on');
        set(handles.std_Y, 'Enable', 'on');
        set(handles.std_Z, 'Enable', 'on');
        set(handles.text_std_X, 'Enable', 'on');
        set(handles.text_std_Y, 'Enable', 'on');
        set(handles.text_std_Z, 'Enable', 'on');
        set(handles.text_std_X_unit, 'Enable', 'on');
        set(handles.text_std_Y_unit, 'Enable', 'on');
        set(handles.text_std_Z_unit, 'Enable', 'on');
    end
    if (strcmp(contents{get(hObject,'Value')},'Const. acceleration'))
        order = 3;
    elseif (strcmp(contents{get(hObject,'Value')},'Const. velocity'))
        order = 2;
    else
        order = 1;
    end
end


% --- Executes during object creation, after setting all properties.
function dyn_mod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dyn_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function approx_lat_Callback(hObject, eventdata, handles)
% hObject    handle to approx_lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of approx_lat as text
%        str2double(get(hObject,'String')) returns contents of approx_lat as a double


% --- Executes during object creation, after setting all properties.
function approx_lat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to approx_lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function approx_lon_Callback(hObject, eventdata, handles)
% hObject    handle to approx_lon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of approx_lon as text
%        str2double(get(hObject,'String')) returns contents of approx_lon as a double


% --- Executes during object creation, after setting all properties.
function approx_lon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to approx_lon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function approx_h_Callback(hObject, eventdata, handles)
% hObject    handle to approx_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of approx_h as text
%        str2double(get(hObject,'String')) returns contents of approx_h as a double


% --- Executes during object creation, after setting all properties.
function approx_h_CreateFcn(hObject, eventdata, handles)
% hObject    handle to approx_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function antenna_h_Callback(hObject, eventdata, handles)
% hObject    handle to antenna_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of antenna_h as text
%        str2double(get(hObject,'String')) returns contents of antenna_h as a double


% --- Executes during object creation, after setting all properties.
function antenna_h_CreateFcn(hObject, eventdata, handles)
% hObject    handle to antenna_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in go_button.
function go_button_Callback(hObject, eventdata, handles)
% hObject    handle to go_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%master station coordinates
crs_contents = cellstr(get(handles.crs,'String'));
master_X = str2double(get(handles.master_X,'String'));
master_Y = str2double(get(handles.master_Y,'String'));
master_Z = str2double(get(handles.master_Z,'String'));
master_lat = str2double(get(handles.master_lat,'String'));
master_lon = str2double(get(handles.master_lon,'String'));
master_h = str2double(get(handles.master_h,'String'));
%KF parameters
std_init = str2double(get(handles.std_init,'String'));
std_X = str2double(get(handles.std_X,'String'));
std_Y = str2double(get(handles.std_Y,'String'));
std_Z = str2double(get(handles.std_Z,'String'));
std_vel = str2double(get(handles.std_vel,'String'));
std_code = str2double(get(handles.std_code,'String'));
if (get(handles.toggle_std_phase,'Value'))
    std_phase = str2double(get(handles.std_phase,'String'));
else
    std_phase = 1e30;
end
if (get(handles.toggle_std_dtm,'Value'))
    std_dtm = str2double(get(handles.std_dtm,'String'));
else
    std_dtm = 1e30;
end
min_nsat = str2double(get(handles.min_sat,'String'));
cutoff = str2double(get(handles.cut_off,'String'));
snr_threshold = str2double(get(handles.snr_thres,'String'));
cs_threshold = str2double(get(handles.cs_thresh,'String'));
antenna_h = str2double(get(handles.antenna_h,'String'));
contents_dyn_mod = cellstr(get(handles.dyn_mod,'String'));
flag_stopGOstop = get(handles.stopGOstop,'Value');
%input files
filerootIN = get(handles.gogps_data_input,'String');
filerootOUT = [get(handles.gogps_data_output,'String') '\' get(handles.gogps_data_output_prefix,'String')];
filerootIN(filerootIN == '\') = '/';
filerootOUT(filerootOUT == '\') = '/';
filename_R_obs = get(handles.RINEX_rover_obs,'String');
filename_M_obs = get(handles.RINEX_master_obs,'String');
filename_nav = get(handles.RINEX_nav,'String');
ref_path = get(handles.ref_path, 'Value');
filename_ref = get(handles.ref_path_input,'String');
dtm_dir = get(handles.dtm_path,'String');
%serial communication
% global COMportR
contents = cellstr(get(handles.com_select_0,'String'));
COMportR0 = contents{get(handles.com_select_0,'Value')};
% contents = cellstr(get(handles.com_select_1,'String'));
% COMportR1 = contents{get(handles.com_select_1,'Value')};
% contents = cellstr(get(handles.com_select_2,'String'));
% COMportR2 = contents{get(handles.com_select_2,'Value')};
% contents = cellstr(get(handles.com_select_3,'String'));
% COMportR3 = contents{get(handles.com_select_3,'Value')};
%TCPIP / NTRIP
flag_NTRIP = get(handles.use_ntrip,'Value');
master_ip = get(handles.IP_address,'String');
master_port = str2double(get(handles.port,'String'));
ntrip_mountpoint = get(handles.mountpoint,'String');
%functioning mode
mode = select_mode(handles);

ready = 1;

%check if everything is OK before starting
if isempty(dir(get(handles.gogps_data_output,'String')))
    msgbox('Output folder does not exist. Please browse to an existing folder.'); ready = 0;
end

if ~isempty(dir(get(handles.gogps_data_output,'String')))
    i = 1;
    j = length(filerootOUT);
    while (~isempty(dir([filerootOUT '_????_rover.bin'])) | ...
            ~isempty(dir([filerootOUT '_master*.bin'])) | ...
            ~isempty(dir([filerootOUT '_????_obs*.bin'])) | ...
            ~isempty(dir([filerootOUT '_????_eph*.bin'])) | ...
            ~isempty(dir([filerootOUT '_????_dyn*.bin'])) | ...
            ~isempty(dir([filerootOUT '_sat*.bin'])) | ...
            ~isempty(dir([filerootOUT '_kal*.bin'])) | ...
            ~isempty(dir([filerootOUT '_dt*.bin'])) | ...
            ~isempty(dir([filerootOUT '_conf*.bin'])) | ...
            ~isempty(dir([filerootOUT '_dop*.bin'])) | ...
            ~isempty(dir([filerootOUT '_ECEF*.txt'])) | ...
            ~isempty(dir([filerootOUT '_geod*.txt'])) | ...
            ~isempty(dir([filerootOUT '_plan*.txt'])) | ...
            ~isempty(dir([filerootOUT '_????_NMEA*.txt'])) | ...
            ~isempty(dir([filerootOUT '.kml'])) )
        
        if (i == 100)
            msgbox('Automatic numbering of output file has reached the maximum (99). Please provide a different output prefix, select a different output folder or remove some files from the current output folder.');
            ready = 0;
            break
        end
        
        filerootOUT(j+1:j+3) = ['_' num2str(i,'%02d')];
        i = i + 1;
    end
end
if (mode < 10) %if post-processing
    if (get(handles.file_type, 'SelectedObject') == handles.gogps_data & (isempty(dir([filerootIN '_obs*.bin'])) | isempty(dir([filerootIN '_eph*.bin']))))
        msgbox('Input goGPS binary files not found (both *_obs_* and *_eph_* files are needed).'); ready = 0;
    elseif (get(handles.file_type, 'SelectedObject') == handles.rinex_files & isempty(dir(filename_R_obs)))
        msgbox('Input rover observation RINEX file not found.'); ready = 0;
    elseif (mod(mode,2) & get(handles.file_type, 'SelectedObject') == handles.rinex_files & isempty(dir(filename_M_obs)))
        msgbox('Input master observation RINEX file not found.'); ready = 0;
    elseif (get(handles.file_type, 'SelectedObject') == handles.rinex_files & isempty(dir(filename_nav)))
        msgbox('Input navigation RINEX file not found.'); ready = 0;
    end
end

if (mode < 11) %if not rover and/or master monitor
    if (ref_path & isempty(dir(filename_ref)))
        msgbox('Input reference path file not found.'); ready = 0;
    elseif ((mode == 1 | mode == 2) & get(handles.toggle_std_dtm,'Value'))
        try
            load([dtm_dir '/tiles/tile_header'], 'tile_header');
            load([dtm_dir '/tiles/tile_georef'], 'tile_georef');
        catch
            msgbox('DTM directory does not contain valid data.'); ready = 0;
        end
    elseif (mod(mode,2) & strcmp(get(handles.crs,'Enable'),'on') & ...
            strcmp(crs_contents{get(handles.crs,'Value')},'ECEF (X,Y,Z)') & ...
            ((isnan(master_X) | isnan(master_Y) | isnan(master_Z)) | ...
            (master_X == 0 & master_Y == 0 & master_Z == 0)))
        msgbox('Please provide valid values for the master station ECEF coordinates.'); ready = 0;
    elseif (mod(mode,2) & strcmp(get(handles.crs,'Enable'),'on') & ...
            strcmp(crs_contents{get(handles.crs,'Value')},'Geodetic coord. (lat,lon,h)') & ...
            ((isnan(master_lat) | isnan(master_lon) | isnan(master_h)) | ...
            (abs(master_lat) > 90 | abs(master_lon) > 180)))
        msgbox('Please provide valid values for the master station geodetic coordinates.'); ready = 0;
    elseif (isnan(std_init) | std_init < 0)
        msgbox('Please provide a valid value for the initial state error standard deviation.'); ready = 0;
    elseif (isnan(std_X) | std_X < 0)
        msgbox('Please provide a valid value for the East coordinate error standard deviation.'); ready = 0;
    elseif (isnan(std_Y) | std_Y < 0)
        msgbox('Please provide a valid value for the North coordinate error standard deviation.'); ready = 0;
    elseif (isnan(std_Z) | std_Z < 0)
        msgbox('Please provide a valid value for the Up coordinate error standard deviation.'); ready = 0;
    elseif (isnan(std_vel) | std_vel < 0)
        msgbox('Please provide a valid value for the velocity coordinate error standard deviation.'); ready = 0;
    elseif (isnan(std_code) | std_code < 0)
        msgbox('Please provide a valid value for the code error standard deviation.'); ready = 0;
    elseif (isnan(std_phase) | std_phase < 0)
        msgbox('Please provide a valid value for the phase error standard deviation.'); ready = 0;
    elseif (isnan(std_dtm) | std_dtm < 0)
        msgbox('Please provide a valid value for the DTM error standard deviation.'); ready = 0;
    elseif (isnan(min_nsat) | min_nsat < 0)
        msgbox('Please provide a valid value for the minimum number of satellites.'); ready = 0;
    elseif (isnan(cutoff) | cutoff < 0 | cutoff > 90)
        msgbox('Please provide a valid value for the cutoff (between 0 and 90 degrees).'); ready = 0;
    elseif (isnan(snr_threshold) | snr_threshold < 0 | snr_threshold > 60)
        msgbox('Please provide a valid value for the SNR threshold (between 0 and 60 dB).'); ready = 0;
    elseif (isnan(cs_threshold) | cs_threshold < 0)
        msgbox('Please provide a valid value for the cycle slip threshold.'); ready = 0;
    elseif (isnan(antenna_h) | antenna_h < 0)
        msgbox('Please provide a valid value for the antenna height.'); ready = 0;
    end
end

%check if the dataset was surveyed with a variable dynamic model
d = dir([filerootIN '_dyn_00.bin']);
if (mode < 10 & (flag_stopGOstop | strcmp(contents_dyn_mod{get(handles.dyn_mod,'Value')},'Variable')) & isempty(d))
    msgbox('The selected dataset was not surveyed with a variable dynamic model: please select another dynamic model.'); ready = 0;
end

if (mode == 11 | mode == 12 | mode == 14) %if a COM connection to the rover is required
    if(strcmp(COMportR0, 'NA'))
        msgbox('Please select an existing COM port.'); ready = 0;
    end
end

if (mode == 11 | mode == 13 | mode == 14) %if a TCP/IP connection to the master is required
    if (isempty(master_ip))
        msgbox('Please provide an IP address for the connection to the master.'); ready = 0;
    elseif (isnan(master_port) | master_port < 0 | master_port > 65535)
        msgbox('Please provide a valid port number for the connection to the master (between 0 and 65535).'); ready = 0;
    end
    if (flag_NTRIP) %if a NTRIP connection is required
        if (isempty(ntrip_mountpoint))
            msgbox('Please provide a mountpoint for the NTRIP connection.'); ready = 0;
        end
    end
end

if (ready)
    saveState(handles,'../data/settings/last_settings.mat');
    uiresume(handles.main_panel);
end

% --- Executes on button press in load_button.
function load_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%allow the user to choose which settings to load
[filename, pathname] = uigetfile('*.mat', 'Choose file with saved settings','../data/settings');

if pathname == 0 %if the user pressed cancelled, then we exit this callback
    return
end

%construct the path name of the file to be loaded
loadDataName = fullfile(pathname,filename);

%load the settings, which creates a new gui
loadState(handles, loadDataName);

% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%allow the user to specify where to save the settings file
[filename,pathname] = uiputfile('*.mat','Save your GUI settings','../data/settings');

if pathname == 0 %if the user pressed cancelled, then we exit this callback
    return
end
%construct the path name of the save location
saveDataName = fullfile(pathname,filename);

%saves the gui data
saveState(handles, saveDataName);

% --- Executes when selected object is changed in weight_select.
function weight_select_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in weight_select
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in com_select_0.
function com_select_0_Callback(hObject, eventdata, handles)
% hObject    handle to com_select_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns com_select_0 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from com_select_0
% contents1 = cellstr(get(handles.com_select_1,'String'));
% contents2 = cellstr(get(handles.com_select_2,'String'));
% contents3 = cellstr(get(handles.com_select_3,'String'));
% serialInfo = instrhwinfo('serial');
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents1{get(handles.com_select_1,'Value')})) = [];
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents2{get(handles.com_select_2,'Value')})) = [];
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents3{get(handles.com_select_3,'Value')})) = [];
% set(hObject, 'String', serialInfo.AvailableSerialPorts);

% contents  = cellstr(get(hObject,'String'));
% COMportR0 = contents{get(hObject,'Value')};



% --- Executes during object creation, after setting all properties.
function com_select_0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to com_select_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    serialInfo = instrhwinfo('serial');
catch
end
if (~isempty(serialInfo.AvailableSerialPorts))
    set(hObject, 'String', serialInfo.AvailableSerialPorts);
else
    notAvailable = cell(1);
    notAvailable{1} = 'NA';
    set(hObject, 'String', notAvailable);
end

% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)


% --- Executes on button press in no_skyplot_snr.
function no_skyplot_snr_Callback(hObject, eventdata, handles)
% hObject    handle to no_skyplot_snr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of no_skyplot_snr


% --- Executes on key press with focus on password and none of its controls.
function password_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to password (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
show = get(handles.show_password,'Value');
password = get(hObject,'Userdata');
mask = get(hObject,'String');
key = eventdata.Key;

switch key
    case 'backspace'
        password = password(1:end-1); % Delete the last character in the password
        if(show)
            set(hObject,'String',password)% Set the text in the password edit box to the password
        else
            mask = mask(1:end-1); % Delete the last asterisk
            set(hObject,'String',mask)% Set the text in the password edit box to the asterisk string
        end
    case 'delete'
        password = password(2:end); % Delete the first character in the password
        if(show)
            set(hObject,'String',password)% Set the text in the password edit box to the password
        else
            mask = mask(2:end); % Delete the first asterisk
            set(hObject,'String',mask)% Set the text in the password edit box to the asterisk string
        end
    otherwise
        % If pressed key produces a printable character
        if (uint8(eventdata.Character) > 32)
            password = [password eventdata.Character]; % Add the typed character to the password
            pause(0.001)%to avoid unwanted character output before the cursor
            if(show)
                set(hObject,'String',password)% Set the text in the password edit box to the password
            else
                mask = [mask '*']; % Add an asterisk
                set(hObject,'String',mask)% Set the text in the password edit box to the asterisk string
            end
        end
end

set(hObject,'Userdata',password) % Store the password in its current state


% --- Executes on button press in show_password.
function show_password_Callback(hObject, eventdata, handles)
% hObject    handle to show_password (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of show_password
if get(hObject,'Value')
    password = get(handles.password,'Userdata');
    set(handles.password,'String',password);
else
    SizePass = size(get(handles.password,'Userdata')); % Find the number of asterisks
    if SizePass(2) > 0
        asterisk(1,1:SizePass(2)) = '*'; % Create a string of asterisks the same size as the password
        set(handles.password,'String',asterisk) % Set the text in the password edit box to the asterisk string
    else
        set(handles.password,'String','')
    end
end


% --------------------------------------------------------------------
function menu_tools_Callback(hObject, eventdata, handles)
% hObject    handle to menu_tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_about_Callback(hObject, eventdata, handles)
% hObject    handle to menu_about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function decode_streams_unix_Callback(hObject, eventdata, handles)
% hObject    handle to decode_streams_unix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_decode_stream_unix;


% --------------------------------------------------------------------
function merge_goGPS_bin_unix_Callback(hObject, eventdata, handles)
% hObject    handle to merge_goGPS_bin_unix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_merge_goGPSbin_unix;


% --------------------------------------------------------------------
function RINEX2goGPSbin_unix_Callback(hObject, eventdata, handles)
% hObject    handle to RINEX2goGPSbin_unix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_RINEX2goGPSbin_unix;


% --- Executes on button press in plotproc.
function plotproc_Callback(hObject, eventdata, handles)
% hObject    handle to plotproc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotproc
if (get(hObject,'Value'))
    set(handles.no_skyplot_snr, 'Enable', 'on');
    set(handles.google_earth, 'Enable', 'on');
    set(handles.err_ellipse, 'Enable', 'on');
    set(handles.plot_master, 'Enable', 'on');
    check_mode = cellstr(get(handles.mode,'String'));
    check_phase = cellstr(get(handles.code_dd_sa,'String'));
    if (~strcmp(check_mode{get(handles.mode,'Value')},'Real-time') & (strcmp(check_phase{get(handles.code_dd_sa,'Value')}, ...
        'Code and phase double difference') | strcmp(check_phase{get(handles.code_dd_sa,'Value')},'Code and phase stand-alone')))
        set(handles.plot_amb, 'Enable', 'on');
        plot_amb_Callback(handles.plot_amb, [], handles);
    end
else
    set(handles.no_skyplot_snr, 'Enable', 'off');
    set(handles.google_earth, 'Enable', 'off');
    set(handles.err_ellipse, 'Enable', 'off');
    set(handles.plot_master, 'Enable', 'off');
    set(handles.plot_amb, 'Enable', 'off');
end


% --- Executes on button press in stopGOstop.
function stopGOstop_Callback(hObject, eventdata, handles)
% hObject    handle to stopGOstop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of stopGOstop
check_nav = cellstr(get(handles.nav_mon,'String'));
check_KF = cellstr(get(handles.kalman_ls,'String'));
if (get(hObject,'Value') | strcmp(check_nav{get(handles.nav_mon,'Value')},'Master monitor') | ~strcmp(check_KF{get(handles.kalman_ls,'Value')},'Kalman filter'))
    set(handles.dyn_mod, 'Enable', 'off');
else
    set(handles.dyn_mod, 'Enable', 'on');
end


% --------------------------------------------------------------------
function polyline_unix_Callback(hObject, eventdata, handles)
% hObject    handle to polyline_unix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_polyline_simplification_unix;


% --------------------------------------------------------------------
function help_Callback(hObject, eventdata, handles)
% hObject    handle to help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function about_unix_Callback(hObject, eventdata, handles)
% hObject    handle to about_unix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_about_unix;


% --- Executes on selection change in protocol_select_0.
function protocol_select_0_Callback(hObject, eventdata, handles)
% hObject    handle to protocol_select_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns protocol_select_0 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from protocol_select_0


% --- Executes during object creation, after setting all properties.
function protocol_select_0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to protocol_select_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function mode = select_mode(handles)
contents_mode = cellstr(get(handles.mode,'String'));
contents_nav_mon = cellstr(get(handles.nav_mon,'String'));
contents_kalman_ls = cellstr(get(handles.kalman_ls,'String'));
contents_code_dd_sa = cellstr(get(handles.code_dd_sa,'String'));
if (strcmp(contents_mode{get(handles.mode,'Value')},'Post-processing'))
    if (strcmp(contents_kalman_ls{get(handles.kalman_ls,'Value')},'Kalman filter'))
        if (strcmp(contents_code_dd_sa{get(handles.code_dd_sa,'Value')},'Code and phase double difference'))
            mode = 1;
        elseif (strcmp(contents_code_dd_sa{get(handles.code_dd_sa,'Value')},'Code and phase stand-alone'))
            mode = 2;
        elseif (strcmp(contents_code_dd_sa{get(handles.code_dd_sa,'Value')},'Code double difference'))
            mode = 5;
        else %Code stand-alone
            mode = 6;
        end
    else %Least squares
        if (strcmp(contents_code_dd_sa{get(handles.code_dd_sa,'Value')},'Code double difference'))
            mode = 3;
        else %Code stand-alone
            mode = 4;
        end
    end
else %Real-time
    if (strcmp(contents_nav_mon{get(handles.nav_mon,'Value')},'Navigation'))
        mode = 11;
    elseif (strcmp(contents_nav_mon{get(handles.nav_mon,'Value')},'Rover monitor'))
        mode = 12;
    elseif (strcmp(contents_nav_mon{get(handles.nav_mon,'Value')},'Master monitor'))
        mode = 13;
    else %Rover and Master monitor
        mode = 14;
    end
end


% --- Executes on button press in flag_doppler.
function flag_doppler_Callback(hObject, eventdata, handles)
% hObject    handle to flag_doppler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of flag_doppler


% --- Executes on selection change in amb_select.
function amb_select_Callback(hObject, eventdata, handles)
% hObject    handle to amb_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns amb_select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from amb_select
global amb_restart_method
contents = cellstr(get(hObject,'String'));
selection = contents{get(hObject,'Value')};
if (strcmp(selection, 'Observed code - phase difference'))
    amb_restart_method = 0;
elseif (strcmp(selection, 'Kalman-predicted code - phase difference'))
    amb_restart_method = 1;
else
    amb_restart_method = 2;
end

% --- Executes during object creation, after setting all properties.
function amb_select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to amb_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in com_select_1.
function com_select_1_Callback(hObject, eventdata, handles)
% hObject    handle to com_select_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns com_select_1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from com_select_1

% contents0 = cellstr(get(handles.com_select_0,'String'));
% contents2 = cellstr(get(handles.com_select_2,'String'));
% contents3 = cellstr(get(handles.com_select_3,'String'));
% serialInfo = instrhwinfo('serial');
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents0{get(handles.com_select_0,'Value')})) = [];
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents2{get(handles.com_select_2,'Value')})) = [];
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents3{get(handles.com_select_3,'Value')})) = [];
% set(hObject, 'String', serialInfo.AvailableSerialPorts);

% contents1  = cellstr(get(hObject,'String'));
% COMportR1 = contents1{get(hObject,'Value')};

% --- Executes during object creation, after setting all properties.
function com_select_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to com_select_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    serialInfo = instrhwinfo('serial');
catch
end
if (~isempty(serialInfo.AvailableSerialPorts))
%     serialInfo.AvailableSerialPorts(1) = [];
%     serialInfo.AvailableSerialPorts(3) = [];
%     serialInfo.AvailableSerialPorts(4) = [];
    set(hObject, 'String', serialInfo.AvailableSerialPorts);
end
if (isempty(serialInfo.AvailableSerialPorts))
    notAvailable = cell(1);
    notAvailable{1} = 'NA';
    set(hObject, 'String', notAvailable);
end


% --- Executes on selection change in protocol_select_1.
function protocol_select_1_Callback(hObject, eventdata, handles)
% hObject    handle to protocol_select_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns protocol_select_1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from protocol_select_1


% --- Executes during object creation, after setting all properties.
function protocol_select_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to protocol_select_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in com_select_2.
function com_select_2_Callback(hObject, eventdata, handles)
% hObject    handle to com_select_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns com_select_2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from com_select_2
% contents0 = cellstr(get(handles.com_select_0,'String'));
% contents1 = cellstr(get(handles.com_select_1,'String'));
% contents3 = cellstr(get(handles.com_select_3,'String'));
% serialInfo = instrhwinfo('serial');
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents0{get(handles.com_select_0,'Value')})) = [];
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents1{get(handles.com_select_1,'Value')})) = [];
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents3{get(handles.com_select_3,'Value')})) = [];
% set(hObject, 'String', serialInfo.AvailableSerialPorts);

% contents  = cellstr(get(hObject,'String'));
% COMportR2 = contents{get(hObject,'Value')};

% --- Executes during object creation, after setting all properties.
function com_select_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to com_select_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    serialInfo = instrhwinfo('serial');
catch
end
if (~isempty(serialInfo.AvailableSerialPorts))
%     serialInfo.AvailableSerialPorts(1) = [];
%     serialInfo.AvailableSerialPorts(2) = [];
%     serialInfo.AvailableSerialPorts(4) = [];
    set(hObject, 'String', serialInfo.AvailableSerialPorts);
end
if (isempty(serialInfo.AvailableSerialPorts))
    notAvailable = cell(1);
    notAvailable{1} = 'NA';
    set(hObject, 'String', notAvailable);
end


% --- Executes on selection change in protocol_select_2.
function protocol_select_2_Callback(hObject, eventdata, handles)
% hObject    handle to protocol_select_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns protocol_select_2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from protocol_select_2


% --- Executes during object creation, after setting all properties.
function protocol_select_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to protocol_select_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in com_select_3.
function com_select_3_Callback(hObject, eventdata, handles)
% hObject    handle to com_select_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns com_select_3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from com_select_3
% contents0 = cellstr(get(handles.com_select_0,'String'));
% contents1 = cellstr(get(handles.com_select_1,'String'));
% contents2 = cellstr(get(handles.com_select_2,'String'));
% serialInfo = instrhwinfo('serial');
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents0{get(handles.com_select_0,'Value')})) = [];
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents1{get(handles.com_select_1,'Value')})) = [];
% serialInfo.AvailableSerialPorts(strcmp(serialInfo.AvailableSerialPorts,contents2{get(handles.com_select_2,'Value')})) = [];
% set(hObject, 'String', serialInfo.AvailableSerialPorts);

% contents  = cellstr(get(hObject,'String'));
% COMportR3 = contents{get(hObject,'Value')};


% --- Executes during object creation, after setting all properties.
function com_select_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to com_select_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    serialInfo = instrhwinfo('serial');
catch
end
if (~isempty(serialInfo.AvailableSerialPorts))
%     serialInfo.AvailableSerialPorts(1) = [];
%     serialInfo.AvailableSerialPorts(2) = [];
%     serialInfo.AvailableSerialPorts(3) = [];
    set(hObject, 'String', serialInfo.AvailableSerialPorts);
end
if (isempty(serialInfo.AvailableSerialPorts))
    notAvailable = cell(1);
    notAvailable{1} = 'NA';
    set(hObject, 'String', notAvailable);
end


% --- Executes on selection change in protocol_select_3.
function protocol_select_3_Callback(hObject, eventdata, handles)
% hObject    handle to protocol_select_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns protocol_select_3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from protocol_select_3


% --- Executes during object creation, after setting all properties.
function protocol_select_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to protocol_select_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in num_receivers.
function num_receivers_Callback(hObject, eventdata, handles)
% hObject    handle to num_receivers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns num_receivers contents as cell array
%        contents{get(hObject,'Value')} returns selected item from num_receivers
contents = cellstr(get(hObject,'String'));
check_mode = cellstr(get(handles.mode,'String'));
contents_nav_mon = cellstr(get(handles.nav_mon,'String'));
if (strcmp(check_mode{get(handles.mode,'Value')},'Real-time') & strcmp(contents_nav_mon{get(handles.nav_mon,'Value')},'Rover monitor'))
    set(handles.text_num_receivers, 'Enable', 'on');
    set(handles.num_receivers, 'Enable', 'on');
    if (size(contents,1) >= get(hObject,'Value'))
        if (strcmp(contents{get(hObject,'Value')},'1'))
            set(handles.com_select_0, 'Enable', 'on');  set(handles.protocol_select_0, 'Enable', 'on');
            set(handles.com_select_1, 'Enable', 'off'); set(handles.protocol_select_1, 'Enable', 'off');
            set(handles.com_select_2, 'Enable', 'off'); set(handles.protocol_select_2, 'Enable', 'off');
            set(handles.com_select_3, 'Enable', 'off'); set(handles.protocol_select_3, 'Enable', 'off');
        elseif (strcmp(contents{get(hObject,'Value')},'2'))
            set(handles.com_select_0, 'Enable', 'on');  set(handles.protocol_select_0, 'Enable', 'on');
            set(handles.com_select_1, 'Enable', 'on');  set(handles.protocol_select_1, 'Enable', 'on');
            set(handles.com_select_2, 'Enable', 'off'); set(handles.protocol_select_2, 'Enable', 'off');
            set(handles.com_select_3, 'Enable', 'off'); set(handles.protocol_select_3, 'Enable', 'off');
        elseif (strcmp(contents{get(hObject,'Value')},'3'))
            set(handles.com_select_0, 'Enable', 'on');  set(handles.protocol_select_0, 'Enable', 'on');
            set(handles.com_select_1, 'Enable', 'on');  set(handles.protocol_select_1, 'Enable', 'on');
            set(handles.com_select_2, 'Enable', 'on');  set(handles.protocol_select_2, 'Enable', 'on');
            set(handles.com_select_3, 'Enable', 'off'); set(handles.protocol_select_3, 'Enable', 'off');
        elseif (strcmp(contents{get(hObject,'Value')},'4'))
            set(handles.com_select_0, 'Enable', 'on');  set(handles.protocol_select_0, 'Enable', 'on');
            set(handles.com_select_1, 'Enable', 'on');  set(handles.protocol_select_1, 'Enable', 'on');
            set(handles.com_select_2, 'Enable', 'on');  set(handles.protocol_select_2, 'Enable', 'on');
            set(handles.com_select_3, 'Enable', 'on');  set(handles.protocol_select_3, 'Enable', 'on');
        end
    else
        set(hObject,'Value',1);
        set(handles.com_select_0, 'Enable', 'on');  set(handles.protocol_select_0, 'Enable', 'on');
        set(handles.com_select_1, 'Enable', 'off'); set(handles.protocol_select_1, 'Enable', 'off');
        set(handles.com_select_2, 'Enable', 'off'); set(handles.protocol_select_2, 'Enable', 'off');
        set(handles.com_select_3, 'Enable', 'off'); set(handles.protocol_select_3, 'Enable', 'off');
    end
elseif (strcmp(check_mode{get(handles.mode,'Value')},'Real-time') & strcmp(contents_nav_mon{get(handles.nav_mon,'Value')},'Master monitor'))
    set(handles.text_num_receivers, 'Enable', 'off');
    set(handles.num_receivers, 'Enable', 'off');
    set(handles.com_select_0, 'Enable', 'off');
    set(handles.com_select_1, 'Enable', 'off');
    set(handles.com_select_2, 'Enable', 'off');
    set(handles.com_select_3, 'Enable', 'off');
    set(handles.protocol_select_0, 'Enable', 'off');
    set(handles.protocol_select_1, 'Enable', 'off');
    set(handles.protocol_select_2, 'Enable', 'off');
    set(handles.protocol_select_3, 'Enable', 'off');
else
    set(hObject,'Value',1);
    set(handles.com_select_0, 'Enable', 'on');  set(handles.protocol_select_0, 'Enable', 'on');
    set(handles.com_select_1, 'Enable', 'off'); set(handles.protocol_select_1, 'Enable', 'off');
    set(handles.com_select_2, 'Enable', 'off'); set(handles.protocol_select_2, 'Enable', 'off');
    set(handles.com_select_3, 'Enable', 'off'); set(handles.protocol_select_3, 'Enable', 'off');
end
    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------


% --- Executes during object creation, after setting all properties.
function num_receivers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_receivers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

contents = {'1';'2';'3';'4'};
serialInfo = instrhwinfo('serial');
num_ports = size(serialInfo.AvailableSerialPorts,1);
if num_ports == 0
    set(hObject,'String','1');
elseif num_ports <= size(contents,1);
    set(hObject,'String',contents(1:num_ports));
else
    set(hObject,'String',contents);
end
