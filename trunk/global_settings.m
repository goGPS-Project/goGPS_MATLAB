% SYNTAX:
%   global_settings;
%
% DESCRIPTION:
%   User-defined global settings.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 pre-alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Media Center, Osaka City University, Japan
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

%-------------------------------------------------------------------------------
% INPUT/OUTPUT FILENAME PREFIX
%-------------------------------------------------------------------------------

filerootIN  = '../data/data_goGPS/ultresc';
filerootOUT = '../data/out';

i = 1;
j = length(filerootOUT);
while (~isempty(dir([filerootOUT '_rover*.bin'])) | ...
       ~isempty(dir([filerootOUT '_master*.bin'])) | ...
       ~isempty(dir([filerootOUT '_obs*.bin'])) | ...
       ~isempty(dir([filerootOUT '_eph*.bin'])) | ...
       ~isempty(dir([filerootOUT '_sat*.bin'])) | ...
       ~isempty(dir([filerootOUT '_kal*.bin'])) | ...
       ~isempty(dir([filerootOUT '_dt*.bin'])) | ...
       ~isempty(dir([filerootOUT '_conf*.bin'])) | ...
       ~isempty(dir([filerootOUT '_geod*.txt'])) | ...
       ~isempty(dir([filerootOUT '_plan*.txt'])) | ...
       ~isempty(dir([filerootOUT '_NMEA*.txt'])) )
   
   filerootOUT(j+1:j+3) = ['_' num2str(i,'%02d')];
   i = i + 1;
end

%-------------------------------------------------------------------------------
% INPUT FILENAMES (RINEX data mode)
%-------------------------------------------------------------------------------

if (mode_data == 0)

    filename_R_obs    = '../data/data_RINEX/basket/perim2.08o';
    filename_R_nav    = '../data/data_RINEX/basket/COMO1190.08n';
    filename_M_obs    = '../data/data_RINEX/basket/COMO1190.08o';
    filename_M_nav    = '../data/data_RINEX/basket/COMO1190.08n';

    %---------------------------------------------------------------------------
    % FILE CHECK
    %---------------------------------------------------------------------------

    %check the file existence
    if (fopen(filename_R_obs,'r') == -1)
        error('Warning: ROVER observation file not found.');
    end

    if (fopen(filename_R_nav,'r') == -1)
        error('Warning: ROVER navigation file not found.');
    end

    if (fopen(filename_M_obs,'r') == -1)
        error('Warning: MASTER observation file not found.');
    end

    if (fopen(filename_M_nav,'r') == -1)
        error('Warning: MASTER navigation file not found.');
    end

    %close all temporarily opened files
    fclose('all');

end

%-------------------------------------------------------------------------------
% REFERENCE PATH FILENAME
%-------------------------------------------------------------------------------

filename_ref = '../data/data_RINEX/basket/refBASKET.mat';

%-------------------------------------------------------------------------------
% MASTER STATION POSITION
%-------------------------------------------------------------------------------

global XM YM ZM
global phiM lamM hM
global EST_M NORD_M
global NM MM RM

%Como permanent station (marker, January 2008):
%XM = 4398306.2420;
%YM =  704149.9120;
%ZM = 4550154.7080;

%Como permanent station (base of antenna, January 2008):
XM = 4398306.3892;
YM =  704149.9356;
ZM = 4550154.8609;

%Como permanent station antenna - TPSCR3_GGD CONE - phase center
% L1: EAST=+0.0008, NORTH=-0.0004, h=0.0615
% L2: EAST=+0.0003, NORTH=+0,      h=0.0948

%Master position in geodetic coordinates
%'a' and 'e' are global vars, but should be passed as parameters
[phiM, lamM, hM] = cart2geod(XM, YM, ZM);

%Master position in UTM coordinates (East, North, h)
[EST_M, NORD_M] = geod2plan(phiM, lamM);

NM = a / sqrt(1 - e^2 * (sin(phiM))^2);
MM = NM * (1 - e^2) / (1 - e^2 * (sin(phiM))^2);
RM = sqrt(NM*MM);

%-------------------------------------------------------------------------------
% KALMAN FILTER
%-------------------------------------------------------------------------------

global sigmaq0 sigmaq_velx sigmaq_vely sigmaq_velz sigmaq_vel
global sigmaq_cod1 sigmaq_cod2 sigmaq_ph sigmaq0_N sigmaq_dtm
global min_nsat cutoff cutoff_init snr_threshold weights order o1 o2 o3

%variance of initial state
sigmaq0 = 1e2;

%variance of velocity coordinates [m^2/s^2]
sigmaq_velx = 1e-1;
sigmaq_vely = 1e-1;
sigmaq_velz = 1e-1;
sigmaq_vel = 1e-0;

%variance of code observations [m^2]
% sigmaq_cod1 = 0.36;
sigmaq_cod1 = 9;
sigmaq_cod2 = 0.16;

%variance of phase observations [m^2]
%(maximize to obtain a code-based solution)
% sigmaq_ph = 0.000004;
sigmaq_ph = 0.001;
%sigmaq_ph = 0.001e30;

%variance of ambiguity combinations [cycles]
sigmaq0_N = 10;

%variance of DEM height [m^2]
%(maximize to disable DEM usage)
% sigmaq_dtm = 0.09;
sigmaq_dtm = 1e30;

%minimum number of satellites to be used in the Kalman filter
min_nsat = 2;

%cut-off [degrees]
cutoff = 15;

%initialization cut-off [degrees]
cutoff_init = 15;

%signal-to-noise ratio threshold [dB]
snr_threshold = 0;

%parameter used to select the weight mode for GPS observations
%          - weights=0: same weight for all the observations
%          - weights=1: weight based on satellite elevation
%          - weights=2: weight based on signal-to-noise ratio
%          - weights=3: weight based on combined elevation and signal-to-noise ratio
weights = 3;

%order of the dynamic model polynomial
order = 2;

%useful values to index matrices
o1 = order;
o2 = order*2;
o3 = order*3;

%-------------------------------------------------------------------------------
% RECEIVER
%-------------------------------------------------------------------------------

global h_antenna

%antenna height from the ground [m]
h_antenna = 0;
% h_antenna = 1.07;

%-------------------------------------------------------------------------------
% DTM (SET PATH AND LOAD PARAMETER FILES)
%-------------------------------------------------------------------------------

%parameters common to all DTM tiles
global tile_header

%parameters used to georeference every DTM tile
global tile_georef

%path to DTM folder containing DTM files
global dtm_dir

%folder containing DTM files
dtm_dir = '../data/dtm/20m';

fid = fopen([dtm_dir '/tiles/tile_header.mat'],'r');
if (fid ~= -1)
    load([dtm_dir '/tiles/tile_header'], 'tile_header');
    load([dtm_dir '/tiles/tile_georef'], 'tile_georef');
    fclose(fid);
else
    tile_header.nrows = 0;
    tile_header.ncols = 0;
    tile_header.cellsize = 0;
    tile_header.nodata = 0;
    tile_georef = zeros(1,1,4);
end

%-------------------------------------------------------------------------------
% MATLAB DISPLAY
%-------------------------------------------------------------------------------

global p_max pid
global window
global x_circle id_ellipse

%maximum number of drawn trajectory points
p_max = 200;

%trajectory point id
pid = zeros(p_max,1);

%dimension of the time windows (ambiguity plot)
window = 20;

%reference circle (for covariance plot)
try
    % if Statistics Toolbox is installed
    alpha = 0.05;                       % significance level
    r = sqrt(chi2inv(1-alpha,2));       % circle radius
catch
    % if Statistics Toolbox is not installed,
    % then alpha = 0.05, chi2inv(1-alpha,2) = 5.991
    r = sqrt(5.991);
end
theta = (0 : 2*pi/100 : 2*pi)';     % angle values
x_circle(:,1) = r * cos(theta);     % x-coordinate of the circle
x_circle(:,2) = r * sin(theta);     % y-coordinate of the circle

%error ellipse identifier
id_ellipse = [];

%-------------------------------------------------------------------------------
% COMMUNICATION INTERFACES
%-------------------------------------------------------------------------------

global COMportR
global master_ip master_port ntrip_user ntrip_pw ntrip_mountpoint
global server_delay
global nmea_init

if (mode == 11 | mode == 12) & flag_COM == 1
   %detect u-blox COM port
   COMportR = ublox_COM_find();

   if (isempty(COMportR))
        %Override rover data input port
        COMportR = 'COM8';
   end
else
    COMportR = 'COM8';
end

%MASTER/NTRIP connection parameters
master_ip = 'xxx.xxx.xxx.xxx';
master_port = 2101;

%NTRIP parameters
ntrip_user = 'uuuuuu';
ntrip_pw = 'ppppp';
ntrip_mountpoint = 'mmmmmmm';

%server waiting time (to check if packet transmission is finished)
% ( 0 < waiting time < 1, depending on server speed)
server_delay = 0.05;

%Initial NMEA sentence required by some NTRIP casters
nmea_init = NMEA_string_generator([XM YM ZM],5);

% %approximate coordinates to initialize NMEA string
% phiApp = 45.00;
% lamApp = 9.00;
% hApp = 200;
%
% [XApp,YApp,ZApp] = geod2cart (phiApp*pi/180, lamApp*pi/180, hApp, a, f);
%
% nmea_init = NMEA_string_generator([XApp YApp ZApp],0);

%-------------------------------------------------------------------------------
% INTERNET CONNECTION
%-------------------------------------------------------------------------------

global connection_delay

%waiting time for the Internet connection to be established
connection_delay = 5;

%-------------------------------------------------------------------------------
% GOOGLE EARTH
%-------------------------------------------------------------------------------

global GE_path GE_append
global link_filename kml_filename

%path to link to Google Earth executable
current_path = pwd;
current_path(current_path == '\') = '/';
GE_path = ['"' current_path '/../data/google_earth/googleearth.exe.lnk"'];

%KML file append(1) or re-write(0)
GE_append = 0;

%files used by Google Earth
link_filename = '../data/google_earth/link.kml';
kml_filename = '../data/google_earth/goGPS.kml';
