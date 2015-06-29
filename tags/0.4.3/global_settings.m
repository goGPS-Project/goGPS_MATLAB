% SYNTAX:
%   global_settings;
%
% DESCRIPTION:
%   User-defined global settings.

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

%-------------------------------------------------------------------------------
% INPUT/OUTPUT FILENAME PREFIX
%-------------------------------------------------------------------------------

folderIN  = '../data/data_goGPS';
prefixIN  = 'yamatogawa';

folderOUT = '../data/out/batch';
prefixOUT = 'mila_UBLX_14031';

filerootIN  = [folderIN '/' prefixIN];
filerootOUT = [folderOUT '/' prefixOUT];

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
       ~isempty(dir([filerootOUT '_dop*.bin'])) | ...
       ~isempty(dir([filerootOUT '_ECEF*.txt'])) | ...
       ~isempty(dir([filerootOUT '_geod*.txt'])) | ...
       ~isempty(dir([filerootOUT '_plan*.txt'])) | ...
       ~isempty(dir([filerootOUT '_NMEA*.txt'])) | ...
       ~isempty(dir([filerootOUT '_ublox_NMEA*.txt'])) | ...
       ~isempty(dir([filerootOUT '.kml'])) )

   filerootOUT(j+1:j+4) = ['_' num2str(i,'%03d')];
   i = i + 1;
end

%-------------------------------------------------------------------------------
% MASTER STATION POSITION
%-------------------------------------------------------------------------------

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

%set master station position manually
pos_M_man = [XM; YM; ZM];

%-------------------------------------------------------------------------------
% KALMAN FILTER
%-------------------------------------------------------------------------------

global sigmaq0 sigmaq_vE sigmaq_vN sigmaq_vU sigmaq_vel %#ok<*TLEV>
global sigmaq_cod1 sigmaq_cod2 sigmaq_ph sigmaq0_N sigmaq_dtm
global min_nsat cutoff snr_threshold cs_threshold weights snr_a snr_0 snr_1 snr_A order o1 o2 o3
global amb_restart_method

%variance of initial state
sigmaq0 = 1;

%variance of velocity coordinates [m^2/s^2]
sigmaq_vE = 0.5^2;
sigmaq_vN = 0.5^2;
sigmaq_vU = 0.1^2;
sigmaq_vel = 0.1^2;

%variance of code observations [m^2]
sigmaq_cod1 = 0.3^2;
sigmaq_cod2 = 0.4^2;

%variance of phase observations [m^2]
%(maximize to obtain a code-based solution)
sigmaq_ph = 0.003^2;
% sigmaq_ph = 0.001e30;

%variance of ambiguity combinations [cycles]
sigmaq0_N = 1000;

%variance of DEM height [m^2]
%(maximize to disable DEM usage)
% sigmaq_dtm = 0.09;
sigmaq_dtm = 1e30;

%minimum number of satellites to be used in the Kalman filter
min_nsat = 2;

%cut-off [degrees]
cutoff = 15;

%initialization cut-off [degrees]
% cutoff_init = 15;

%signal-to-noise ratio threshold [dB]
snr_threshold = 0;

%cycle slip threshold [cycles]
cs_threshold = 1;

%parameter used to select the weight mode for GPS observations
%          - weights=0: same weight for all the observations
%          - weights=1: weight based on satellite elevation (sin)
%          - weights=2: weight based on signal-to-noise ratio
%          - weights=3: weight based on combined elevation and signal-to-noise ratio
%          - weights=4: weight based on satellite elevation (exp)
weights = 1;

%weight function parameters
snr_a = 30;
snr_0 = 10;
snr_1 = 50;
snr_A = 30;

%order of the dynamic model polynomial
order = 1;

%useful values to index matrices
o1 = order;
o2 = order*2;
o3 = order*3;

%ambiguity restart method
amb_restart_method = 2;

%-------------------------------------------------------------------------------
% INTEGER AMBIGUITY RESOLUTION
%-------------------------------------------------------------------------------
global IAR_method P0 mu flag_auto_mu

%choose Integer Least Squares estimator
%IAR_method = 0; %ILS method with numeration in search (LAMBDA2)
IAR_method = 1; %ILS method with shrinking ellipsoid during search (LAMBDA3)
%IAR_method = 2; %ILS method with numeration in search (LAMBDA3)
%IAR_method = 3; %integer rounding method (LAMBDA3)
%IAR_method = 4; %integer bootstrapping method (LAMBDA3)
%IAR_method = 5; %Partial Ambiguity Resolution (PAR) (LAMBDA3)

%user defined fixed failure rate (for methods 1,2) or
%minimum required success rate (for method 5)
P0 = 0.001;

%user defined threshold for ratio test
mu = 0.5;

%flag for enabling the automatic determination of mu
flag_auto_mu = 1;

%flag for enabling the default value for P0
flag_default_P0 = 1;

%-------------------------------------------------------------------------------
% RECEIVER
%-------------------------------------------------------------------------------

global h_antenna

%antenna height from the ground [m]
h_antenna = 0;

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
dtm_dir = '../data/dtm';

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

%-------------------------------------------------------------------------------
% COMMUNICATION INTERFACES
%-------------------------------------------------------------------------------

global COMportR
global master_ip master_port ntrip_user ntrip_pw ntrip_mountpoint %#ok<NUSED>
global nmea_init

manCOMport = 'COM8';

COMportR = manCOMport;

% %MASTER/NTRIP connection parameters
% master_ip = 'xxx.xxx.xxx.xxx';
% master_port = 2101;
%
% %NTRIP parameters
% ntrip_user = 'uuuuuu';
% ntrip_pw = 'ppppp';
% ntrip_mountpoint = 'mmmmmmm';

%set approximate coordinates manually to initialize NMEA string
% phiApp = 34.5922;
% lamApp = 135.5059;
% hApp = 20;
% 
% [XApp,YApp,ZApp] = geod2cart (phiApp*pi/180, lamApp*pi/180, hApp, a_GPS, f_GPS);

XApp = -3749409.0399;
YApp = 3683775.1374;
ZApp = 3600727.7577;

%Initial NMEA sentence required by some NTRIP casters
nmea_init = NMEA_GGA_gen([XApp YApp ZApp],10);

%-------------------------------------------------------------------------------
% INI file
%-------------------------------------------------------------------------------
iniFile = './settings/Milano_daily_test_VRS_InputFiles.ini';

%initialize INI file reading
global goIni;
goIni = goIniReader;
goIni.setFileName(iniFile);

%extract user-defined settings from INI file
data_path = goIni.getData('Receivers','data_path');
file_name = goIni.getData('Receivers','file_name');
filename_R_obs = [data_path file_name];
data_path = goIni.getData('Master','data_path');
file_name = goIni.getData('Master','file_name');
filename_M_obs = [data_path file_name];
data_path = goIni.getData('Navigational','data_path');
file_name = goIni.getData('Navigational','file_name');
filename_nav = [data_path file_name];
flag_SP3 = goIni.getData('Navigational','isSP3');
if isempty(flag_SP3)
    if (strcmpi(filename_nav(end-3:end),'.sp3'))
        flag_SP3 = 1;
    else
        flag_SP3 = 0;
    end
end
data_path = goIni.getData('RefPath','data_path');
file_name = goIni.getData('RefPath','file_name');
filename_ref = [data_path file_name];
