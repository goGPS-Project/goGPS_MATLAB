% SYNTAX:
%   global_settings;
%
% DESCRIPTION:
%   User-defined global settings.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.2 alpha
%
% Copyright (C) 2009-2010 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Graduate School for Creative Cities, Osaka City University, Japan
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
folderOUT = '../data';

prefixIN  = 'yamatogawa';
prefixOUT = 'out';

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
    if (~exist(filename_R_obs,'file'))
        error('Warning: ROVER observation file not found.');
    end

    if (~exist(filename_R_nav,'file'))
        error('Warning: ROVER navigation file not found.');
    end

    if (~exist(filename_M_obs,'file'))
        error('Warning: MASTER observation file not found.');
    end

    if (~exist(filename_M_nav,'file'))
        error('Warning: MASTER navigation file not found.');
    end

end

%-------------------------------------------------------------------------------
% REFERENCE PATH FILENAME
%-------------------------------------------------------------------------------

filename_ref = '../data/data_RINEX/basket/refBASKET.mat';

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

global sigmaq0 sigmaq_velx sigmaq_vely sigmaq_velz sigmaq_vel
global sigmaq_cod1 sigmaq_cod2 sigmaq_ph sigmaq0_N sigmaq_dtm
global min_nsat cutoff snr_threshold cs_threshold weights snr_a snr_0 snr_1 snr_A order o1 o2 o3

%variance of initial state
sigmaq0 = 9;

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
cutoff = 0;

%initialization cut-off [degrees]
% cutoff_init = 15;

%signal-to-noise ratio threshold [dB]
snr_threshold = 0;

%cycle slip threshold [cycles]
cs_threshold = 10;

%parameter used to select the weight mode for GPS observations
%          - weights=0: same weight for all the observations
%          - weights=1: weight based on satellite elevation
%          - weights=2: weight based on signal-to-noise ratio
%          - weights=3: weight based on combined elevation and signal-to-noise ratio
weights = 2;

%weight function parameters
snr_a = 30;
snr_0 = 10;
snr_1 = 50;
snr_A = 30;

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
h_antenna = 1;

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

if (mode == 11 | mode == 12) & flag_COM == 1
   %detect u-blox COM port
   COMportR = ublox_COM_find();

   if (isempty(COMportR))
        %Override rover data input port
        COMportR = manCOMport;
   end
else
    COMportR = manCOMport;
end

% %MASTER/NTRIP connection parameters
% master_ip = 'xxx.xxx.xxx.xxx';
% master_port = 2101;
%
% %NTRIP parameters
% ntrip_user = 'uuuuuu';
% ntrip_pw = 'ppppp';
% ntrip_mountpoint = 'mmmmmmm';

%set approximate coordinates manually to initialize NMEA string
phiApp = 34.5922;
lamApp = 135.5059;
hApp = 20;

[XApp,YApp,ZApp] = geod2cart (phiApp*pi/180, lamApp*pi/180, hApp, a, f);

%Initial NMEA sentence required by some NTRIP casters
nmea_init = NMEA_GGA_gen([XApp YApp ZApp],10);

%-------------------------------------------------------------------------------
% AMBIGUITY ESTIMATION
%-------------------------------------------------------------------------------

global flag_LS_N_estim

%use least squares ambiguity estimation when new satellites
% are available and when cycle slips occur (0 = NO; 1 = YES)
%
% NOTE: LS amb. estimation is automatically switched off when the number of
% satellites with phase available is not sufficient (< 4 incl. pivot)
flag_LS_N_estim = 0;
