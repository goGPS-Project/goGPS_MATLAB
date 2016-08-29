function goGPS_batch(iniFile, mode, year_start, year_end, doy_start, doy_end, markerR, sessionR, extR, markerM, sessionM, extM, idN, sessionN, extN, folder) %#ok<INUSL>

% SYNTAX:
%   goGPS_batch(iniFile, mode, year_start, year_end, doy_start, doy_end, markerR, sessionR, extR, markerM, sessionM, extM, idN, sessionN, extN, folder);
%
% INPUT:
%   iniFile    = full path to the goGPS .ini file                       [string]
%   mode       = functioning mode (see goGNSS section 'goGPS MODES')    [integer]
%   year_start = year (start)                                           [integer]
%   year_end   = year (end)                                           [integer]
%   doy_start  = day-of-year (start)                                    [integer]
%   doy_end    = day-of-year (end)                                      [integer]
%   markerR    = rover marker name (e.g. 'UBLX')                        [string]
%   sessionR   = rover session id (e.g. '0')                            [string]
%   extR       = rover observation id in file extension (e.g. 'o')      [string]
%   markerM    = master marker name (e.g. 'mila')                       [string]
%   sessionM   = master session id (e.g. '_15s')                        [string]
%   extM       = master observation id in file extension (e.g. 'O')     [string]
%   idN        = navigation id in filename (e.g. 'brdc', 'igs', ...)    [string]
%   sessionN   = navigation session id (e.g. '0')                       [string]
%   extN       = navigation id in file extension (e.g. 'n', 'sp3', ...) [string]
%   folder     = name of the directory (under ../data/batch) where to store results [string]
%
% EXAMPLE:
%   rover RINEX filename:      UBLX0510.14o
%   master RINEX filename:     mila051_15s.14O
%   navigation RINEX filename: brdc0510.14n
%
%   call example:
%   goGPS_batch('./settings/test_Polimi_InputFiles.ini', goGNSS.MODE_PP_KF_CP_DD, 2014, 51, 57, 'UBLX', '0', 'o', 'mila', '_15s', 'O', 'brdc', '0', 'n');
%
% DESCRIPTION:
%   Wrapper function to run goGPS multiple times.

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
%global min_ambfloatRMS
addpath(genpath(pwd));

if (~exist('folder', 'var'))
    folder = '';
end

is_batch = 1; %#ok<*NASGU>
folderOUT = ['../data/out/batch/' folder];
if (exist(folderOUT,'dir') ~= 7)
    mkdir('../data/out/',['batch/' folder]);
end

%-------------------------------------------------------------------------------
% Read INI file
%-------------------------------------------------------------------------------

global goIni;
goIni = goIniReader;
goIni.setFileName(iniFile);

%extract user-defined settings from INI file
path_R_obs = goIni.getData('Receivers','data_path');
% file_name = goIni.getData('Receivers','file_name');
% filename_R_obs = [path_R_obs file_name];

path_M_obs = goIni.getData('Master','data_path');
% file_name = goIni.getData('Master','file_name');
% filename_M_obs = [path_M_obs file_name];

path_N_obs = goIni.getData('Navigational','data_path');
file_name = goIni.getData('Navigational','file_name');
filename_nav = [path_N_obs file_name];

data_path = goIni.getData('RefPath','data_path');
file_name = goIni.getData('RefPath','file_name');
filename_ref = [data_path file_name];

data_path = goIni.getData('PCO_PCV_file','data_path');
file_name = goIni.getData('PCO_PCV_file','file_name');
filename_pco = [data_path file_name];

data_path = goIni.getData('OCEAN_LOADING_file','data_path');
file_name = goIni.getData('OCEAN_LOADING_file','file_name');
filename_blq = [data_path file_name];

iono = goIni.getData('ATM_model','iono');
tropo = goIni.getData('ATM_model','tropo');
if (isempty(iono))
    iono_model = 1;
else
    iono_model = iono;
end
if (isempty(tropo))
    tropo_model = 1;
else
    tropo_model = tropo;
end

%-------------------------------------------------------------------------------
% PROCESSING OPTIONS
%-------------------------------------------------------------------------------
GPS_flag = 1;
GLO_flag = 0;
GAL_flag = 0;
BDS_flag = 0;
QZS_flag = 0;
SBS_flag = 0;
[constellations] = goGNSS.initConstellation(GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag);

flag_SBAS = 0;          % apply SBAS corrections --> no=0, yes=1
flag_IAR = 1;           % try to solve integer ambiguities by LAMBDA method --> no=0, yes=1

min_epoch = 1440;       % minimum number of observed epoch to process

seamless_proc = 1;      % seamless RINEX processing (i.e. do not re-initialize the Kalman filter at each DOY)

%-------------------------------------------------------------------------------
% MASTER STATION POSITION
%-------------------------------------------------------------------------------

flag_ms_pos = 0;     % read master station position from RINEX header
pos_M_crd=[];
pos_R_crd=[];

if (~flag_ms_pos)
    sta_coord_file = '../data/stations/stations.crd';
    
    %parse file containing station coordinates
    [markers, coords_X, coords_Y, coords_Z, flags] = textread(sta_coord_file,'%s%f%f%f%d'); %#ok<DTXTRD>
    
    % master
    %find the correct marker
    marker_idx = find(strcmp(markers, markerM));
    
    %extract the corresponding coordinates
    XM = coords_X(marker_idx);
    YM = coords_Y(marker_idx);
    ZM = coords_Z(marker_idx);
    
    %set master station position
    pos_M_crd = [XM; YM; ZM];
    
    flag_XM = flags(marker_idx);

    % rover
    %find the correct marker
    marker_idx = find(strcmp(markers, markerR));
    
    if ~isempty(marker_idx)
        %extract the corresponding coordinates
        XR = coords_X(marker_idx);
        YR = coords_Y(marker_idx);
        ZR = coords_Z(marker_idx);
        
        %set rover station position
        pos_R_crd = [XR; YR; ZR];
        
        flag_XR = flags(marker_idx);
    end
end

%-------------------------------------------------------------------------------
% KALMAN FILTER
%-------------------------------------------------------------------------------

global sigmaq0 sigmaq_vE sigmaq_vN sigmaq_vU sigmaq_vel
global sigmaq_cod1 sigmaq_cod2 sigmaq_codIF sigmaq_ph sigmaq_phIF sigmaq0_N sigmaq_dtm  sigmaq0_tropo sigmaq_tropo sigmaq_rclock
global min_nsat cutoff snr_threshold cs_threshold_preprocessing cs_threshold weights snr_a snr_0 snr_1 snr_A order o1 o2 o3
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
sigmaq_codIF = 1.2^2;

%variance of phase observations [m^2]
%(maximize to obtain a code-based solution)
sigmaq_ph = 0.003^2;
% sigmaq_ph = 0.001e30;
sigmaq_phIF = 0.009^2;

%variance of ambiguity combinations [cycles]
sigmaq0_N = 1000;

%variance of DEM height [m^2]
%(maximize to disable DEM usage)
% sigmaq_dtm = 0.09;
sigmaq_dtm = 1e30;

%variance of a priori tropospheric delay
sigmaq0_tropo = 1e-2;

%variance of tropospheric delay
sigmaq_tropo = 2.0834e-07;

%variance of receiver clock
sigmaq_rclock = 1e3;

%minimum number of satellites to be used in the Kalman filter
min_nsat = 2;

%cut-off [degrees]
cutoff = 10;

%initialization cut-off [degrees]
% cutoff_init = 15;

%signal-to-noise ratio threshold [dB]
snr_threshold = 0;

%cycle slip threshold (pre-processing) [cycles]
cs_threshold_preprocessing = 0.5;

%cycle slip threshold (processing) [cycles]
cs_threshold = 0.5;

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
% THRESHOLDS
%-------------------------------------------------------------------------------
global SPP_threshold max_code_residual max_phase_residual

SPP_threshold = 4;         %threshold on the code point-positioning least squares
                           %   estimation error [m]
                          
max_code_residual = 30;    %threshold on the maximum residual of code
                           %   observations [m]
                          
max_phase_residual = 0.05; %threshold on the maximum residual of phase
                           %   observations [m]

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
% ITERATIVE SETTING OF INPUT/OUTPUT FILENAME PREFIX AND goGPS LAUNCH
%-------------------------------------------------------------------------------


if (~isempty(markerM))
    markerM_undersc = [markerM '_'];
else
    markerM_undersc = [];
end
fid_extract = fopen([folderOUT '/' markerM_undersc markerR '_' num2str(year_start,'%04d') num2str(doy_start,'%03d') '-' num2str(year_end,'%04d') num2str(doy_end,'%03d') '_extraction.txt'],'w');

fid_extract_TRP = fopen([folderOUT '/' markerM_undersc markerR '_' num2str(year_start,'%04d') num2str(doy_start,'%03d') '-' num2str(year_end,'%04d') num2str(doy_end,'%03d') '_troposphere.txt'],'w');

fid_extract_POS = fopen([folderOUT '/' markerM_undersc markerR '_' num2str(year_start,'%04d') num2str(doy_start,'%03d') '-' num2str(year_end,'%04d') num2str(doy_end,'%03d') '_position.txt'],'w');
fprintf(fid_extract_POS,' yyyy-ddd   date          time           UTM east         UTM north      ellips. height        ZTD\n'); 
fprintf(fid_extract_POS,'+--------+----------+----------------+----------------+----------------+----------------+----------------+\n');

fid_extract_OBS = fopen([folderOUT '/' markerM_undersc markerR '_' num2str(year_start,'%04d') num2str(doy_start,'%03d') '-' num2str(year_end,'%04d') num2str(doy_end,'%03d') '_qualityOBS.txt'],'w');
fprintf(fid_extract_OBS,' yy-ddd  Rover observation file            Rate  #Sat   #Epoch    #Frq   #C1/P1  #C2/P2     #L1     #L2   #DOP1   #DOP2  %%Epoch %%L2/L1    Master observation file            Rate  #Sat   #Epoch    #Frq   #C1/P1  #C2/P2     #L1     #L2   #DOP1   #DOP2  %%Epoch %%L2/L1\n'); 
fprintf(fid_extract_OBS,'+------+------------------------------+--------+-----+--------+-------+--------+-------+-------+-------+-------+-------+-------+------+---------------------------------+--------+-----+--------+-------+--------+-------+-------+-------+-------+-------+-------+------+\n');

n_session=1;
if (seamless_proc)
    kalman_initialized = 0;
end

for year = year_start : 1 : year_end
    year4 = year;
    year2 = two_digit_year(year);
    
    if (year ~= year_start)
        doy_s = 1;
    else
        doy_s = doy_start;
    end
    
    if (year ~= year_end)
        doy_e = 366;
    else
        doy_e = doy_end;
    end
    
    for doy = doy_s : 1 : doy_e
        for session_i=1:n_session
            delete('../data/EMS/*.ems');
            if n_session > 1
                sessionR=num2str(session_i,'%d');
                sessionM=num2str(session_i,'%d');
            end
            if (strcmp(path_R_obs(end),'/') || strcmp(path_R_obs(end),'\'))
                delim = '';
            else
                delim = '/';
            end
            filename_R_obs = [path_R_obs delim markerR num2str(doy,'%03d') sessionR '.' num2str(year2,'%02d') extR];
            filename_M_obs = [path_M_obs delim markerM num2str(doy,'%03d') sessionM '.' num2str(year2,'%02d') extM];
            if (strcmpi(extN,'sp3'))
                [gps_week, ~, gps_dow] = date2gps(datevec(doy2date(doy,year4)));
                filename_nav = [path_N_obs delim idN num2str(gps_week,'%04d') num2str(gps_dow,'%01d') '.' extN];
            else
                filename_nav = [path_N_obs delim idN num2str(doy,'%03d') sessionN '.' num2str(year2,'%02d') extN];
            end
            prefixOUT = [markerM_undersc markerR '_' num2str(year2,'%02d') num2str(doy,'%03d') num2str(session_i,'%d')];
            if exist(filename_R_obs,'file')>0 && exist(filename_nav,'file')>0 && (exist(filename_M_obs,'file')>0 || isempty(markerM))
                
                %-------------------------------------------------------------------------------
                % INPUT/OUTPUT FILENAME PREFIX
                %-------------------------------------------------------------------------------
                
                %dummy values for filerootIN (not used for batch processing)
                folderIN  = '../data/data_goGPS';
                prefixIN  = 'yamatogawa';
                filerootIN  = [folderIN '/' prefixIN];
                filerootOUT = [folderOUT '/' prefixOUT];
                
                i = 1;
                
                j = length(filerootOUT);
                while (~isempty(dir([filerootOUT '_rover*.bin'])) || ...
                        ~isempty(dir([filerootOUT '_master*.bin'])) || ...
                        ~isempty(dir([filerootOUT '_obs*.bin'])) || ...
                        ~isempty(dir([filerootOUT '_eph*.bin'])) || ...
                        ~isempty(dir([filerootOUT '_sat*.bin'])) || ...
                        ~isempty(dir([filerootOUT '_kal*.bin'])) || ...
                        ~isempty(dir([filerootOUT '_dt*.bin'])) || ...
                        ~isempty(dir([filerootOUT '_conf*.bin'])) || ...
                        ~isempty(dir([filerootOUT '_dop*.bin'])) || ...
                        ~isempty(dir([filerootOUT '_ECEF*.txt'])) || ...
                        ~isempty(dir([filerootOUT '_geod*.txt'])) || ...
                        ~isempty(dir([filerootOUT '_plan*.txt'])) || ...
                        ~isempty(dir([filerootOUT '_NMEA*.txt'])) || ...
                        ~isempty(dir([filerootOUT '_ublox_NMEA*.txt'])) || ...
                        ~isempty(dir([filerootOUT '.kml'])) )
                    
                    filerootOUT(j+1:j+4) = ['_' num2str(i,'%03d')];
                    i = i + 1;
                end
                
                clear X_KAL
                clear Xhat_t_t
                goGPS
                
                idx = size(date_R,1);
                if exist('X_KAL','var') && exist('Xhat_t_t_OUT','var')
                    fprintf(fid_extract,'%04d-%03d  %02d/%02d/%02d    %02d:%02d:%06.3f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f\n', year4, doy, date_R(idx,1), date_R(idx,2), date_R(idx,3), date_R(idx,4), date_R(idx,5), date_R(idx,6), X_KAL(idx), Y_KAL(idx), Z_KAL(idx), EAST_UTM(idx), NORTH_UTM(idx), h_KAL(idx));
                    tropo_vec = nan(1,86400/interval); tropo_vec(1,1:length(Xhat_t_t_OUT(end-1,:))) = Xhat_t_t_OUT(end-1,:); %#ok<NODEF>
                    fprintf(fid_extract_TRP,'%.6f ', tropo_vec);
                    fprintf(fid_extract_TRP,'\n');
                    for e = 1 : idx
                        fprintf(fid_extract_POS,' %04d-%03d  %02d/%02d/%02d    %02d:%02d:%06.3f %16.6f %16.6f %16.6f %15.6f\n', year4, doy, date_R(e,1), date_R(e,2), date_R(e,3), date_R(e,4), date_R(e,5), date_R(e,6), EAST_UTM(e), NORTH_UTM(e), h_KAL(e), Xhat_t_t_OUT(end-1,e));
                    end
                    delete([filerootOUT '_*.bin']);
                else
                    fprintf(fid_extract,'%04d-%03d\n', year4, doy);
                end
                
                % append report information in the database
                fid_rep_i=fopen([filerootOUT,'_report.txt'],'rt');
                if fid_rep_i~=-1
                    line=fgetl(fid_rep_i);
                    while (isempty(strfind(line,'Observations (RAW)              Start time')))
                        line = fgetl(fid_rep_i);
                    end
                    line=fgetl(fid_rep_i);
                    line1=fgetl(fid_rep_i);
                    line1(33:80)='';
                    line2=fgetl(fid_rep_i);
                    if ~isempty(line2)
                        line2(33:80)='';
                    end
                    
                    fprintf(fid_extract_OBS,' %s-%s  %s     %s\n',num2str(year2,'%02d'),num2str(doy,'%03d'), line1, line2);
                    fclose(fid_rep_i);
                end
            end
        end
    end
end

fclose(fid_extract);
fclose(fid_extract_TRP);
fclose(fid_extract_POS);
fclose(fid_extract_OBS);

figure
tropo = load([folderOUT '/' markerM_undersc markerR '_' num2str(year_start,'%04d') num2str(doy_start,'%03d') '-' num2str(year_end,'%04d') num2str(doy_end,'%03d') '_troposphere.txt']);
plot(reshape(tropo',1,size(tropo,1)*size(tropo,2)));
