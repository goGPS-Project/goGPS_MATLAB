% SYNTAX:
%   global_init;
%
% DESCRIPTION:
%   Global variables initialization.

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
% CONSTANTS
%-------------------------------------------------------------------------------

% global v_light
% global f1 f2
% global lambda1 lambda2
global a_GPS a_GLO f_GPS e_GPS
global GM_GPS GM_GLO GM_GAL GM_BDS GM_QZS
global Omegae_dot_GPS Omegae_dot_GLO Omegae_dot_GAL Omegae_dot_BDS Omegae_dot_QZS
global J2_GLO
global circle_rad

%velocity of light in the void
% v_light = 299792458; % [m/s]

%GPS carriers frequencies
% f1 = 1575420000; % [1/s]
% f2 = 1227600000; % [1/s]

%GPS carriers wavelengths
% lambda1 = v_light / f1; % [m]
% lambda2 = v_light / f2; % [m]

%CRS parameters, according to each GNSS system CRS definition
% (ICD document in brackets):
%
% *_GPS --> WGS-84   (IS-GPS200E)
% *_GLO --> PZ-90    (GLONASS-ICD 5.1)
% *_GAL --> GTRF     (Galileo-ICD 1.1)
% *_BDS --> CGCS2000 (BeiDou-ICD 1.0)
% *_QZS --> WGS-84   (IS-QZSS 1.5D)

%ellipsoid semi-major axis [m]
a_GPS = 6378137;
a_GLO = 6378136;
% a_GAL = 6378137;
% a_BDS = 6378137;
% a_QZS = 6378137;

%ellipsoid flattening
f_GPS = 1/298.257223563;

%eccentricity
e_GPS = sqrt(1-(1-f_GPS)^2);

%gravitational constant (mass of Earth) [m^3/s^2]
GM_GPS = 3.986005e14;    
GM_GLO = 3.9860044e14;
GM_GAL = 3.986004418e14; 
GM_BDS = 3.986004418e14;
GM_QZS = 3.986005e14;

%angular velocity of the Earth rotation [rad/s]
Omegae_dot_GPS = 7.2921151467e-5;
Omegae_dot_GLO = 7.292115e-5;
Omegae_dot_GAL = 7.2921151467e-5;
Omegae_dot_BDS = 7.292115e-5;
Omegae_dot_QZS = 7.2921151467e-5;

%second zonal harmonic of the geopotential
J2_GLO = 1.0826257e-3;

%pi value used for orbit computation
pi_orbit   = 3.1415926535898;
circle_rad = 2*pi_orbit;

%-------------------------------------------------------------------------------
% SATELLITE CONFIGURATION
%-------------------------------------------------------------------------------

global azR elR distR
global azM elM distM
global elea

%azimuth, elevation and distance of satellites with respect to the ROVER
azR = [];
elR = [];
distR = [];

%azimuth, elevation and distance of satellites with respect to the MASTER
azM = [];
elM = [];
distM = [];

elea  = 10; % default value for the exponential elevation weight function

%-------------------------------------------------------------------------------
% DILUTION OF PRECISION
%-------------------------------------------------------------------------------

global PDOP HDOP VDOP KPDOP KHDOP KVDOP

PDOP = [];
HDOP = [];
VDOP = [];
KPDOP = [];
KHDOP = [];
KVDOP = [];

%-------------------------------------------------------------------------------
% LAMBDA STATISTIC
%-------------------------------------------------------------------------------

global success

success = [];

%-------------------------------------------------------------------------------
% KALMAN FILTER
%-------------------------------------------------------------------------------

%transition matrix
global T

%identity matrix
global I

%state estimation at time t
global Xhat_t_t

%state estimation at time t+1 (using dynamics only)
global X_t1_t

%estimation error covariance matrix at time t
global Cee

%number of visible satellites at time t
global nsat  %%% Developer's note: OBSOLETE %%%

%satellite configuration (1: visible, 0: not visible) at time t
global conf_sat

%cycle-slip configuration (1: cs, 0: no cs) at time t
global conf_cs

%index of the current pivot satellite
global pivot

%index of the previous pivot satellite
global pivot_old

%number of unknown phase ambiguities
global nN

%method used to estimate phase ambiguities
%          - amb_estim_method=0: observed code - phase difference
%          - amb_estim_method=1: Kalman-predicted code - phase difference
%          - amb_estim_method=2: Least squares adjustment
global amb_estim_method

%interval between epochs
global interval

%initialization
T = [];
I = [];
Xhat_t_t = [];
X_t1_t = [];
Cee = [];
nsat = [];
conf_sat = [];
conf_cs = [];
pivot = [];
pivot_old = [];
nN = [];
amb_estim_method = 0;
interval = 1; %default 1 Hz (to avoid problems with real-time modes)

%-------------------------------------------------------------------------------
% KALMAN FILTER (CONSTRAINED VERSION)
%-------------------------------------------------------------------------------

%constrained trajectory
global s0 %#ok<NUSED>
global ax ay az %#ok<NUSED>
global Yhat_t_t Y_t1_t %#ok<NUSED>

%-------------------------------------------------------------------------------
% REAL-TIME MANAGEMENT
%-------------------------------------------------------------------------------
global master rover %#ok<NUSED>

%-------------------------------------------------------------------------------
% RECEIVER
%-------------------------------------------------------------------------------

global rec_clock_error
global flag_doppler_cs
global doppler_pred_range1_R
global doppler_pred_range2_R
global doppler_pred_range1_M
global doppler_pred_range2_M

%receiver clock error
rec_clock_error = 0;

%flag to enable Doppler-based cycle slip detection
flag_doppler_cs = 0;

%Doppler-predicted range (ROVER)
doppler_pred_range1_R = zeros(nSatTot,1);
doppler_pred_range2_R = zeros(nSatTot,1);

%Doppler-predicted range (MASTER)
doppler_pred_range1_M = zeros(nSatTot,1);
doppler_pred_range2_M = zeros(nSatTot,1);

%-------------------------------------------------------------------------------
% MASTER STATION
%-------------------------------------------------------------------------------

global nmea_update_rate

%NMEA update rate (waiting time for sending a new $GGA string to NTRIP caster)
nmea_update_rate = 10; %[sec]

%-------------------------------------------------------------------------------
% UBX POLL RATES
%-------------------------------------------------------------------------------

global hui_poll_rate

%AID-HUI message poll rate
hui_poll_rate = 3600; %[sec]

%-------------------------------------------------------------------------------
% SERVER CONNECTION
%-------------------------------------------------------------------------------

global server_delay

%server waiting time (to check if packet transmission is finished)
server_delay = 0.05;

%-------------------------------------------------------------------------------
% MATLAB DISPLAY
%-------------------------------------------------------------------------------

global p_max pid
global satid labid pivid
global window
global EAST_O NORTH_O
global x_circle id_ellipse

%maximum number of drawn trajectory points
p_max = 200;

%trajectory point id
pid = zeros(p_max,1);

%sky-plot id
satid = zeros(nSatTot,1);
labid = zeros(nSatTot,1);
pivid = 0;

%master station point id
msid = [];

%dimension of the time windows (ambiguity plot)
window = 20;

%reference origin
EAST_O = 0;
NORTH_O = 0;

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
% GOOGLE EARTH
%-------------------------------------------------------------------------------

global link_filename kml_filename

%files used by Google Earth
link_filename = '../data/google_earth/link.kml';
kml_filename = '../data/google_earth/goGPS.kml';

%-------------------------------------------------------------------------------
% THRESHOLDS
%-------------------------------------------------------------------------------
global clock_delay_thresh
global cond_num_threshold

clock_delay_thresh = 100;
cond_num_threshold = 1e6; %threshold on the condition number on the
                          % eigenvalues of the N matrix (least squares)

%-------------------------------------------------------------------------------
% PHASE-SMOOTHED CODE
%-------------------------------------------------------------------------------

global sm_weight

%weight for code smoothing algorithm
sm_weight = 1;

%-------------------------------------------------------------------------------
% GEOID GRID
%-------------------------------------------------------------------------------

global geoid

try
    load ../data/geoid/geoid_EGM2008_05.mat
    %geoid grid and parameters
    geoid.grid = N_05x05;
    geoid.cellsize = 0.5;
    geoid.Xll = -180;
    geoid.Yll = -90;
    geoid.ncols = 720;
    geoid.nrows = 360;

    clear N_05x05
catch
    %geoid unavailable
    geoid.grid = 0;
    geoid.cellsize = 0;
    geoid.Xll = 0;
    geoid.Yll = 0;
    geoid.ncols = 0;
    geoid.nrows = 0;
end

%-------------------------------------------------------------------------------
% LAMBDA
%-------------------------------------------------------------------------------

global ratiotest mutest succ_rate fixed_solution

ratiotest = [];
mutest = [];
succ_rate = [];
fixed_solution = [];

%-------------------------------------------------------------------------------
% MULTI-CONSTELLATION
%-------------------------------------------------------------------------------

global n_sys

n_sys = 1;

%-------------------------------------------------------------------------------
% Fisher test table
%-------------------------------------------------------------------------------
global FTABLE;
FTABLE=finv(0.9995,1,1:200)';
