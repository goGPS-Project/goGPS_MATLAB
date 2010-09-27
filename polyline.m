%-----------------------------------------------------------
% PARAMETERS
%-----------------------------------------------------------

% input data file names
rootname = 'campetto_ublox_COMO_cv';
dat_filename = ['../data/' rootname '_position.txt'];
cov_filename = ['../data/' rootname '_cov.txt'];

tab_filename = ['../data/' rootname '_table.txt'];
nod_filename = ['../data/' rootname '_node.txt'];

dN1 = 0; % disregarded points at the beginning
dN2 = 5; % disregarded points at the end

% threshold on the angle between arcs
angle_threshold = 170 * pi/180;

% threshold on the distance between nodes (AGNES method)
dist_threshold_AGNES = 2;

% bounding box dimension
delta_iter0 = 1;
delta_iter1 = 0.5; % (typically to be reduced)

% threshold on the distance between old/new nodes
dist_threshold_update_iter0 = 2;
dist_threshold_update_iter1 = 2; % (typically the same)

% weighted observations
flag_iter0 = 1;   % (0/1) = (independent/weighted)
flag_iter1 = 1;

%-----------------------------------------------------------
% STEP 1 - ITERATION 0
%-----------------------------------------------------------

[nodes] = polyline_nodesDetection (dat_filename, dN1, dN2, angle_threshold, dist_threshold_AGNES);

%-----------------------------------------------------------
% STEP 2 - ITERATION 0
%-----------------------------------------------------------

[table, nodes] = polyline_arcsClustering (dat_filename, cov_filename, tab_filename, nodes, dN1, dN2, delta_iter0);

%-----------------------------------------------------------
% STEP 3 - ITERATION 0
%-----------------------------------------------------------

[nodes] = polyline_leastSquaresFit (tab_filename, nod_filename, nodes, flag_iter0, dist_threshold_update_iter0);

%-----------------------------------------------------------
% STEP 2 - ITERATION 1
%-----------------------------------------------------------

[table, nodes] = polyline_arcsClustering (dat_filename, cov_filename, tab_filename, nodes, dN1, dN2, delta_iter1);

%-----------------------------------------------------------
% STEP 3 - ITERATION 1
%-----------------------------------------------------------

[nodes] = polyline_leastSquaresFit (tab_filename, nod_filename, nodes, flag_iter1, dist_threshold_update_iter1);

%-----------------------------------------------------------