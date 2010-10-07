function polyline(filerootIN, angle_threshold)

% SYNTAX:
%   polyline(filerootIN, angle_threshold);
%
% INPUT:
%   filerootIN = input file root (rover data, binary stream)
%   angle_threshold = threshold on the angle between arcs [degrees]
%
% DESCRIPTION:
%   Polyline simplification algorithm.

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

%-----------------------------------------------------------
% PARAMETERS
%-----------------------------------------------------------

% input data file names
dat_filename = [filerootIN '_position.txt'];
cov_filename = [filerootIN '_cov.txt'];

% intermediate processing file names
tab_filename = [filerootIN '_table.txt'];
nod_filename = [filerootIN '_node.txt'];

dN1 = 0; % disregarded points at the beginning
dN2 = 5; % disregarded points at the end

% threshold on the angle between arcs
angle_threshold = angle_threshold * pi/180;

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