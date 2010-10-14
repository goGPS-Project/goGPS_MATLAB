function polyline(filerootIN, angle_threshold, dist_threshold_AGNES, dN1, dN2, ...
    delta_iter0, delta_iter1, dist_threshold_update_iter0, dist_threshold_update_iter1, ...
    flag_iter0, flag_iter1, min_nodes)

% SYNTAX:
%   polyline(filerootIN, angle_threshold, dN1, dN2, ...
%   delta_iter0, delta_iter1, dist_threshold_update_iter0, dist_threshold_update_iter1, ...
%   flag_iter0, flag_iter1, min_nodes);
%
% INPUT:
%   filerootIN = input file root (rover data, binary stream)
%   angle_threshold = threshold on the angle between arcs [degrees]
%   dist_threshold_AGNES = threshold on the distance between nodes (AGNES method)
%   dN1 = number of neglected points at the beginning of the path
%   dN2 = number of neglected points at the end of the path
%   delta_iter0 = bounding box dimension (iteration 0)
%   delta_iter1 = bounding box dimension (iteration 1) (typically to be reduced)
%   dist_threshold_update_iter0 = threshold on the distance between old/new nodes (iteration 0)
%   dist_threshold_update_iter1 = threshold on the distance between old/new nodes (iteration 1)
%   flag_iter0 = weighted observations flag (0/1) = (independent/weighted) (iteration 0)
%   flag_iter1 = weighted observations flag (0/1) = (independent/weighted) (iteration 1)
%   min_nodes = minimum number of nodes (threshold on "compression level")
%
% DESCRIPTION:
%   Polyline simplification algorithm.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.2 alpha
%
% Copyright (C) 2009-2010 Mirko Reguzzoni, Eugenio Realini
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

% threshold on the angle between arcs
angle_threshold = angle_threshold * pi/180;

%-----------------------------------------------------------
% STEP 1 - ITERATION 0
%-----------------------------------------------------------

[nodes, utm_zone] = polyline_nodesDetection (dat_filename, dN1, dN2, angle_threshold, dist_threshold_AGNES, min_nodes);

%-----------------------------------------------------------
% STEP 2 - ITERATION 0
%-----------------------------------------------------------

[table, nodes] = polyline_arcsClustering (dat_filename, cov_filename, tab_filename, nodes, dN1, dN2, delta_iter0); %#ok<ASGLU>

%-----------------------------------------------------------
% STEP 3 - ITERATION 0
%-----------------------------------------------------------

[nodes] = polyline_leastSquaresFit (tab_filename, nod_filename, nodes, flag_iter0, dist_threshold_update_iter0);

%-----------------------------------------------------------
% STEP 2 - ITERATION 1
%-----------------------------------------------------------

[table, nodes] = polyline_arcsClustering (dat_filename, cov_filename, tab_filename, nodes, dN1, dN2, delta_iter1); %#ok<ASGLU>

%-----------------------------------------------------------
% STEP 3 - ITERATION 1
%-----------------------------------------------------------

[nodes] = polyline_leastSquaresFit (tab_filename, nod_filename, nodes, flag_iter1, dist_threshold_update_iter1);

%-----------------------------------------------------------

%-------------------------------------------------------------------------
% print the result in a KML file

line_colorR = 'ff00ffff';
line_width = 2;

pos = find(filerootIN == '/');
kml_name = filerootIN(pos(end)+1:end);

fid_kml = fopen([filerootIN '_polyline.kml'], 'wt');
    fprintf(fid_kml, '<?xml version="1.0" encoding="UTF-8"?>\n');
    fprintf(fid_kml, '<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">\n');
    fprintf(fid_kml, '<Document>\n');
    fprintf(fid_kml, '\t<name>%s</name>\n', [kml_name '_polyline.kml']);
    fprintf(fid_kml, '\t<snippet>created by goGPS</snippet>\n');
    fprintf(fid_kml, '\t<Style id="goLine1">\n');
    fprintf(fid_kml, '\t\t<LineStyle>\n');
    fprintf(fid_kml, '\t\t\t<color>%s</color>\n',line_colorR);
    fprintf(fid_kml, '\t\t\t<width>%d</width>\n' ,line_width);
    fprintf(fid_kml, '\t\t</LineStyle>\n');
    fprintf(fid_kml, '\t</Style>\n');
    fprintf(fid_kml, '\t<Placemark>\n');
    fprintf(fid_kml, '\t\t<name>Rover track (simplified polyline)</name>\n');
    fprintf(fid_kml, '\t\t<styleUrl>#goLine1</styleUrl>\n');
    fprintf(fid_kml, '\t\t<LineString>\n');
    fprintf(fid_kml, '\t\t\t<coordinates>\n\t\t\t\t');
    for i = 1 : length(nodes)
        [nodes_lat, nodes_lon] = utm2deg(nodes(i,1), nodes(i,2), utm_zone(i,:));
        fprintf(fid_kml, '%.8f,%.8f,0 ',nodes_lon,nodes_lat);
    end
    fprintf(fid_kml, '\n\t\t\t</coordinates>\n');
    fprintf(fid_kml, '\t\t</LineString>\n');
    fprintf(fid_kml, '\t</Placemark>\n');
    fprintf(fid_kml, '</Document>\n</kml>');
    fclose(fid_kml);

end

%-------------------------------------------------------------------------