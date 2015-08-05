function [nodes, utm_zone] = polyline_nodesDetection (filename, dN1, dN2, angle_threshold, dist_threshold, min_nodes)

% SYNTAX:
%   [nodes, utm_zone] = polyline_nodesDetection (filename, dN1, dN2, angle_threshold, dist_threshold, min_nodes);
%
% INPUT:
%   filename = input data file name
%   dN1 = disregarded points at the beginning
%   dN2 = disregarded points at the end
%   angle_threshold = threshold on the angle between arcs [radiants]
%   dist_threshold = threshold on the distance between nodes [meters]
%   min_nodes = minimum number of nodes
%
% OUTPUT:
%   nodes = coordinates of the selected nodes
%   utm_zone = UTM zone of the node coordinates
%
% DESCRIPTION:
%   Determine the nodes using an agglomerative method.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Lisa Pertusini, Alemu Befkadu
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
% loading the data
%-----------------------------------------------------------

fid = fopen(filename,'rt');      % open file
fgets(fid);                      % jump the header
data = fscanf(fid,'%d/%d/%d %d:%d:%f %d %f %f %f %f %f %f %f %f %f %f %4c %f %f %f %f %f %d %f',[28 inf])';
fclose(fid);
clear filename

%-----------------------------------------------------------
% assigning the data
%-----------------------------------------------------------

x0 = data(:,16);       % x coordinates
y0 = data(:,15);       % y coordinates
utm_zone = char([data(:,18) data(:,19) data(:,20) data(:,21)]); % UTM zones
clear data

%-----------------------------------------------------------
% translate the coordinates in order to have better numbers
%-----------------------------------------------------------

x0_start = x0(1,1);   % initial point (no translation)
y0_start = y0(1,1);   % final point (no translation)

x0 = x0 - x0_start;    % new x coordinates
y0 = y0 - y0_start;    % new y coordinates

%-----------------------------------------------------------
% data selection
%-----------------------------------------------------------

N = length(x0); % number of data

N1 = 1+dN1;     % disregard initial data
N2 = N-dN2;     % disregard final data
N = N2-N1+1;    % new number of data

%-----------------------------------------------------------
% dissimilarity table
%-----------------------------------------------------------

pos0 = (N1   : N2-2)';
pos1 = (N1+1 : N2-1)';
pos2 = (N1+2 : N2  )';

dx1 = x0(pos1) - x0(pos0);
dx2 = x0(pos1) - x0(pos2);
dy1 = y0(pos1) - y0(pos0); 
dy2 = y0(pos1) - y0(pos2);
clear pos0 pos1 pos2

% initialization
table = zeros(N,3);
% index of the segment
table(1:N,1) = (N1 : N2)';
% calculating the angle
table(2:N-1,2) = acos((dx1.*dx2 + dy1.*dy2) ./ (sqrt(dx1.^2 + dy1.^2) .* sqrt(dx2.^2 + dy2.^2))); 

clear dx1 dx2
clear dy1 dy2

%-----------------------------------------------------------
% agglomerative nesting (AGNES)
%-----------------------------------------------------------

% finding the maximum angle
[null_max pos] = max(table(:,2)); %#ok<ASGLU>

while (table(pos,2) > angle_threshold) & (size(table,1) > min_nodes)

    table(pos,:) = []; % merging or combining

    %---------------------------------
    % plotting the data

    % plot(x0(N1:N2),y0(N1:N2),'-b')
    % hold on
    % plot(x0(table(:,1)),y0(table(:,1)), '.-k');
    % hold off
    % axis equal
    % pause
    %---------------------------------

    % updating the dissimilarity table
    if (pos ~= 2) & (pos ~= size(table,1))
        j = pos-2;
        c1 = (x0( table(j+1,1) )-x0( table(j,1) )) * (x0( table(j+1,1) )-x0( table(j+2,1) ));
        c2 = (y0( table(j+1,1) )-y0( table(j,1) )) * (y0( table(j+1,1) )-y0( table(j+2,1) ));
        d1 = sqrt((x0( table(j+1,1) )-x0( table(j,1) ))^2 + (y0( table(j+1,1) )-y0( table(j,1) ))^2);
        d2 = sqrt((x0( table(j+1,1) )-x0( table(j+2,1) ))^2 + (y0( table(j+1,1) )-y0( table(j+2,1) ))^2);
        table(pos-1,2) = acos((c1 + c2) ./ (d1 * d2));

        j = pos-1;
        c1 = (x0( table(j+1,1) )-x0( table(j,1) )) * (x0( table(j+1,1) )-x0( table(j+2,1) ));
        c2 = (y0( table(j+1,1) )-y0( table(j,1) )) * (y0( table(j+1,1) )-y0( table(j+2,1) ));
        d1 = sqrt((x0( table(j+1,1) )-x0( table(j,1) ))^2 + (y0( table(j+1,1) )-y0( table(j,1) ))^2);
        d2 = sqrt((x0( table(j+1,1) )-x0( table(j+2,1) ))^2 + (y0( table(j+1,1) )-y0( table(j+2,1) ))^2);
        table(pos,2) = acos((c1 + c2) ./ (d1 * d2));
    elseif (pos == 2)
        j = pos-1;
        c1 = (x0( table(j+1,1) )-x0( table(j,1) )) * (x0( table(j+1,1) )-x0( table(j+2,1) ));
        c2 = (y0( table(j+1,1) )-y0( table(j,1) )) * (y0( table(j+1,1) )-y0( table(j+2,1) ));
        d1 = sqrt((x0( table(j+1,1) )-x0( table(j,1) ))^2 + (y0( table(j+1,1) )-y0( table(j,1) ))^2);
        d2 = sqrt((x0( table(j+1,1) )-x0( table(j+2,1) ))^2 + (y0( table(j+1,1) )-y0( table(j+2,1) ))^2);
        table(pos,2) = acos((c1 + c2) ./ (d1 * d2));
    else
        j = pos-2;
        c1 = (x0( table(j+1,1) )-x0( table(j,1) )) * (x0( table(j+1,1) )-x0( table(j+2,1) ));
        c2 = (y0( table(j+1,1) )-y0( table(j,1) )) * (y0( table(j+1,1) )-y0( table(j+2,1) ));
        d1 = sqrt((x0( table(j+1,1) )-x0( table(j,1) ))^2 + (y0( table(j+1,1) )-y0( table(j,1) ))^2);
        d2 = sqrt((x0( table(j+1,1) )-x0( table(j+2,1) ))^2 + (y0( table(j+1,1) )-y0( table(j+2,1) ))^2);
        table(pos-1,2) = acos((c1 + c2) ./ (d1 * d2));
    end

    % finding the maximum angle
    [null_max pos] = max(table(:,2)); %#ok<ASGLU>

    %---------------------------------
    % plotting the node

    % hold on
    % plot(x0(table(pos,1)),y0(table(pos,1)),'.r');
    % hold off
    % axis equal
    % pause
    %---------------------------------
end

clear c1 c2
clear d1 d2
clear pos

%-----------------------------------------------------------
% displaying result
%-----------------------------------------------------------

nodes(:,1) = x0(table(:,1));     % node - east coordinates
nodes(:,2) = y0(table(:,1));     % node - north coordinates

figure
plot(x0(N1:N2),y0(N1:N2),'.k')
hold on 
plot(nodes(:,1),nodes(:,2),'*-g')
%axis equal
title('Selected nodes by AGNES');
xlabel('East [m]');
ylabel('North [m]');

clear nodes

%-----------------------------------------------------------
% merging if too close
%-----------------------------------------------------------

i = 1;
while i < (size(table,1))
    j = i+1;
    while (j <= size(table,1) & (sqrt((x0(table(i,1)) - x0(table(j,1)))^2 + (y0(table(i,1)) - y0(table(j,1)))^2) < dist_threshold))
        table(j,:) = [];
    end
    i = i+1;
end

clear i j

%-----------------------------------------------------------
% displaying result
%-----------------------------------------------------------

figure
plot(x0(N1:N2),y0(N1:N2),'.k')
hold on 
plot(x0(table(:,1)),y0(table(:,1)),'*-g')
%axis equal
title('Selected nodes by AGNES after distance test');
xlabel('East [m]');
ylabel('North [m]');

%-----------------------------------------------------------
% output data
%-----------------------------------------------------------

nodes(:,1) = x0(table(:,1)) + x0_start;     % node - east coordinates
nodes(:,2) = y0(table(:,1)) + y0_start;     % node - north coordinates

% num_nodes = length(nodes);     % number of nodes

%-----------------------------------------------------------