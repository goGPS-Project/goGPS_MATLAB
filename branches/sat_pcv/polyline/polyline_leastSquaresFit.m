function [nodes] = polyline_leastSquaresFit (tab_filename, nod_filename, nodes, flag, dist_threshold)

% SYNTAX:
%   [nodes] = polyline_leastSquaresFit (filename, nodes);
%
% INPUT:
%   tab_filename = input data file name (table)
%   nod_filename = input data file name (node)
%   nodes = input coordinates of the nodes
%   flag = independent/weighted observation (0/1)
%   dist_threshold = threshold on the distance between nodes [meters]
%
% OUTPUT:
%   nodes = output coordinates of the nodes
%
% DESCRIPTION:
%   Least-squares adjustment for each arc.

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

% filename = 'Lura_ublox_cv_low_position.txt';
% filename = 'realdata.txt';

%--------------------------------------------------------------------------
% loading the data
%--------------------------------------------------------------------------

data = load(tab_filename);     % east | north | variance | label
clear filename

%--------------------------------------------------------------------------
% assigning the data
%--------------------------------------------------------------------------

x0 = data(:,1);    % x coordinates
y0 = data(:,2);    % y coordinates
vx = data(:,3);    % x data variances
vy = data(:,4);    % y data variances
l = data(:,5);     % data labels

if (flag == 0)
    vx = ones(size(vx));  % no weights in the least-squares adjustment
    vy = ones(size(vy));
end

x0_start = x0(1);      % initial point
y0_start = y0(1);      % final point

x0 = x0 - x0_start;    % new x coordinates
y0 = y0 - y0_start;    % new y coordinates

nodes(:,1) = nodes(:,1) - x0_start;    % new x nodes
nodes(:,2) = nodes(:,2) - y0_start;    % new y nodes

%fprintf('%.3f\n',x0_start)
%fprintf('%.3f\n',y0_start)

%--------------------------------------------------------------------------
% least-squares adjustment
%--------------------------------------------------------------------------

NN = max(l);        % define the number of segments

xm = zeros(NN,1);   % initialization
ym = zeros(NN,1);
par = zeros(NN,2);

for k = 1 : NN;     % k is the index of the segment

    pos = find(l == k);     % indexes of points belonging to the segment k

    xm(k) = mean(x0(pos));  % mean of the observations to find the mid point of the segment
    ym(k) = mean(y0(pos));

    F = (y0(pos) - ym(k)) ./ (x0(pos) - xm(k));
    m_tilde = mean(F);      % approximate value of the angular coefficient

    q = (vy(pos) + m_tilde^2 .* vx(pos)).^(-1);

    % least square solution
    par(k,:) = LSsolver(x0(pos), y0(pos), xm(k), q);

end

clear pos F q
clear m_tilde

%--------------------------------------------------------------------------
% node determination (from least-squares adjustment)
%--------------------------------------------------------------------------

xA = nodes(1:end-1,1);
xB = nodes(2:end,1);

yA = par(:,2) .* (xA-xm) + par(:,1);
yB = par(:,2) .* (xB-xm) + par(:,1);

%--------------------------------------------------------------------------
% node determination (arc intersection or node average)
%--------------------------------------------------------------------------

xNodes = zeros (NN+1,1);    % initialization
yNodes = zeros (NN+1,1);

% for i = 1 : NN-1;
%     xNodes(i+1) = (par(i+1,2) * xm(i+1) - par(i,2) * xm(i) + par(i,1) - par(i+1,1)) / (par(i+1,2) - par(i,2));  % node of each segment in x direction
%     yNodes(i+1) =  par(i,2) * (xNodes(i+1)-xm(i)) + par(i,1);                                                   % node of each segment in y direction
%      if (xNodes(i+1) < min(nodes(i,1),nodes(i+2,1))) | (xNodes(i+1) > max(nodes(i,1),nodes(i+2,1))) | ...
%         (yNodes(i+1) < min(nodes(i,2),nodes(i+2,2))) | (yNodes(i+1) > max(nodes(i,2),nodes(i+2,2)))
%            xNodes(i+1) = xA(i+1);              % node of each segment in x direction
%            yNodes(i+1) = (yA(i+1) + yB(i))/2;  % node of each segment in y direction
%      end
% end

for i = 1 : NN-1;
    xNodes_avg = xA(i+1);              % node of each segment in x direction
    yNodes_avg = (yA(i+1) + yB(i))/2;  % node of each segment in y direction
    xNodes_int = (par(i+1,2) * xm(i+1) - par(i,2) * xm(i) + par(i,1) - par(i+1,1)) / (par(i+1,2) - par(i,2));  % node of each segment in x direction
    yNodes_int = par(i,2) * (xNodes_int-xm(i)) + par(i,1);                                                   % node of each segment in y direction
    % if sqrt((xNodes_avg - xNodes_int)^2 + (yNodes_avg - yNodes_int)^2) > dist_threshold
    if sqrt((nodes(i+1,1) - xNodes_int)^2 + (nodes(i+1,2) - yNodes_int)^2) > dist_threshold
        xNodes(i+1) = xNodes_avg;
        yNodes(i+1) = yNodes_avg;
    else
        xNodes(i+1) = xNodes_int;
        yNodes(i+1) = yNodes_int;
    end
end

xNodes(1) = xA(1);       % initial node
yNodes(1) = yA(1);       % initial node
xNodes(NN+1) = xB(end);  % final node
yNodes(NN+1) = yB(end);  % final node

%--------------------------------------------------------------------------
% plotting the data
%--------------------------------------------------------------------------

figure;
plot(x0,y0,'.k');
hold on
plot(nodes(:,1), nodes(:,2), '*g')

for i = 1 : NN;
    plot([xA(i) xB(i)], [yA(i) yB(i)], '.-b');
end

for i = 1 : NN+1;
    plot(xNodes(i),yNodes(i),'*r');
end

%axis equal
title('Least-squares fitting');
xlabel('East [m]');
ylabel('North [m]');

%------------------------------

figure;
plot(x0,y0,'.k');
hold on
plot(xNodes,yNodes,'*-r');
hold off

%axis equal
title('Final solution');
xlabel('East [m]');
ylabel('North [m]');

% for i = 1 : 60
%     pos = find(l == i);
%     plot(realdata(pos,1),realdata(pos,2),'.b');
%     hold on
%     plot([xA(i) xB(i)], [yA(i); yB(i)], '.-r')
%     plot(xNodes(i),yNodes(i) ,'*c');
%     plot(xNodes(i+1),yNodes(i+1) ,'*c');
%
%     pos = find(l == i+1);
%     plot(table(pos,1),table(pos,2),'.c');
%     hold on
%     plot([xA(i+1) xB(i+1)], [yA(i+1); yB(i+1)], '.-r')
%     plot(xNodes(i+1),yNodes(i+1) ,'*c');
%     plot(xNodes(i+2),yNodes(i+2) ,'*c');
%
%     hold off
%     pause
% end

%-----------------------------------------------------------
% output data

clear nodes
nodes(:,1) = xNodes + x0_start;     % node - east coordinates
nodes(:,2) = yNodes + y0_start;     % node - north coordinates

%-------------------------------------------------------------------------
% print the result in a txt file

fid = fopen(nod_filename,'wt');
fprintf(fid,'%.3f %.3f\n',[nodes(:,1), nodes(:,2)]');
fclose (fid);

%-------------------------------------------------------------------------
