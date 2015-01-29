function [table, nodes] = polyline_arcsClustering (dat_filename, cov_filename, filename, nodes, dN1, dN2, delta)

% SYNTAX:
%   [table, nodes] = polyline_arcsClustering (dat_filename, cov_filename, filename, nodes, dN1, dN2, delta);
%
% INPUT:
%   dat_filename = input data file name
%   cov_filename = input variances file name
%   filename = output file name
%   nodes = input coordinates of the nodes
%   dN1 = disregarded points at the beginning
%   dN2 = disregarded points at the end
%   angle_threshold = threshold on the angle between arcs [radiants]
%   dist_threshold = threshold on the distance between nodes [meters]
%
% OUTPUT:
%   table = output coordinates of the points with labels
%   nodes = output coordinates of the nodes
%
% DESCRIPTION:
%   Classify data belonging to each arc.

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

fid = fopen(dat_filename,'rt');      % open file
fgets(fid);                          % jump the header
data = fscanf(fid,'%d/%d/%d %d:%d:%f %d %f %f %f %f %f %f %f %f %f %f %4c %f %f %f %f %f %d %f',[28 inf])';
fclose(fid);
clear dat_filename

x0 = data(:,16);    % x coordinates
y0 = data(:,15);    % y coordinates
clear data

x0_start = x0(1);    % initial point
y0_start = y0(1);    % final point

x0 = x0 - x0_start;    % new x coordinates
y0 = y0 - y0_start;    % new y coordinates

nodes(:,1) = nodes(:,1) - x0_start;    % new x nodes
nodes(:,2) = nodes(:,2) - y0_start;    % new y nodes

%fprintf('%.3f\n',x0_start)
%fprintf('%.3f\n',y0_start)

%-----------------------------------------------------------
% loading the variances
%-----------------------------------------------------------

d = dir(cov_filename);

if ~isempty(d)

    fid = fopen(cov_filename,'rt');      % open file
    fgets(fid);                          % jump the header
    variances = fscanf(fid,'%f %f %f %f %f %f',[6 inf])';
    fclose(fid);
    clear cov_filename

    sigma2x = variances(:,4);
    sigma2y = variances(:,1);
    clear variances

else

    sigma2x = ones(size(x0));
    sigma2y = ones(size(y0));

end

%-----------------------------------------------------------
% labelling the data
%-----------------------------------------------------------

% boolean variable
loop = 1;

while loop

    %---------------------------------
    % plotting the data/nodes

    figure
    plot(x0,y0,'.k');
    hold on;
    plot(nodes(:,1),nodes(:,2),'*g'); % nodes
    plot(nodes(:,1),nodes(:,2),'-g'); % lines from one node to the other
    hold off
    %axis equal
    title('Bounding boxes');
    xlabel('East [m]');
    ylabel('North [m]');
    
    %---------------------------------

    NN = length(nodes)-1;        % define the number of segments

    % compute the angular coefficient
    Mo = zeros(NN,1);
    for i = 1 : NN;
        Mo(i) =( nodes(i+1,2) - nodes(i,2) ) / ( nodes(i+1,1) - nodes(i,1) );
    end

    % find the intersection point of the segment with the y-axis
    q = zeros(NN,1);
    for i = 1 : NN;
        q(i) = nodes(i,2) - Mo(i,1) * nodes(i,1); 
    end

    %---------------------------------------------------------------------

    % quadratic parameter and their solution for each segment (backward)
    a = zeros(NN,1);
    b = zeros(NN,1);
    c = zeros(NN,1);
    b0 = zeros(NN,1);
    yb0 = zeros(NN,1);

    for i = 1 : NN;

        a(i) = (1 + Mo(i)^2);
        b(i) = (2 * Mo(i) * q(i) - 2 * nodes(i,1) - 2 * Mo(i) * nodes(i,2));
        c(i) = nodes(i,1)^2 + q(i)^2 - 2 * q(i) * nodes(i,2) + nodes(i,2)^2 - delta^2;  % remember - delta^2 !!

        b0_1(i) = ( -b(i) + sqrt(b(i)^2 - 4*a(i)*c(i))) / (2*a(i) ); % the first quadratic solution
        b0_2(i) = ( -b(i) - sqrt(b(i)^2 - 4*a(i)*c(i))) / (2*a(i) ); % the second quadratic solution

        yb0_1(i) = Mo(i) * b0_1(i) + q(i);  % the corresponding y value for the first one
        yb0_2(i) = Mo(i) * b0_2(i) + q(i);

        % checking the distance
        dp1b0_1(i) = sqrt(( b0_1(i) - nodes(i+1,1) )^2 + (yb0_1(i) - nodes(i+1,2) )^2);
        dp1b0_2(i) = sqrt(( b0_2(i) - nodes(i+1,1) )^2 + (yb0_2(i) - nodes(i+1,2) )^2);

        % choosing the greatest distance
        if dp1b0_1(i) > dp1b0_2(i) ;
            b0(i) = b0_1(i) ;
            yb0(i) = yb0_1(i);
        else
            b0(i) = b0_2(i) ;
            yb0(i) = yb0_2(i);
        end        
    end

    clear a b c
    clear b0_1 b0_2
    clear yb0_1 yb0_2
    clear dp1b0_1 dp1b0_2

    %---------------------------------------------------------------------

    % quadratic parameter and their solution for each segment (forward)
    a = zeros(NN,1);
    b = zeros(NN,1);
    c = zeros(NN,1);
    b1 = zeros(NN,1);
    yb1 = zeros(NN,1);

    for i = 1 : NN;
        a(i) = (1 + Mo(i)^2);
        b(i) = (2 * Mo(i) * q(i) - 2 * nodes(i+1,1) - 2 * Mo(i) * nodes(i+1,2));
        c(i) = nodes(i+1,1)^2 + q(i)^2 - 2 * q(i) * nodes(i+1,2) + nodes(i+1,2)^2 - delta^2;                    % remember - delta^2 !!

        b1_1(i) = ( -b(i) + sqrt(b(i)^2 - 4*a(i)*c(i))) / (2*a(i) ); % the first quadratic solution
        b1_2(i) = ( -b(i) - sqrt(b(i)^2 - 4*a(i)*c(i))) / (2*a(i) ); % the second quadratic solution

        yb1_1(i) = Mo(i) * b1_1(i) + q(i);  % the coresponding y value for the first one
        yb1_2(i) = Mo(i) * b1_2(i) + q(i);

        % checking the distance
        dp0b1_1(i) = sqrt(( b1_1(i) - nodes(i , 1) )^2 +  (yb1_1(i) - nodes(i , 2) )^2);
        dp0b1_2(i) = sqrt(( b1_2(i) - nodes(i , 1) )^2 +  (yb1_2(i) - nodes(i , 2) )^2);

        %choosing the greatest distance
        if dp0b1_1(i) > dp0b1_2(i) ;
            b1(i) = b1_1(i) ;
            yb1(i) = yb1_1(i);
        else
            b1(i) = b1_2(i) ;
            yb1(i) = yb1_2(i);
        end
    end

    clear a b c
    clear b1_1 b1_2
    clear yb1_1 yb1_2
    clear dp0b1_1 dp0b1_2

    %---------------------------------------------------------------------

    % quadratic parameter and their solution for each segment
    a = zeros(NN,1);
    b = zeros(NN,1);
    c = zeros(NN,1);

    p0 = zeros(NN,1);
    for i = 1 : NN;
        p0(i) = yb0(i) + (1 / Mo(i)) * b0(i);
    end

    for i = 1 : NN;

        a(i) = (1 + (1 / -Mo(i))^2);
        b(i) = (2 * (1 / -Mo(i)) * p0(i) - 2 * b0(i) - 2 * (1 / -Mo(i)) * yb0(i));
        c(i) = b0(i)^2 + p0(i)^2 - 2 * p0(i) * yb0(i) + yb0(i)^2 - delta^2;       % remember - delta^2 !!

        b2_1(i) = ( -b(i) + sqrt(b(i)^2 - 4*a(i)*c(i))) / (2*a(i) ); % the first quadratic solution
        b2_2(i) = ( -b(i) - sqrt(b(i)^2 - 4*a(i)*c(i))) / (2*a(i) ); % the second quadratic solution

        yb2_1(i) = (1 / -Mo(i)) * b2_1(i) + p0(i);  % the corresponding y value for the first one
        yb2_2(i) = (1 / -Mo(i)) * b2_2(i) + p0(i);
    end

    clear a b c
    clear p0

    %---------------------------------------------------------------------

    % quadratic parameter and their solution for each segment
    a = zeros(NN,1);
    b = zeros(NN,1);
    c = zeros(NN,1);

    p1 = zeros(NN,1);
    for i = 1 : NN;
        p1(i) = yb1(i) + (1 / Mo(i)) * b1(i);
    end

    for i = 1 : NN;

        a(i) = (1 + (1 / -Mo(i))^2);
        b(i) = (2 * (1 / -Mo(i)) * p1(i) - 2 * b1(i) - 2 * (1 / -Mo(i)) * yb1(i));
        c(i) =  b1(i)^2 + p1(i)^2 - 2 * p1(i) * yb1(i) + yb1(i)^2 - delta^2;      % remember - delta^2 !!

        b3_1(i) = ( -b(i) + sqrt(b(i)^2 - 4*a(i)*c(i))) / (2*a(i)); % the first quadratic solution
        b3_2(i) = ( -b(i) - sqrt(b(i)^2 - 4*a(i)*c(i))) / (2*a(i)); % the second quadratic solution

        yb3_1(i) = (1 / -Mo(i)) * b3_1(i) + p1(i);  % the coresponding y value for the first one
        yb3_2(i) = (1 / -Mo(i)) * b3_2(i) + p1(i);
    end

    clear a b c
    clear p1

    %---------------------------------
    % hold on
    % plot(b0,yb0,'bo'); % plotting the backward nodes
    % plot(b1,yb1,'ro'); % plotting the forward nodes
    % plot([b2_1 , b2_2] , [yb2_1 , yb2_2] , 'g*'); % plotting the bounding nodes
    % plot([b3_1 , b3_2] , [yb3_1 , yb3_2] , 'g*'); % plotting the bounding nodes
    % hold off
    %---------------------------------

    hold on
    for i = 1 : NN
        plot([b2_1(i), b2_2(i)], [yb2_1(i), yb2_2(i)], '-m'); % plotting the bounding nodes
        plot([b2_2(i), b3_2(i)], [yb2_2(i), yb3_2(i)], '-m');
        plot([b3_2(i), b3_1(i)], [yb3_2(i), yb3_1(i)], '-m');
        plot([b3_1(i), b2_1(i)], [yb3_1(i), yb2_1(i)], '-m');
    end
    hold off

    %---------------------------------------------------------------------
    % PREPARE THE TABLE
    %---------------------------------------------------------------------

    N = length(x0);    % number of data

    N1 = 1+dN1;        % disregard initial data
    N2 = N-dN2;        % disregard final data
    Nt = N2-N1+1;      % new number of data

    segm = zeros(Nt,NN);         % one cell for each segment
    sum_segm = zeros(Nt,1);      % summing cell
    table = [x0(N1:N2), y0(N1:N2), sigma2x(N1:N2), sigma2y(N1:N2), segm , sum_segm];
    clear segm
    clear sum_segm

    %---------------------------------------------------------------------
    % define the rotation angle
    
    alpha = zeros(NN,1);
    for o = 1 : NN;
        alpha(o)= atan(Mo(o));
    end

    %---------------------------------------------------------------------
    % apply the rotation

    for i = 1 : NN

        T  = zeros(2,Nt);  % a matrix contained [x_rotated ; y0_rotated]
        B2_1 = zeros(2,1);
        B2_2 = zeros(2,1);
        B3_1 = zeros(2,1);
        B3_2 = zeros(2,1);

        %first translate in the origin
        T (1,:) = x0(N1:N2)' - b2_1(i);
        T (2,:) = y0(N1:N2)' - yb2_1(i);
        B2_1 (1) = b2_1(i) - b2_1(i);
        B2_1 (2) = yb2_1(i) - yb2_1(i);
        B2_2 (1) = b2_2(i) - b2_1(i);
        B2_2 (2) = yb2_2(i) - yb2_1(i);
        B3_1 (1) = b3_1(i) - b2_1(i);
        B3_1 (2) = yb3_1(i) - yb2_1(i);
        B3_2 (1) = b3_2(i) - b2_1(i);
        B3_2 (2) = yb3_2(i) - yb2_1(i);

        %then appply the rotation
        R = [cos(alpha(i)) sin(alpha(i)); -sin(alpha(i)) cos(alpha(i))];  %rotation matrix
        T_rot = R * T;
        B2_1_rot = R * B2_1;
        B2_2_rot = R * B2_2;
        B3_1_rot = R * B3_1;
        B3_2_rot = R * B3_2; %#ok<NASGU>

        clear R T
        clear B2_1 B2_2 B3_1 B3_2
                
        % to know which point belongs to which segment
        T_rot = T_rot';
        for j = 1:length(T_rot)

            right = max(B2_1_rot(1), B3_1_rot(1));
            left = min(B2_1_rot(1), B3_1_rot(1));
            up = max(B2_1_rot(2), B2_2_rot(2));
            down = min(B2_1_rot(2), B2_2_rot(2));

            if (T_rot(j,1) >= left) & (T_rot(j,1) <= right) & (T_rot(j,2) >= down) & (T_rot(j,2) <= up)
                table(j,4+i) = 1;
            else
                table(j,4+i) = 0;
            end
        end
    end

    clear b0 b1
    clear yb0 yb1
    clear b2_1 b2_2 b3_1 b3_2
    clear yb2_1 yb2_2 yb3_1 yb3_2
    clear B2_1_rot B2_2_rot
    clear B3_1_rot B3_2_rot
    clear T_rot
    clear left right
    clear up down

    %---------------------------------------------------------------------    
    % check whether there are points that are classified twice

    table(:,end) = sum(table(:,5:end-1),2);

    for i = 1 : Nt
        if table(i,end) > 1

            % d_point_line = ones(NN,1).*10e4;        % a trick to consider minimum distances correctly
            % j = find(table(i,5:end-1)==1);
            % d_point_line(j) = abs((nodes(j+1,1) - nodes(j,1)).*(nodes(j,2) - table(i,2)) - (nodes(j,1) - table(i,1)).*(nodes(j+1,2) - nodes(j,2)))./(sqrt((nodes(j+1,1) - nodes(j,1)).^2+(nodes(j+1,2) - nodes(j,2)).^2));
            % k = find(d_point_line == min(d_point_line));

            d_point_nodes = ones(NN,1).*10e4;        % a trick to consider minimum distances correctly
            j = find(table(i,5:end-1)==1);
            d_point_nodes(j) = sqrt((table(i,1)-nodes(j,1)).^2 + (table(i,2)-nodes(j,2)).^2) + sqrt((table(i,1)-nodes(j+1,1)).^2 + (table(i,2)-nodes(j+1,2)).^2);
            k = find(d_point_nodes == min(d_point_nodes));

            table(i,5:end-1) = 0;
            table(i,4+k) = 1;
        end
    end

    table(:,end) = 0;   % the field will be used for labelling

    %---------------------------------------------------------------------    
    % labelling

    for i = 1 : Nt
        for j = 1 : NN;
            if find(table(i,4+j)==1)
                table(i,end) = j;
            end
        end
    end

    %---------------------------------------------------------------------    
    % check whether there are segments containing less than two data

    h = hist(table(:,end),1:NN);
    pos = find(h<2);
    clear h

    if isempty(pos)
        loop = 0;
    else
        nodes(pos,:) = [];
    end

end

%-------------------------------------------------------------------------
% check whether there are points that are not classified

pos = table(:,end) == 0;
table(pos,:) = [];  % delete data with the label 0
clear pos

%-------------------------------------------------------------------------
% debugging
%-------------------------------------------------------------------------

% for i = 1 : 60
%     pos = find(label == i);
%     plot(table(pos,1),table(pos,2),'.b');
%     hold on
%     plot(nodes(i,1), nodes(i,2), '.r')
%     plot(nodes(i+1,1), nodes(i+1,2), '.r')
%     hold off
%     pause
% end

% for i = 1 : 60
%     pos = find(label == i);
%     plot(table(pos,1),table(pos,2),'.b');
%     hold on
%     plot(nodes(i,1), nodes(i,2), '.r')
%     plot(nodes(i+1,1), nodes(i+1,2), '.r')
%     
%     pos = find(label == i+1);
%     plot(table(pos,1),table(pos,2),'.c');
%     hold on
%     plot(nodes(i+1,1), nodes(i+1,2), '.m')
%     plot(nodes(i+2,1), nodes(i+2,2), '.m')
%     
%     hold off
%     pause
% end

%-------------------------------------------------------------------------
% print the result in a txt file

nodes(:,1) = nodes(:,1) + x0_start;
nodes(:,2) = nodes(:,2) + y0_start;

table(:,1) = table(:,1) + x0_start;
table(:,2) = table(:,2) + y0_start;
data = [table(:,1:4), table(:,end)];

fid = fopen(filename,'wt');
fprintf(fid,'%12.4f %12.4f %7.4f %7.4f %5.0d\n',data');
fclose (fid);

%-------------------------------------------------------------------------
