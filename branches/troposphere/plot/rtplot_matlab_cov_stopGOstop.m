function rtplot_matlab_cov_stopGOstop (t, pos_R, pos_M, covpos_R, P1, P2, origin, ref, matrix, flag_dyn, flag_amb)

% SYNTAX:
%   rtplot_matlab_cov_stopGOstop (t, pos_R, pos_M, covpos_R, P1, P2, origin, ref, matrix, flag_dyn, flag_amb);
%
% INPUT:
%   t = survey time (t=1,2,...)
%   pos_R = ROVER position (X,Y,Z)
%   pos_M = MASTER station position (X,Y,Z)
%   covpos_R = ROVER position covariance matrix
%   P1 = initial point of the straight line
%   P2 = final point of the straight line
%   origin = boolean variable for displaying the axes origin or not
%   ref = reference path
%   matrix = adjacency matrix
%   flag_dyn = integer variable for dynamical model
%   flag_amb = boolean variable for displaying ambiguities or not
%
% DESCRIPTION:
%   Real-time plot of the assessed ROVER path with respect to
%   a reference path.

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

global msid pid
global EAST_O NORTH_O
global x_circle id_ellipse

%-------------------------------------------------------------------------------
% REFERENCE DISPLAY (DEMO)
%-------------------------------------------------------------------------------

if (nargin == 12 & flag_amb)
    if (t == 1)
        figure('Units','normalized','Position',[0 0 1 1])
    end
    subplot(5,3,[1 2 3 4 5 6])
else
    subplot(2,3,[1 2 4 5])
end

if (t == 1)
    xlabel('[m]'); ylabel('[m]');
    title('stop(blue) - go(green)');
    grid; axis equal
    hold on

    % fixing of the origin O
    %EAST_O = 0; %NORTH_O = 0;
    %EAST_O = EST_M; %NORTH_O = NORD_M;
    [EAST_O, NORTH_O] = cart2plan(pos_R(1), pos_R(2), pos_R(3));

    if ~isempty(ref) % & ~isempty(matrix)
        [EAST_ref, NORTH_ref, h_ref] = cart2plan(ref(:,1), ref(:,2), ref(:,3)); %#ok<NASGU>

        plot(EAST_ref-EAST_O, NORTH_ref-NORTH_O, 'm', 'LineWidth', 2);
        for i = 1 : length(EAST_ref)-1
            for j = i+1 : length(EAST_ref)
                if (matrix(i,j) == 1)
                    plot([EAST_ref(i)-EAST_O,EAST_ref(j)-EAST_O],[NORTH_ref(i)-NORTH_O,NORTH_ref(j)-NORTH_O],'-m', 'LineWidth', 2);
                end
            end
        end
    end
end

if (origin & sum(abs(pos_M)) ~= 0)

    %Master position in UTM coordinates (East, North, h)
    [EST_M, NORD_M] = cart2plan(pos_M(1), pos_M(2), pos_M(3));

    %master station plot
    if (isempty(msid))
        msid = plot(EST_M-EAST_O, NORD_M-NORTH_O, 'xm', 'LineWidth', 2);
    else
        set(msid, 'XData', EST_M-EAST_O, 'YData', NORD_M-NORTH_O);
    end
end

%-------------------------------------------------------------------------------
% REAL-TIME DISPLAY
%-------------------------------------------------------------------------------

X = pos_R(1);
Y = pos_R(2);
Z = pos_R(3);

%conversion into metric coordinates
[EAST, NORTH] = cart2plan(X, Y, Z);

%plot the estimated points (static in blue, kinematic in green) 
if (flag_dyn == 1) && (pid(1) == 0)
    pcol = 'b';
    pid(1) = plot(EAST-EAST_O, NORTH-NORTH_O, ['.' pcol]);
elseif (flag_dyn == 1) && (pid(1) ~= 0)
    pcol = 'b';
    set(pid(1), 'XData', EAST-EAST_O, 'YData', NORTH-NORTH_O, 'Color', pcol);
elseif (flag_dyn == 2)
    pcol = 'g';
    plot(EAST-EAST_O, NORTH-NORTH_O, ['.' pcol]);
elseif (flag_dyn == 3) && (pid(2) == 0)
    pcol = 'b';
    pid(2) = plot(EAST-EAST_O, NORTH-NORTH_O, ['.' pcol]);
elseif (flag_dyn == 3) && (pid(2) ~= 0)
    pcol = 'b';
    set(pid(2), 'XData', EAST-EAST_O, 'YData', NORTH-NORTH_O, 'Color', pcol);
end

%draw the straight line
if (t == 1)
    pid(3) = plot([P1(1),P2(1)]-EAST_O, [P1(2),P2(2)]-NORTH_O, '.-r');
    set(pid(3), 'LineWidth', 1.5);
else
    set(pid(3), 'XData', [P1(1),P2(1)]-EAST_O, 'YData', [P1(2),P2(2)]-NORTH_O, 'Color', 'r');
end

%make the straight line visible or not
if (P1(1) == P2(1)) & (P1(2) == P2(2))
    set(pid(3), 'Visible', 'off');
else
    set(pid(3), 'Visible', 'on');
end

%put the straight line in foreground
pid_all = get(gca, 'children');
set(gca, 'children', [pid(3); setdiff(pid_all, pid(3))]);

%-------------------------------------------------------------------------------

%covariance propagation
covpos_R = global2localCov(covpos_R, pos_R);

T = chol(covpos_R(1:2,1:2));        % Cholesky decomposition
for j = 1 : size(x_circle,1)        % ellipse computation
    x_ellipse(j,:) = x_circle(j,:) * T + [EAST-EAST_O, NORTH-NORTH_O];
end

%ellipse plot
if ~isempty(id_ellipse), delete(id_ellipse); end;
id_ellipse = plot(x_ellipse(:,1),x_ellipse(:,2), pcol);

%-------------------------------------------------------------------------------