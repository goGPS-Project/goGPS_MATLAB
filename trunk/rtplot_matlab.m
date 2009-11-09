function rtplot_matlab (t, pos_R, check_on, check_off, check_pivot, check_cs, origin, ref, matrix)

% SYNTAX:
%   rtplot_matlab (t, pos_R, check_on, check_off, check_pivot, check_cs, origin, ref, matrix);
%
% INPUT:
%   t = survey time (t=1,2,...)
%   pos_R = ROVER assessed position (X,Y,Z)
%   check_on = boolean variable for satellite birth
%   check_off = boolean variable for satellite death
%   check_pivot = boolean variable for change of pivot
%   check_cs = boolean variable for cycle-slip
%   origin = boolean variable for displaying the axes origin or not
%   ref = reference path
%   matrix = adjacency matrix
%
% DESCRIPTION:
%   Real-time plot of the assessed ROVER path with respect to 
%   a reference path.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 pre-alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Media Center, Osaka City University, Japan
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

global pivot
global EST_M NORD_M
global pid p_max

%-------------------------------------------------------------------------------
% REFERENCE DISPLAY (DEMO)
%-------------------------------------------------------------------------------

subplot(2,3,[1 2 4 5])

if (t == 1)
    if (origin == 1)
        plot(0, 0, 'xm', 'LineWidth', 2);
    end
    
    if ~isempty(ref) % & ~isempty(matrix)
       hold on
       [EST_ref, NORD_ref, h_ref] = cart2plan(ref(:,1), ref(:,2), ref(:,3));
       EST_ref = EST_ref - EST_M;
       NORD_ref = NORD_ref - NORD_M;

       plot(EST_ref, NORD_ref, 'm', 'LineWidth', 2);
       for i = 1 : length(EST_ref)-1
           for j = i+1 : length(EST_ref)
               if (matrix(i,j) == 1)
                   plot([EST_ref(i),EST_ref(j)],[NORD_ref(i),NORD_ref(j)],'-m', 'LineWidth', 2);
               end
           end
       end
       hold off

    end
    xlabel('[m]'); ylabel('[m]');
    title('b:ok k:birth r:death m:pivot g:cycle-slip y:only-dyn');
    grid; axis equal
    hold on
end

%-------------------------------------------------------------------------------
% REAL-TIME DISPLAY (SIMULATED)
%-------------------------------------------------------------------------------

X = pos_R(1);
Y = pos_R(2);
Z = pos_R(3);

%conversion from cartesian to geodetic coordinates
[phi, lam, h] = cart2geod(X, Y, Z);

%conversion into metric coordinates
[EST, NORD] = geod2plan(phi,lam);

%permanent station origin translation
EST = EST - EST_M;
NORD = NORD - NORD_M;

%color choise
if (pivot == 0)
    pcol = 'y';
elseif check_cs
    pcol = 'g';
elseif check_pivot
    pcol = 'm';
elseif check_off
    pcol = 'r';
elseif check_on
    pcol = 'k';
else
    pcol = 'b';
end

%plotting  of path with Kalman filter
if (t <= p_max)
    pid(t) = plot(EST, NORD, ['.' pcol]);
else
    i = mod(t-1,p_max) + 1;
    set(pid(i), 'XData', EST, 'YData', NORD, 'Color', pcol);
end

%-------------------------------------------------------------------------------