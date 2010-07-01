function rtplot_matlab (t, pos_R, pos_M, check_on, check_off, check_pivot, check_cs, origin, ref, matrix, flag_amb)

% SYNTAX:
%   rtplot_matlab (t, pos_R, pos_M, check_on, check_off, check_pivot, check_cs, origin, ref, matrix, flag_amb);
%
% INPUT:
%   t = survey time (t=1,2,...)
%   pos_R = ROVER position (X,Y,Z)
%   pos_M = MASTER station position (X,Y,Z)
%   check_on = boolean variable for satellite birth
%   check_off = boolean variable for satellite death
%   check_pivot = boolean variable for change of pivot
%   check_cs = boolean variable for cycle-slip
%   origin = boolean variable for displaying the axes origin or not
%   ref = reference path
%   matrix = adjacency matrix
%   flag_amb = boolean variable for displaying ambiguities or not
%
% DESCRIPTION:
%   Real-time plot of the assessed ROVER path with respect to
%   a reference path.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 alpha
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

global pivot
global msid pid p_max
global EST_O NORD_O

%-------------------------------------------------------------------------------
% REFERENCE DISPLAY (DEMO)
%-------------------------------------------------------------------------------

if (nargin == 11 & flag_amb)
    if (t == 1)
        figure('Units','normalized','Position',[0 0 1 1])
    end
    subplot(5,3,[1 2 3 4 5 6])
else
    subplot(2,3,[1 2 4 5])
end

if (t == 1)
    xlabel('[m]'); ylabel('[m]');
    title('b:ok k:birth r:death m:pivot g:cycle-slip y:only-dyn');
    grid; axis equal
    hold on
    
    % fixing of the origin O
    %EST_O = 0; %NORD_O = 0;
    %EST_O = EST_M; %NORD_O = NORD_M;
    [EST_O, NORD_O] = cart2plan(pos_R(1), pos_R(2), pos_R(3));
    
    if ~isempty(ref) % & ~isempty(matrix)
        [EST_ref, NORD_ref, h_ref] = cart2plan(ref(:,1), ref(:,2), ref(:,3)); %#ok<NASGU>
        
        plot(EST_ref-EST_O, NORD_ref-NORD_O, 'm', 'LineWidth', 2);
        for i = 1 : length(EST_ref)-1
            for j = i+1 : length(EST_ref)
                if (matrix(i,j) == 1)
                    plot([EST_ref(i)-EST_O,EST_ref(j)-EST_O],[NORD_ref(i)-NORD_O,NORD_ref(j)-NORD_O],'-m', 'LineWidth', 2);
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
        msid = plot(EST_M-EST_O, NORD_M-NORD_O, 'xm', 'LineWidth', 2);
    else
        set(msid, 'XData', EST_M-EST_O, 'YData', NORD_M-NORD_O);
    end
end

%-------------------------------------------------------------------------------
% REAL-TIME DISPLAY
%-------------------------------------------------------------------------------

X = pos_R(1);
Y = pos_R(2);
Z = pos_R(3);

%conversion into metric coordinates
[EST, NORD] = cart2plan(X, Y, Z);

%color choice
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

%Kalman filter path plot
if (t <= p_max)
    pid(t) = plot(EST-EST_O, NORD-NORD_O, ['.' pcol]);
else
    i = mod(t-1,p_max) + 1;
    set(pid(i), 'XData', EST-EST_O, 'YData', NORD-NORD_O, 'Color', pcol);
end

%-------------------------------------------------------------------------------