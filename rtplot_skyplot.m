function rtplot_skyplot (az, el, obs, pivot)

% SYNTAX:
%   rtplot_skyplot (az, el, obs, pivot);
%
% INPUT:
%   az  = azimuth               [degrees]
%   el  = elevation             [degrees]
%   obs = kind of observation
%            0 = not used
%           +1 = code & phase
%           -1 = only code
%   pivot = pivot satellite
%
% DESCRIPTION:
%   Real time skyplot.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 pre-alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini*
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
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

global cutoff

sat = find(el > 0);
az = az(sat);
el = el(sat);
obs = obs(sat);

%-------------------------------------------------------------------------------

subplot(2,3,3)

%theta =   0° --> NORTH
%theta =  90° --> EAST
%theta = 180° --> SOUTH
%theta = -90° --> WEST

theta = (90-az) * pi/180;
rho = 90-el;

x = rho .* cos(theta);
y = rho .* sin(theta);
D = squareform(pdist([x y])) + 180*eye(length(rho));

X = x(:,ones(length(x),1));
Y = y(:,ones(length(y),1));
dX = X-X';
dY = Y-Y';

polar_goGPS(theta,rho,obs,'.',cutoff);

hold on
for i = 1: length(sat)
   if (sat(i) == pivot)
      p = plot(rho(i)*cos(theta(i)),rho(i)*sin(theta(i)),'ok');
      set(p,'MarkerSize',8,'LineWidth',1);
   end
   if ~isempty(D)
      [minD, d] = min(D(:,i));
   else
      d = 1;
   end
   alpha = atan2(dY(d,i),dX(d,i));
   %[sat(i) sat(d) minD alpha*180/pi]
   if (cos(alpha) < 0)
      h = text(rho(i)*cos(theta(i))-10*cos(alpha)-6,rho(i)*sin(theta(i))-10*sin(alpha),num2str(sat(i)));
   else
      if (sat(i) < 10)
         h = text(rho(i)*cos(theta(i))-10*cos(alpha),rho(i)*sin(theta(i))-10*sin(alpha),num2str(sat(i)));
      else
         h = text(rho(i)*cos(theta(i))-10*cos(alpha)-2,rho(i)*sin(theta(i))-10*sin(alpha),num2str(sat(i)));
      end
   end
   %h = text(rho(i)*cos(theta(i)),rho(i)*sin(theta(i)),num2str(sat(i)));
   set(h,'Color','k')
   set(h,'FontWeight','bold')
end;
title('sky plot')
hold off

%-------------------------------------------------------------------------------