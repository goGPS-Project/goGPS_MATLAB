function rtplot_skyplot (t, az, el, obs, pivot, Eph, SP3)

% SYNTAX:
%   rtplot_skyplot (t, az, el, obs, pivot, Eph, SP3);
%
% INPUT:
%   t   = survey time (t=1,2,...)
%   az  = azimuth     [degrees]
%   el  = elevation   [degrees]
%   obs = observation type
%            0 = not used
%           +1 = code & phase
%           -1 = only code
%   pivot = pivot satellite
%   Eph = matrix containing 33 navigation parameters for each satellite
%   SP3 = structure containing precise ephemeris data
%
% DESCRIPTION:
%   Real time skyplot.

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

global cutoff
global satid labid pivid

num_sat = length(el);

if (isempty(SP3))
    eph_avail = Eph(30,:);
    sys = Eph(31,:);
    prn = Eph(1,:);
else
    eph_avail = SP3.avail';
    sys = SP3.sys';
    prn = SP3.prn';
end

%----------------------------------------------------------------------------------------------
% SKY-PLOT BACKGROUND
%----------------------------------------------------------------------------------------------

% location on the screen
subplot(2,3,3)

if (t == 1)

    % define radial limits and step
    rmin = 0; %#ok<NASGU>
    rmax = 90;
    rinc = 15;

    % make a bounding box
    hhh = plot([-rmax -rmax rmax rmax],[-rmax rmax rmax -rmax]);
    set(gca,'dataaspectratio',[1 1 1],'plotboxaspectratiomode','auto')
    delete(hhh);
    hold on

    % define a circle
    th = 0:pi/50:2*pi;
    xunit = cos(th);
    yunit = sin(th);

    % plot white circle background
    patch('xdata',xunit*rmax,'ydata',yunit*rmax, ...
        'edgecolor','k','facecolor','w', ...
        'handlevisibility','off');

    % draw radial circles
    for i=0:rinc:90
        hhh = plot(xunit*i,yunit*i,':','color','k', ...
            'linewidth',1,'handlevisibility','off');
    end
    set(hhh,'linestyle','-') % make outer circle solid

    hhh = plot(xunit*(90-cutoff),yunit*(90-cutoff),':','color','r', ...
        'linewidth',1,'handlevisibility','off');
    set(hhh,'linestyle','-') % make cutoff circle in red

    % plot spokes
    th = (1:4)*pi/4;
    cst = cos(th);
    snt = sin(th);
    cs = [-cst; cst];
    sn = [-snt; snt];
    plot(rmax*cs,rmax*sn,':','color','k','linewidth',1,...
        'handlevisibility','off')

    % annotate spokes in degrees
    rt = 1.1*rmax;
    loc{1} = char('NE');  loc{2} = char('SW');
    loc{3} = char('N');   loc{4} = char('S');
    loc{5} = char('NW');  loc{6} = char('SE');
    loc{7} = char('W');   loc{8} = char('E');
    for i = 1:length(th)
        text(rt*cst(i),rt*snt(i),loc{2*i-1},...
            'horizontalalignment','center',...
            'handlevisibility','off');
        text(-rt*cst(i),-rt*snt(i),loc{2*i},...
            'horizontalalignment','center',...
            'handlevisibility','off')
    end

    % set axis limits
    axis(rmax*[-1 1 -1.15 1.15]);

    % no background
    set(gca,'dataaspectratio',[1 1 1])
    axis off

    % write the title
    title('sky plot')

end

%----------------------------------------------------------------------------------------------
% COMPUTATION OF POLAR/CARTESIAN COORDINATES
%----------------------------------------------------------------------------------------------

%theta =   0° --> NORTH
%theta =  90° --> EAST
%theta = 180° --> SOUTH
%theta = -90° --> WEST

% polar coordinates
theta = (90-az) * pi/180;
rho = 90-el;

% cartesian coordinates
x = rho .* cos(theta);
y = rho .* sin(theta);

%----------------------------------------------------------------------------------------------
% PLOT OF SATELLITE POINTS
%----------------------------------------------------------------------------------------------

% visible satellites
sat = [];

for i = 1 : num_sat

    if (el(i) > 0)

        sat = [sat; i];
        
        if (satid(i) == 0)
            
            satid(i) = plot(x(i), y(i), '.');
            set(satid(i), 'MarkerSize', 15);
            
            switch sys(i)
                case 'G' %GPS
                    set(satid(i), 'Color', 'b');
                case 'R' %GLONASS
                    set(satid(i), 'Color', 'g');
                case 'E' %Galileo
                    set(satid(i), 'Color', 'y');
                case 'C' %BeiDou
                    set(satid(i), 'Color', 'c');
                case 'J' %QZSS
                    set(satid(i), 'Color', 'm');
            end
        else
            set(satid(i), 'XData', x(i), 'YData', y(i));
        end

        if (obs(i) == 1)
            set(satid(i), 'Marker', '.');
            set(satid(i), 'MarkerSize', 15);
        elseif (obs(i) == -1)
            set(satid(i), 'Marker', 'o');
            set(satid(i), 'MarkerSize', 6);
        else % if (obs(i) == 0)
            set(satid(i), 'Marker', '^');
            set(satid(i), 'MarkerSize', 6);
        end

        if (i == pivot)
            if (pivid ~= 0)
                delete(pivid);
            end
            pivid = plot(x(i), y(i), 'ok');
            set(pivid, 'MarkerSize', 8, 'LineWidth', 1);
        end

    elseif (satid(i) > 0)

        delete(satid(i))
        satid(i) = 0;

    end
end

%----------------------------------------------------------------------------------------------
% SELECTION OF SATELLITE INFORMATION
%----------------------------------------------------------------------------------------------

rho = rho(sat); %#ok<NASGU>
theta = theta(sat); %#ok<NASGU>

x = x(sat);
y = y(sat);

%----------------------------------------------------------------------------------------------
% PLOT OF SATELLITE LABELS
%----------------------------------------------------------------------------------------------

% delete previous labels and re-initialize them
delete(labid(labid > 0));
labid = zeros(num_sat,1);

% identify GNSS systems and prepare system IDs and PRNs/slot numbers
[~, idx1, idx2] = intersect(eph_avail, sat);
sat = sat(idx2);
sys = char(sys(idx1));
prn = prn(idx1);

try
    % if Statistics Toolbox is installed
    D = squareform(pdist([x y])) + 180*eye(length(sat));
catch
    % if Statistics Toolbox is not installed
    D = [];
end

X = x(:,ones(length(sat),1));
Y = y(:,ones(length(sat),1));
dX = X-X';
dY = Y-Y';

for i = 1 : length(sat)

    if ~isempty(D)
        [null, d] = min(D(:,i)); %#ok<ASGLU>
    else
        d = 1;
    end

    alpha = atan2(dY(d,i),dX(d,i));
    %[sat(i) sat(d) minD alpha*180/pi]

    if (cos(alpha) < 0)
        labid(sat(i)) = text(x(i)-10*cos(alpha)-6,y(i)-10*sin(alpha),[sys(i) sprintf('%d',prn(i))]);
    else
        if (sat(i) < 10)
            labid(sat(i)) = text(x(i)-10*cos(alpha),y(i)-10*sin(alpha),[sys(i) sprintf('%d',prn(i))]);
        else
            labid(sat(i)) = text(x(i)-10*cos(alpha)-2,y(i)-10*sin(alpha),[sys(i) sprintf('%d',prn(i))]);
        end
    end
    %labid(sat(i)) = text(x(i),y(i),num2str(sat(i)));
    set(labid(sat(i)),'Color','k')
    set(labid(sat(i)),'FontSize',9)
    set(labid(sat(i)),'FontWeight','bold')
end
