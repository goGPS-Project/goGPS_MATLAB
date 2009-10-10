function polar_goGPS(theta,rho,flag,line_style,cutoff)
%POLAR  Polar coordinate plot.
%   POLAR(THETA, RHO) makes a plot using polar coordinates of
%   the angle THETA, in radians, versus the radius RHO.
%   POLAR(THETA,RHO,S) uses the linestyle specified in string S.
%   See PLOT for a description of legal linestyles.
%
%   See also PLOT, LOGLOG, SEMILOGX, SEMILOGY.
%
%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.17 $  $Date: 1998/09/30 15:25:05 $
%
%   Adapted by Mirko Reguzzoni, Eugenio Realini

%---------------------------------------------------------------

% get hold state
cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;

% get x-axis text color so grid is in same color
tc = get(cax,'xcolor');
ls = get(cax,'gridlinestyle');

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
fWeight = get(cax, 'DefaultTextFontWeight');
fUnits  = get(cax, 'DefaultTextUnits');
set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
         'DefaultTextFontName',   get(cax, 'FontName'), ...
         'DefaultTextFontSize',   get(cax, 'FontSize'), ...
         'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
         'DefaultTextUnits','data')

% only do grids if hold is off
if ~hold_state

% make a radial grid
    hold on;
    maxrho = max(abs(rho(:)));
    hhh=plot([-maxrho -maxrho maxrho maxrho],[-maxrho maxrho maxrho -maxrho]);
    set(gca,'dataaspectratio',[1 1 1],'plotboxaspectratiomode','auto')
    v = [get(cax,'xlim') get(cax,'ylim')];
    ticks = sum(get(cax,'ytick')>=0);
    delete(hhh);

% check radial limits and ticks
    rmin = 0;
    rmax = 90;
    rinc = 15;
    rticks = (rmax-rmin)/rinc + 1;

% define a circle
    th = 0:pi/50:2*pi;
    xunit = cos(th);
    yunit = sin(th);

% now really force points on x/y axes to lie on them exactly
    inds = 1:(length(th)-1)/4:length(th);
    xunit(inds(2:2:4)) = zeros(2,1);
    yunit(inds(1:2:5)) = zeros(3,1);

% plot background if necessary
    if ~isstr(get(cax,'color')),
       patch('xdata',xunit*rmax,'ydata',yunit*rmax, ...
             'edgecolor',tc,'facecolor',get(gca,'color'),...
             'handlevisibility','off');
    end

% draw radial circles
    c82 = cos(82*pi/180);
    s82 = sin(82*pi/180);
    rinc = 15;
    for i=0:rinc:90
        hhh = plot(xunit*i,yunit*i,ls,'color',tc,'linewidth',1,...
                   'handlevisibility','off');
%        text((i+rinc/20)*c82,(i+rinc/20)*s82, ...
%            ['  ' num2str(i)],'verticalalignment','bottom',...
%            'handlevisibility','off')
    end
    set(hhh,'linestyle','-') % Make outer circle solid

    hhh = plot(xunit*(90-cutoff),yunit*(90-cutoff),ls,'color','r','linewidth',1,...
               'handlevisibility','off');
    set(hhh,'linestyle','-') % Cutoff

% plot spokes
    th = (1:4)*pi/4;
    cst = cos(th); snt = sin(th);
    cs = [-cst; cst];
    sn = [-snt; snt];
    plot(rmax*cs,rmax*sn,ls,'color',tc,'linewidth',1,...
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

% set view to 2-D
    view(2);
% set axis limits
    axis(rmax*[-1 1 -1.15 1.15]);
end

% Reset defaults.
set(cax, 'DefaultTextFontAngle',  fAngle , ...
         'DefaultTextFontName',   fName , ...
         'DefaultTextFontSize',   fSize, ...
         'DefaultTextFontWeight', fWeight, ...
         'DefaultTextUnits',      fUnits );

% transform data to Cartesian coordinates.
xx = rho.*cos(theta);
yy = rho.*sin(theta);

% plot data on top of grid
f = find(flag == 1);
p = plot(xx(f),yy(f),[line_style 'b']);
set(p,'MarkerSize',15);
hold on
f = find(flag == -1);
p = plot(xx(f),yy(f),[line_style 'g']);
set(p,'MarkerSize',15);
f = find(flag == 0);
p = plot(xx(f),yy(f),[line_style 'm']);
set(p,'MarkerSize',15);
hold off

if ~hold_state
    set(gca,'dataaspectratio',[1 1 1]), axis off; set(cax,'NextPlot',next);
end
set(get(gca,'xlabel'),'visible','on')
set(get(gca,'ylabel'),'visible','on')
