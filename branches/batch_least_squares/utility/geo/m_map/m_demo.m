function m_demo(index);
% M_DEMO  Demonstration program showing various maps in M_Map package
%         Dig into this to look for examples of things you want to do.
%
%         M_DEMO runs all the demos.
%         M_DEMO(NUM) runs examples NUM (1<=NUM<=10), pausing between examples.
%
%         Some demos may require you to install GSHHS or TerrainBase datafiles
%         (see documentation)

% Rich Pawlowicz (rich@ocgy.ubc.ca) 7/May/1997
% (thanks to Art Newhall for putting these examples into an m-file,
% and the Chuck Denham for enhancing the interface)
%
% 27/July/98 - more examples.
% 17/Aug/98     "
% 15/Nov/98  - another example, better interface.
% 23/Dec/98  - another example.

%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.


N_EXAMPLES=13;

if nargin==0,
 index=1:N_EXAMPLES;
end; 

for i=index,

clf;
switch i,

  case 1,

    m_proj('ortho','lat',48','long',-123');
    m_coast('patch','r');
    m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',[-80:40:80]);
    xlabel('Orthographic Projection','visible','on');

  case 2,

    m_proj('lambert','long',[-160 -40],'lat',[30 80]);
    m_coast('patch',[1 .85 .7]);
    m_elev('contourf',[500:500:6000]);
    m_grid('box','fancy','tickdir','in');
    colormap(flipud(copper));
    xlabel('Conic Projection of North America with elevations','visible','on');
   
  case 3,

    m_proj('stereographic','lat',90,'long',30,'radius',25);
    m_elev('contour',[-3500:1000:-500],'edgecolor','b');
    m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linestyle','-');
    m_coast('patch',[.7 .7 .7],'edgecolor','r');
    xlabel('Polar Stereographic Projection with bathymetry','visible','on');

  case 4,
  
    subplot(211);
    Slongs=[-100 0;-75 25;-5 45; 25 145;45 100;145 295;100 290];
    Slats= [  8 80;-80  8; 8 80;-80   8; 8  80;-80   0;  0  80];
    for l=1:7,
     m_proj('sinusoidal','long',Slongs(l,:),'lat',Slats(l,:));
     m_grid('fontsize',6,'xticklabels',[],'xtick',[-180:30:360],...
            'ytick',[-80:20:80],'yticklabels',[],'linestyle','-','color',[.9 .9 .9]);
     m_coast('patch','g');
    end;
    xlabel('Interrupted Sinusoidal Projection of World Oceans');
    % In order to see all the maps we must undo the axis limits set by m_grid calls:
    set(gca,'xlimmode','auto','ylimmode','auto');

    subplot(212);
    Slongs=[-100 43;-75 20; 20 145;43 100;145 295;100 295];
    Slats= [  0  90;-90  0;-90   0; 0  90;-90   0;  0  90];
    for l=1:6,
     m_proj('mollweide','long',Slongs(l,:),'lat',Slats(l,:));
     m_grid('fontsize',6,'xticklabels',[],'xtick',[-180:30:360],...
            'ytick',[-80:20:80],'yticklabels',[],'linestyle','-','color','k');
     m_coast('patch',[.6 .6 .6]);
    end;
    xlabel('Interrupted Mollweide Projection of World Oceans');
    set(gca,'xlimmode','auto','ylimmode','auto');

  case 5,
  
    %% Nice looking data
    [lon,lat]=meshgrid([-136:2:-114],[36:2:54]);
    u=sin(lat/6);
    v=sin(lon/6);

    m_proj('oblique','lat',[56 30],'lon',[-132 -120],'aspect',.8);

    subplot(121);
    m_coast('patch',[.9 .9 .9],'edgecolor','none');
    m_grid('tickdir','out','yaxislocation','right',...
	   'xaxislocation','top','xlabeldir','end','ticklength',.02);
    hold on;
    m_quiver(lon,lat,u,v);
    xlabel('Simulated surface winds');

    subplot(122);
    m_coast('patch',[.9 .9 .9],'edgecolor','none');
    m_grid('tickdir','out','yticklabels',[],...
	   'xticklabels',[],'linestyle','none','ticklength',.02);
    hold on;
    [cs,h]=m_contour(lon,lat,sqrt(u.*u+v.*v));
    clabel(cs,h,'fontsize',8);
    xlabel('Simulated something else');

  case 6,
  
    % Plot a circular orbit
    lon=[-180:180];
    lat=atan(tan(60*pi/180)*cos((lon-30)*pi/180))*180/pi;

    m_proj('miller','lat',82);
    m_coast('color',[0 .6 0]);
    m_line(lon,lat,'linewidth',3,'color','r');
    m_grid('linestyle','none','box','fancy','tickdir','out');


  case 7,
  
    m_proj('lambert','lon',[-10 20],'lat',[33 48]);
    m_tbase('contourf');
    m_grid('linestyle','none','tickdir','out','linewidth',3);
    colormap(jet);

  case 8,
  
    m_vec;
    
  case 9,

    % Example showing the default coastline and all of the GSHHS coastlines.

    axes('position',[.35 .6 .37 .37]);
    m_proj('albers equal-area','lat',[40 60],'long',[-90 -50],'rect','on');
    m_coast('patch',[0 1 0]);
    m_grid('linestyle','none','linewidth',2,'tickdir','out','xaxisloc','top','yaxisloc','right');
    m_text(-69,41,'Standard coastline','color','r','fontweight','bold');

    axes('position',[.09 .5 .37 .37]);
    m_proj('albers equal-area','lat',[40 54],'long',[-80 -55],'rect','on');
    m_gshhs_c('patch',[.2 .8 .2]);
    m_grid('linestyle','none','linewidth',2,'tickdir','out','xaxisloc','top');
    m_text(-80,52.5,'GSHHS\_C (crude)','color','m','fontweight','bold','fontsize',14);

    axes('position',[.13 .2 .37 .37]);
    m_proj('albers equal-area','lat',[43 48],'long',[-67 -59],'rect','on');
    m_gshhs_l('patch',[.4 .6 .4]);
    m_grid('linestyle','none','linewidth',2,'tickdir','out');
    m_text(-66.5,43.5,'GSHHS\_L (low)','color','m','fontweight','bold','fontsize',14);

    axes('position',[.35 .05 .37 .37]);
    m_proj('albers equal-area','lat',[45.8 47.2],'long',[-64.5 -62],'rect','on');
    m_gshhs_i('patch',[.5 .6 .5]);
    m_grid('linestyle','none','linewidth',2,'tickdir','out','yaxisloc','right');
    m_text(-64.4,45.9,'GSHHS\_I (intermediate)','color','m','fontweight','bold','fontsize',14);

    axes('position',[.55 .23 .37 .37]);
    m_proj('albers equal-area','lat',[46.375 46.6],'long',[-64.2 -63.7],'rect','on');
    m_gshhs_h('patch',[.6 .6 .6]);
    m_grid('linestyle','none','linewidth',2,'tickdir','out','xaxisloc','top','yaxisloc','right');
    m_text(-64.18,46.58,'GSHHS\_H (high)','color','m','fontweight','bold','fontsize',14);

  case 10,
  
    % Example showing a trackline plot

    clf
    m_proj('UTM','long',[-72 -68],'lat',[40 44]);
    m_gshhs_i('color','k');
    m_grid('box','fancy','tickdir','in');

    % fake up a trackline
    lons=[-71:.1:-67];
    lats=60*cos((lons+115)*pi/180);
    dates=datenum(1997,10,23,15,1:41,zeros(1,41));

    m_track(lons,lats,dates,'ticks',0,'times',4,'dates',8,...
           'clip','off','color','r','orient','upright');  
           
  case 11,
  
    % example showing range rings
    
    clf
    m_proj('hammer','clong',170);
    m_grid('xtick',[],'ytick',[],'linestyle','-');
    m_coast('patch','g');
    m_line(100.5,13.5,'marker','s','color','r');
    m_range_ring(100.5,13.5,[1000:1000:15000],'color','b','linewidth',2);
    xlabel('1000km range rings from Bangkok');
    
  case 12,
  
    % Example showing speckle
    
    bndry_lon=[-128.8 -128.8 -128.3 -128 -126.8 -126.6 -128.8];
    bndry_lat=[49      50.33  50.33  50   49.5   49     49];

    clf;
    m_proj('lambert','long',[-130 -121.5],'lat',[47 51]);
    m_gshhs_i('color','k');
    m_gshhs_i('speckle','color','k');
    m_line(bndry_lon,bndry_lat,'linewidth',2,'color','k');     % Area outline ...
    m_hatch(bndry_lon,bndry_lat,'single',30,5,'color','k'); % ...with hatching added.

    m_grid('linewidth',2,'linestyle','none');
    title('Speckled Boundaries for nice B&W presentation (best in postscript format)','fontsize',14);
    m_text(-128,48,5,{'Pacific','Ocean'},'fontsize',18);
        
  case 13,
  
    % Colouring the ocean blue
    
    clf
    m_proj('miller','lat',[-75 75]);
    set(gca,'color',[.9 .99 1]);
    m_coast('patch',[.7 1 .7],'edgecolor','none');
    m_grid('box','fancy','linestyle','none');
       
    cities={'Cairo','Washington','Buenos Aires'};
    lons=[ 30+2/60  -77-2/60   -58-22/60];
    lats=[ 31+21/60  38+53/60  -34-45/60];
    
    for k=1:3,
      [range,ln,lt]=m_lldist([-123-6/60 lons(k)],[49+13/60  lats(k)],40);
      m_line(ln,lt,'color','r','linewidth',2);
      m_text(ln(end),lt(end),sprintf('%s - %d km',cities{k},round(range)));
    end;
      
    title('Great Circle Routes','fontsize',14,'fontweight','bold');
       
  end;
  
 if i<length(index), 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   disp('  hit return to continue');
   pause
   disp('        ...drawing');
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 end;
end;

  
