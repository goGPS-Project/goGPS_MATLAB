function h=m_ruler(posx,posy,varargin);
% M_RULER Draws a distance scalebar for a map
%   M_RULER([X1 X2],Y1) draws a horizontal scale bar between the
%   normalized coordinates (X1,Y1) and (X2,Y1) where both X/Y are
%   in the range 0 to 1.  
%   M_RULER(X1,[Y1 Y2]) draws a vertical scalebar
%   M_RULER(...,NINTS) draws the scalebar with NINTS intervals
%   if NINTS is a scalar (default 4).  Distances of each interval are 
%   chosen to be 'nice'. If NINTS is a vector it is understood to be
%   the distances to be used (in meters)
%
%   M_SCALEBAR(....,'parameter','value',...) lets you specify
%   extra parameter/value pairs in the usual handle-graphics way.
%   'color' and 'fontsize' are probably the most useful, 'tickdir'
%   'in' and 'out' chooses between different styles.
%
%   Probably BEST to call this AFTER M_GRID otherwise placement might
%   seem a bit odd.
%
% WARNING - the scalebar is probably not useful for any global
%           (i.e. whole-world) or even a significant-part-of-the-globe
%           map, but I won't stop you using it. Caveat user!

% R. Pawlowicz rich@eos.ubc.ca  8/Nov/2006
% 7/Dec/11 - Octave 3.2.3 compatibility

% Recognize Octave
a=ver;
if strcmp(a(1).Name,'Octave'),
 IsOctave=logical(1);
else
 IsOctave=logical(0);
end;


nints=4;
fixticks=0;

if length(varargin)>0,
 if isnumeric(varargin{1}),
   nints=varargin{1};
   varargin=varargin(2:end);
   fixticks=1;
  end;
end;    

 
gcolor='k';
glinestyle='-';
glinewidth=3;
gfontsize=get(gca,'fontsize');
gfontname=get(gca,'fontname');
gticklen=get(gca,'ticklength'); gticklen=gticklen(1); 
gtickdir=get(gca,'tickdir'); 
 
 
% Parse parameter list for options. I really should do some
% error checking here, but...

k=1;
while k<=length(varargin),
  switch lower(varargin{k}(1:3)),
    case 'col',
      gcolor=varargin{k+1};
    case 'lin',
      switch lower(varargin{k}(1:5)),
         case 'linew',
           glinewidth=varargin{k+1};
         case 'lines',
           glinestyle=varargin{k+1};
      end;
    case 'fon',
       switch lower(varargin{k}(1:5)),
         case 'fonts',
           gfontsize=varargin{k+1};
         case 'fontn',
           gfontname=varargin{k+1};
       end;
    case 'tic',
       switch lower(varargin{k}(1:5)),
         case 'tickl',
           gticklen=varargin{k+1};
         case 'tickd',
           gtickdir=varargin{k+1};
       end;
    case {'get','usa'},
      disp('      ''ticklength'',value');
      disp('      ''tickdir'',( ''in'' | ''out'' )');
      disp('      ''color'',colorspec');
      disp('      ''linewidth'', value');
      disp('      ''linestyle'', ( linespec | ''none'' )');
      disp('      ''fontsize'',value');
      disp('      ''fontname'',name');
       return;
    case 'set',
      disp(['      ticklength = ' num2str(gticklen)]);
      disp(['      tickdir = ' gtickdir]);
       disp(['      color = ' gcolor]);
      disp(['      linewidth = ' num2str(glinewidth)]);
      disp(['      linestyle = ' glinestyle]);
      disp(['      fontsize = ' num2str(gfontsize)]);
      disp(['      fontname = ' gfontname]);
       return;
  end;
  k=k+2;
end;     

 
% Need earth radius, in m.
erad=6378137; %m (from WGS-84)


if  ( length(posx)==2 & length(posy)==1),
  posy=[posy posy];
  horiz=1;
elseif  ( length(posx)==1 & length(posy)==2),
  posx=[posx posx];
  horiz=0;
end;




xlm=get(gca,'xlim');
ylm=get(gca,'ylim');

% Get into screen coords
posx=xlm(1) + posx*diff(xlm);
posy=ylm(1) + posy*diff(ylm);


if diff(xlm)>10, % we are probably already in meters, i.e. UTM
   scfac=1;
else
  scfac=erad;
end;



distance=(diff(posx)+diff(posy))*scfac;


niceints=[ 1      2      5  ...
           10     20     50 ...
	   100    200    500  ...
	   1000   2000   5000  ...
	   10000  20000  50000  ...
	   100000 200000 500000 ...
	   1000000 2000000 5000000 ];
	   
 	   
if length(nints)==1,
	  
  exactint=distance/(nints);
  [dun,I]=min(abs(niceints-exactint));
 
  if ~fixticks,
    nints=fix(distance/niceints(I));
  end;
  dist=[0:nints]*niceints(I);

else
  dist=(nints-nints(1));
  nints=length(dist)-1;
end;

if max(log10(dist(2:end)))>=3,
   numfac=1000;
else
   numfac=1;
end;
      
if horiz,
 
 
  if strcmp(gtickdir,'in'),

    line(posx(1)+[0 dist(end)/scfac],posy(1)+[0 0],'color',gcolor,'linewidth',glinewidth,'linestyle',glinestyle,...
	 'clipping','off','tag','m_ruler_x');
	 
    XX=posx(1)+[dist;dist]/scfac;
    YY=posy(1)+diff(ylm)*gticklen*[-1;1]*ones(1,nints+1);
    if IsOctave,
      for k=1:size(XX,2),
        line(XX(:,k),YY(:,k),...
	 'color',gcolor,'linewidth',glinewidth/3,'linestyle',glinestyle,...
	 'clipping','off','tag','m_ruler_y');
      end;
    else  
      line(XX,YY,...
	 'color',gcolor,'linewidth',glinewidth/3,'linestyle',glinestyle,...
	 'clipping','off','tag','m_ruler_y');
    end;	 
  else
   
    patch(posx(1)+[dist(1:end-1);dist(1:end-1);dist(2:end);dist(2:end)]/scfac,...
           posy(1)+diff(ylm)*gticklen*[-1;1;1;-1]*ones(1,nints),...
	   repmat(bitmax  ,4,nints),...
           reshape(rem(rem(0:nints*3-1,nints),2)==0,1,nints,3),...
	 'linewidth',glinewidth/3,'linestyle',glinestyle,...
         'clipping','off','tag','m_ruler');

  end;
   
  if nints>1,
    for k=1:nints,
      text(posx(1)+dist(k)/scfac,posy(1)-diff(ylm)*gticklen*2,sprintf('%d',dist(k)/numfac), ...
          'fontsize',gfontsize,'fontname',gfontname,...
         'verticalalignment','top','horizontalalignment','center','tag','m_ruler_label');
    end;
    text(posx(1)+dist(nints+1)/scfac,posy(1)-diff(ylm)*gticklen*2,sprintf('%d km',dist(end)/numfac),...
          'fontsize',gfontsize,'fontname',gfontname,...
          'verticalalignment','top','horizontalalignment','center','tag','m_ruler_label');
  else 
    text(posx(1)+mean(dist)/scfac,posy(1)-diff(ylm)*gticklen*2,sprintf('%d km',dist(end)/numfac),...
          'fontsize',gfontsize,'fontname',gfontname,...
          'verticalalignment','top','horizontalalignment','center','tag','m_ruler_label');
  end;
  
else,

  if strcmp(gtickdir,'in'),

    line(posx(1)+[0 0],posy(1)+[0 dist(end)/scfac],'color',gcolor,'linewidth',glinewidth,'linestyle',glinestyle,...
	 'clipping','off','tag','m_scalebar_x');
	 
    XX=posx(1)+diff(xlm)*gticklen*[-1;1]*ones(1,nints+1);
    YY=posy(1)+[dist;dist]/scfac;
    
    if IsOctave,
      for k=1:size(XX,2),
       line(XX(:,k),YY(:,k),...
	 'color',gcolor,'linewidth',glinewidth/3,'linestyle',glinestyle,...
	 'clipping','off','tag','m_ruler_y');
      end;
    else  
       line(XX,YY,...
	 'color',gcolor,'linewidth',glinewidth/3,'linestyle',glinestyle,...
	 'clipping','off','tag','m_ruler_y');
    end;	 
  else
   
    patch(posx(1)+diff(xlm)*gticklen*[-1;1;1;-1]*ones(1,nints),...
          posy(1)+[dist(1:end-1);dist(1:end-1);dist(2:end);dist(2:end)]/scfac,...
	   repmat(bitmax  ,4,nints),...
           reshape(rem(rem(0:nints*3-1,nints),2)==0,1,nints,3),...
	 'linewidth',glinewidth/3,'linestyle',glinestyle,...
         'clipping','off','tag','m_ruler');

  end;

  if nints>1,
    for k=1:nints,
      text(posx(1)+diff(xlm)*gticklen*2,posy(1)+dist(k)/scfac,sprintf('%d',dist(k)/numfac), ...
          'fontsize',gfontsize,'fontname',gfontname,...
         'verticalalignment','middle','horizontalalignment','left','tag','m_ruler_label');
    end;
    text(posx(1)+diff(xlm)*gticklen*2,posy(1)+dist(nints+1)/scfac,sprintf('%d km',dist(end)/numfac),...
          'fontsize',gfontsize,'fontname',gfontname,...
          'verticalalignment','middle','horizontalalignment','left','tag','m_ruler_label');
  else  
    text(posx(1)+diff(xlm)*gticklen*2,posy(1)+mean(dist)/scfac,sprintf('%d km',dist(end)/numfac),...
          'fontsize',gfontsize,'fontname',gfontname,...
          'verticalalignment','middle','horizontalalignment','left','tag','m_ruler_label');
  end;
       
end;







 







