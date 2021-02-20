function m_utmgrid(varargin)
% M_UTMGRID Draws a UTM grid on a map.
%        M_UTMGRID('parameter','value',...) with any number (or no)
%        optional parameters is used to draw a UTM grid (i.e. with
%        grid lines at 1 km intervals) for a map using the UTM
%        coordinate system. This can be used instead of, or in addition
%        to, the latitude/longitude grid drawn by M_GRID.
%
%        The optional parameters allow the user
%        to control the look of the grid. These parameters are listed
%        by M_UTMGRID('get'), with default parameters in M_UTMGRID('set');
%
%        This is probably useful only for maps less than a few km across.
%
%        see also M_PROJ, M_GRID

% Rich Pawlowicz (rich@eoas.ubc.ca) 25/April/2018

% These structures are initialized by m_proj()

global MAP_PROJECTION MAP_VAR_LIST

% Have to have initialized a map first

if isempty(MAP_PROJECTION)
  disp('No Map Projection initialized - call M_PROJ(''UTM'',... first!');
  return;
elseif ~strcmp(MAP_PROJECTION.name,'UTM')  
  disp('Not a UTM projection');
  return;
elseif strcmp(MAP_VAR_LIST.ellipsoid,'normal')
  disp('Cannot show UTM coordinates using the ``normal'' (earth radius=1) ellipsoid in the m_proj call');
  return;
end
  
gcolor=get(gca,'gridcolor');
gxcolor=get(gca,'xcolor');
gycolor=get(gca,'ycolor');
glinestyle=get(gca,'gridlinestyle');
glinewidth=get(gca,'linewidth');
gfontsize=get(gca,'fontsize');
gfontname=get(gca,'fontname');
gxaxisloc=get(gca,'xaxislocation'); 
gyaxisloc=get(gca,'yaxislocation');
gtickdir=get(gca,'tickdir'); 
gticklen=get(gca,'ticklength'); gticklen=gticklen(1); 
 
k=1;
while k<=length(varargin)
  switch lower(varargin{k}(1:3))
    case 'gri'
      gcolor=varargin{k+1};      
    case 'lin'
      switch lower(varargin{k}(1:5))
         case 'linew'
           glinewidth=varargin{k+1};
         case 'lines'
           glinestyle=varargin{k+1};
      end
    case 'fon'
       switch lower(varargin{k}(1:5))
         case 'fonts'
           gfontsize=varargin{k+1};
         case 'fontn'
           gfontname=varargin{k+1};
       end
    case 'xco'
      gxcolor=varargin{k+1};
    case 'yco'
      gycolor=varargin{k+1};
    case 'xax'
      gxaxisloc=varargin{k+1};
    case 'yax'
      gyaxisloc=varargin{k+1};
    case 'tic'
      switch lower(varargin{k}(1:5))
        case 'tickl'
           gticklen=varargin{k+1};
        case 'tickd'
           gtickdir=varargin{k+1};
      end
    case {'get','usa'}

      disp('      ''ticklength'',value');
      disp('      ''tickdir'',( ''in'' | ''out'' )');
      disp('      ''gridcolor'',colorspec');
      disp('      ''linewidth'', value');
      disp('      ''linestyle'', ( linespec | ''none'' )');
      disp('      ''fontsize'',value');
      disp('      ''fontname'',name');
      disp('      ''Xcolor'',colorspec');
      disp('      ''Ycolor'',colorspec');
      disp('      ''XaxisLocation'',( ''bottom'' | ''top'' ) ');
      disp('      ''YaxisLocation'',( ''left'' | ''right'' ) ');
      return;
    case 'set'
      disp(['      ticklength = ' num2str(gticklen)]);
      disp(['      tickdir = ' gtickdir]);
      disp(['      gridcolor = ' gcolor]);
      disp(['      linewidth = ' num2str(glinewidth)]);
      disp(['      linestyle = ' glinestyle]);
      disp(['      fontsize = ' num2str(gfontsize)]);
      disp(['      fontname = ' gfontname]);
      disp(['      Xcolor = ' gxcolor]);
      disp(['      Ycolor = ' gycolor]);
      disp(['      XaxisLocation = ' gxaxisloc]);
      disp(['      YaxisLocation = ' gyaxisloc]);
      return;
  end
  k=k+2;
end   





% Bring the grid back!

set(gca,'visible','on','layer','top','tickdir',gtickdir,...
   'gridcolor',gcolor,'gridlinestyle',glinestyle,'linewidth',glinewidth,...
   'xaxislocation',gxaxisloc,'yaxislocation',gyaxisloc,'xcolor',gxcolor',...
   'ycolor',gycolor,'fontname',gfontname,'ticklength',[gticklen .025],...
   'fontsize',gfontsize,'xgrid','on','ygrid','on');
 

fnt=gfontsize;
fnt2=fnt*0.8;

% Do the x axis

tk=get(gca,'xtick');
% Nothing smaller than 1 km squares
tk(rem(tk,1e3)~=0)=[];

tk1=fix(tk/1e5);
tk2=fix((tk-1e5*tk1)/1e3);
tk3=rem(tk,1e3);
 
for k=1:length(tk)-1,
    lab{k}=sprintf('\\fontsize{%d}%02d',fnt,tk2(k));
end
k=length(tk);
lab{k}=sprintf('{\\fontsize{%d}%2d}\\fontsize{%d}%02d{\\fontsize{%d}%03d} E',fnt2,tk1(k),fnt,tk2(k),fnt2,tk3(k));
set(gca,'xtick',tk,'xticklabel',lab);


% Do the y axis

tk=get(gca,'ytick');
% Nothing smaller than 1 km squares
tk(rem(tk,1e3)~=0)=[];

tk1=fix(tk/1e5);
tk2=fix((tk-1e5*tk1)/1e3);
tk3=rem(tk,1e3);
 
for k=1:length(tk)-1,
    lab{k}=sprintf('\\fontsize{%d}%02d',fnt,tk2(k));
end
k=length(tk);
lab{k}=sprintf('{\\fontsize{%d}%2d}\\fontsize{%d}%02d{\\fontsize{%d}%03d} N',fnt2,tk1(k),fnt,tk2(k),fnt2,tk3(k));
set(gca,'ytick',tk,'yticklabel',lab','yticklabelrotation',90);

set(gca,'tag','m_utmgrid');
