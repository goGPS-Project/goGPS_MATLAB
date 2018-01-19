function m_grid(varargin)
% M_GRID make a grid on a map.
%        M_GRID('parameter','value',...) with any number (or no)
%        optional parameters is used to draw a lat/long grid for a
%        previously initialized map projection.
%
%        The optional parameters allow the user
%        to control the look of the grid. These parameters are listed
%        by M_GRID('get'), with default parameters in M_GRID('set');
%
%        see also M_PROJ

% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.
%
%  19/6/97 - set visibility of titles and so forth to 'on' (they
%            default to 'off' when axes visibility is turned off)
%  2/11/97 - for M5.1, the old way of making the patch at the bottom (i.e.
%            by rearranging the axes children) instead causes matlab to lose
%            track of titles. Try a different fix.
% 11/01/98 - Added way of making longitude lines cut off to prevent crowding near poles (you have
%            to specify a vector for allowabale latitudes for this to work).
% 16/02/98 - Made a little fudge to allow the user to fully specify grid location
%            without getting the edge points. It doesn't quite work if only *one* edge
%            point is desired....but I hope it will be OK that way.
% 19/02/98 - PC-users complain about layers getting out of order! Their fault for using
%            such an awful OS...however (with help from Eric Firing) I think I have
%            a fix.
%  7/04/98 - Another fix to grid locations to not automatically add edge points
%            (as requested by EF)
%  7/05/98 - Added 'fancy' outline box.
% 14/11/98 - Changed tag names from m_* to m_grid_*.
% 11/07/99 - Apparently fontname changing didn't work (thanks to Dave McCollum)
% 28/04/04 - Changed m_grid_color code so it works right under unix; old
%            way retained for windows (ugh).
% 16/10/05 - Kirk Ireson discovered that the way to fix those annoying 'cut-throughs'
%            in fancy_box was to add a 'large' zdata...so I've adapted his fix in
%            fancybox and fancybox2.
% 21/11/06 - added 'backcolor'
% 16/4/07  - sorted ticklabels when user-specified (prevents an odd problem near in
%            azimuthal projections).
% 4/DEc/11 - isstr to ischar
% 7/Dec/11 - Octave 3.2.3 compatibility
% 8/Sep/13 - added 'tickstyle' parameter
% 27/Sep/13 - matlab 2013b out, includes graphic bug. Workaround provided by
%             Corinne Bassin.
% 10/Jul/14 - in 2014a BITMAX starts not to be used, changed to FLINTMAX...
% 13/Nov/14 - 2014b graphics changes are biting; a number of version-dependent
%              fixes implemented.
% 17/Nov/17 - All kinds of "new graphics" changes in the background patch
% 27/Nov/17 - made a specialized fix for xaxis tick directions in Robinson and
%             Sinusoidal projections. (problem noted by J. Tsoa, but I didn't like
%             the way he fixed it. Also improved consistency of "lake"
%             handling with non-white background colours, and more Octave
%             compatibility.


% Note that much of the work in generating line data 
% is done by calls to the individual projections - 
% most of M_GRID is concerned with the mechanics of plotting


% These structures are initialized by m_proj()

global MAP_PROJECTION MAP_VAR_LIST

% Have to have initialized a map first

if isempty(MAP_PROJECTION)
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end
  


% Otherwise we are drawing a grid!

% Set current projection to geographic
Currentmap=m_coord('set');
m_coord(MAP_PROJECTION.coordsystem.name);

% Default parameters for grid

xtick=6;
ytick=6;
xlabels=NaN;
ylabels=NaN;
gcolor='k';
ggridcolor=[];
gbackcolor='w'; %%get(gcf,'color');
glinestyle=':';
glinewidth=get(gca,'linewidth');
gbox='on'; 
gfontsize=get(gca,'fontsize');
gfontname=get(gca,'fontname');
gxaxisloc=get(gca,'xaxislocation'); 
gyaxisloc=get(gca,'yaxislocation');
gtickdir=get(gca,'tickdir'); 
gticklen=get(gca,'ticklength'); gticklen=gticklen(1); 
gxticklabeldir='middle';
gyticklabeldir='end';
gtickstyle='dm';

dpatch=5; % interpolation factor for fancy grids

% Parse parameter list for options. I really should do some
% error checking here, but...

k=1;
while k<=length(varargin)
  switch lower(varargin{k}(1:3))
    case 'box'
      gbox=varargin{k+1};
    case 'xti'
      if length(varargin{k})==5
        xtick=sort(varargin{k+1});   % Added 'sort' here for people who put things in
      else                           % a random order near poles
        xlabels=varargin{k+1};
      end
    case 'yti'
      if length(varargin{k})==5
        ytick=sort(varargin{k+1});
      else
        ylabels=varargin{k+1};
      end
    case 'xla'
      gxticklabeldir=varargin{k+1};
    case 'yla'
      gyticklabeldir=varargin{k+1};
    case 'col'
      gcolor=varargin{k+1};
    case 'gri'
      ggridcolor=varargin{k+1};      
    case 'bac'
      gbackcolor=varargin{k+1};
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
	    case 'ticks'
	       gtickstyle=varargin{k+1};   
      end
    case {'get','usa'}
      disp('      ''box'',( ''on'' | ''fancy'' | ''off'' )');
      disp('      ''xtick'',( num | [value1 value2 ...])');
      disp('      ''ytick'',( num | [value1 value2 ...])');
      disp('      ''xticklabels'',[label1;label2 ...]');
      disp('      ''yticklabels'',[label1;label2 ...]');
      disp('      ''xlabeldir'', ( ''middle'' | ''end'' )');
      disp('      ''ylabeldir'', ( ''end'' | ''middle'' )');
      disp('      ''ticklength'',value');
      disp('      ''tickdir'',( ''in'' | ''out'' )');
      disp('      ''tickstyle'',(''dm'' | ''dd'' )');  % deg-min or decimal-deg
      disp('      ''color'',colorspec');
      disp('      ''gridcolor'',colorspec');
      disp('      ''backgroundcolor'',colorspec');
      disp('      ''linewidth'', value');
      disp('      ''linestyle'', ( linespec | ''none'' )');
      disp('      ''fontsize'',value');
      disp('      ''fontname'',name');
      disp('      ''XaxisLocation'',( ''bottom'' | ''middle'' | ''top'' ) ');
      disp('      ''YaxisLocation'',( ''left'' | ''middle'' | ''right'' ) ');
      return;
    case 'set'
      disp(['      box = ' gbox]);
      disp(['      xtick = ' num2str(xtick)]);
      disp(['      ytick = ' num2str(ytick)]);
      disp(['      ticklength = ' num2str(gticklen)]);
      disp(['      tickdir = ' gtickdir]);
      disp(['      tickstyle = ' gtickstyle]);
      disp(['      xlabeldir = ' gxticklabeldir]);
      disp(['      ylabeldir = ' gyticklabeldir]);
      disp(['      color = ' gcolor]);
      disp(['      gridcolor = ' gcolor]);
      disp(['      backgroundcolor - ' gbackcolor]);
      disp(['      linewidth = ' num2str(glinewidth)]);
      disp(['      linestyle = ' glinestyle]);
      disp(['      fontsize = ' num2str(gfontsize)]);
      disp(['      fontname = ' gfontname]);
      disp(['      XaxisLocation = ' gxaxisloc]);
      disp(['      YaxisLocation = ' gyaxisloc]);
      return;
  end
  k=k+2;
end   

if isempty(ggridcolor)
    ggridcolor=gcolor;
end

  
if strcmp(gbox,'fancy')
  if strcmp(MAP_VAR_LIST.rectbox,'on') || strcmp(MAP_VAR_LIST.rectbox,'circle')
   gbox='on';
   warning([' No fancy outline with ''rectbox'' set to ''' MAP_VAR_LIST.rectbox '''']);
  end
end

[X,Y]=feval(MAP_PROJECTION.routine,'box');

% These color parts make Octave fail so I moved them here, otherwise 
% they should be before this if block
    
% If there are lakes around, make them the backgroundcolor..
hh=get(gca,'children');
hh_tags=get(hh,'tag');
if length(hh)==1, hh_tags={hh_tags}; end
for i=1:length(hh_tags)
   if ~isempty(strfind(hh_tags{i},'_lake')) && strcmp(get(hh(i),'type'),'patch')
     set(hh(i),'facecolor',gbackcolor);     
   end
end    

% Set the axes colour so lakes get coloured right (if a coastline added
% next)

set(gca,'color',gbackcolor');
     

if MAP_PROJECTION.IsOctave


  %  Ocatve handling here all provided by Jamie Tsao

  % Note that both gnuplot and fltk doesn't really work well with the z-axis,
  % so here we simply ignore it.

  % The tag is important for M_COAST to catch the correct underlying color.
  % As a result, though, M_COAST should be called after M_GRID.
  
  if strcmp(graphics_toolkit(), 'gnuplot')
    % Unfortunately, gnuplot also has a bug with using patch, where a linestyle
    % of 'none' would cause it to barf. This bug has been reported since 2013,
    % but apparently it is still not fixed. It may have something to do with
    % how Octave interacts with the underlying gnuplot backend.

    % On the other hand, this means I can set the linewidth to 0 and it should
    % work exactly as intended. If there is no problems, then I can refactor
    % this part, assuming ftlk doesn't barf on a zero linewidth.    
    patch('xdata', X(:), 'ydata', Y(:), 'facecolor', gbackcolor, ...
          'edgecolor', 'k', 'linewidth', 0, 'tag', 'm_grid_color');
  else
     
    patch('xdata', X(:), 'ydata', Y(:), 'facecolor', gbackcolor, ...
          'edgecolor', 'k', 'linestyle', 'none', 'tag', 'm_grid_color');
  end

  % Now I set it at the bottom of the children list so it gets drawn first
  % (i.e. doesn't cover anything)

  % A bug in Octave 3.8.1 makes set(gca, 'children', -) insert the new
  % permutation before the old one instead of replacing it. This only occurs
  % when the ShowHiddenHandles flag is set to on instead of off. Although fltk
  % is not affected, gnuplot still is. Combined with fltk's lack of support for
  % unicode and hence the degree sign, this is still important. Furthermore,
  % this bug causes delete() to function incorrectly.
 
  show = get(0, 'ShowHiddenHandles');
  set(0, 'ShowHiddenHandles', 'on');

  hh = get(gca,'children');
  htags = get(hh,'tag');
  k = strmatch('m_grid_color', htags);
  hht = hh;
  hh(k) = [];
  hh = [hh;hht(k)];

  % Continuing on the bug, there is no particular need for why this needs to be
  % placed under even hidden handles. Even so, as a workaround, we can grab all
  % hidden handles, set them to be visible, turn off flag, reorder, then set
  % them back as hidden, as noted below.

  hht_inv = strmatch('off', get(hht, 'HandleVisibility'));
  set(hht(hht_inv), 'HandleVisibility', 'on');
  set(0, 'ShowHiddenHandles', 'off'); 
  set(gca,'children',hh);
  set(hht(hht_inv), 'HandleVisibility', 'off');

  set(0, 'ShowHiddenHandles', show);
 

else   % For MATLAB    
    
  
     
   % This is a very problematic part of the code. It turns out the the interaction between
   % PATCH objects and CONTOURF objects does not work correctly in the Painters renderer -
   % this is true in all versions up to 7.7 at least. Patches with large negative Z just
   % don't get drawn under contourgroup patches.
   % 
   % There are several possible workarounds:
   %
   %  1) Make sure you use the 'v6' option in contourf calls (see m_contourf.m to see
   %     how I have tried to do that for some versions of matlab)
   %      - problem: the 'v6' option is going away soon, also you may want the
   %                contourgroup object that the v6 option destroys.
   %
   %  2) Change the renderer to something else:
   %        set(gcf,'renderer','opengl') or
   %        set(gcf,'renderer','zbuffer')
   %      - problem: These other renderers are not available on all systems, they may also
   %                give less-precise results.
   %
   %  3) Use the painters renderer, but reorder the children so that the patch is at the
   %     bottom of the list (painters cares about child order)
   %      - problem: sometimes the child order is rearranged if you click on the figure,
   %                 also (at least in some versions of matlab) this causes labels to
   %                 disappear.
   %
   % With version 7.4 onwards I have discovered that reordering the children apparently
   % is Mathworks-blessed (c.f. the UISTACK function).  So I am going to try to implement
   % the latter as a default.
   %
   % Now, putting in a white background works under linux (at least) and
   % NOT under windows...I don't know about macs.
   %
   % This broke with matlab 2014b - but eventually (I sent in a "bug
   % report" and found it was a "feature") I learned of the 'sortmethod'
   % tag, which made it work again - as long as I didn't use negative z
   % values in the patch command. I don't know WHY that is...
    
   if ~MAP_PROJECTION.newgraphics
        patch('xdata',X(:),'ydata',Y(:),'zdata',-MAP_PROJECTION.LARGVAL*ones(size(X(:))),'facecolor',gbackcolor,...
	   'edgecolor','k','linestyle','none','tag','m_grid_color');
	  % Now I set it at the bottom of the children list so it gets drawn first (i.e. doesn't
	  % cover anything)
	 show=get(0, 'ShowHiddenHandles');
	 set(0, 'ShowHiddenHandles', 'on');
	 hh=get(gca,'children');
	 htags = get(hh,'tag');
	 k = strmatch('m_grid_color',htags);
	 hht = hh;
	 hh(k) = [];
	 hh = [hh;hht(k)];
	 set(gca,'children',hh);
	 set(0, 'ShowHiddenHandles', show);
     
   else  % after 2014b
       
     patch('xdata',X(:),'ydata',Y(:),'zdata',-MAP_PROJECTION.LARGVAL*ones(size(X(:))),'facecolor',gbackcolor,...
             'edgecolor','k','linestyle','none','tag','m_grid_color');
%     patch('xdata',X(:),'ydata',Y(:),'facecolor',gbackcolor,...
%             'edgecolor','k','linestyle','none','tag','m_grid_color');
	 set(gca,'sortmethod','childorder');
     uistack(findobj(gca,'tag','m_grid_color'),'bottom');

   end

   
end


% X-axis labels and grid

if ~isempty(xtick)

 % Tricky thing - if we are drawing a map with the poles, its nasty when the lines get too close
 % together. So we can sort of fudge this by altering MAP_VAR_LIST.lats to be slightly smaller,
 % and then changing it back again later.
 fudge_north='n';fudge_south='n';
 if ~isempty(ytick) && length(ytick)>1
  if MAP_VAR_LIST.lats(2)==90 
    fudge_north='y';
    MAP_VAR_LIST.lats(2)=ytick(end);
  end
  if MAP_VAR_LIST.lats(1)==-90 
    fudge_south='y';
    MAP_VAR_LIST.lats(1)=ytick(1);
  end
 end

 [X,Y,lg,lgI]=feval(MAP_PROJECTION.routine,'xgrid',xtick,gxaxisloc,gtickstyle);
 [labs,scl]=m_labels('lon',lg,xlabels,gtickstyle);
 
 % Draw the grid. Every time we draw something, I first reshape the matrices into a long
 % row so that a) it plots faster, and b) all lines are given the same handle (which cuts
 % down on the number of children hanging onto the axes).

 [n,m]=size(X);
 line(reshape([X;NaN+ones(1,m)],(n+1)*m,1),reshape([Y;NaN+ones(1,m)],(n+1)*m,1),...
      'linestyle',glinestyle,'color',ggridcolor,'linewidth',0.1,'tag','m_grid_xgrid');

 % Get the tick data
 [ltx,lty,utx,uty]=maketicks(X,Y,gticklen,gtickdir);

 % Draw ticks if labels are on top or bottom (not if they are in the middle)

 if strcmp(gxticklabeldir,'middle')
  if lgI==size(X,1) && strcmp(gxaxisloc,'top')  % Check to see if the projection supports this option.
   vert='bottom';horiz='center';drawticks=1;
   xx=utx(1,:);yy=uty(1,:);rotang=atan2(diff(uty),diff(utx))*180/pi+90;
  elseif lgI==1 && strcmp(gxaxisloc,'bottom')
   vert='top';horiz='center';drawticks=1;
   xx=ltx(1,:);yy=lty(1,:);rotang=atan2(diff(lty),diff(ltx))*180/pi-90;
  else
   vert='middle';horiz='center';lgIp1=lgI+1;drawticks=0;
   xx=X(lgI,:); yy=Y(lgI,:);rotang=atan2(Y(lgIp1,:)-Y(lgI,:),X(lgIp1,:)-X(lgI,:))*180/pi-90;
  end
 else
  if lgI==size(X,1) && strcmp(gxaxisloc,'top')  % Check to see if the projection supports this option.
   vert='middle';horiz='left';drawticks=1;
   xx=utx(1,:);yy=uty(1,:);rotang=atan2(diff(uty),diff(utx))*180/pi+180;
  elseif lgI==1 && strcmp(gxaxisloc,'bottom')
   vert='middle';horiz='right';drawticks=1;
   xx=ltx(1,:);yy=lty(1,:);rotang=atan2(diff(lty),diff(ltx))*180/pi;
  else
   vert='top';horiz='center';lgIp1=lgI+1;drawticks=0;
   xx=X(lgI,:); yy=Y(lgI,:);rotang=atan2(Y(lgIp1,:)-Y(lgI,:),X(lgIp1,:)-X(lgI,:))*180/pi;
  end
 end

 % For some projections it is better if the x-grid labels are horizontal or
 % vertical
 
 if any(strcmp(MAP_PROJECTION.name,{'Robinson','Sinusoidal'}) )
     if strcmp(gxaxisloc,'bottom') || strcmp(gxaxisloc,'top')
         if strcmp(gxticklabeldir,'middle')
            rotang = zeros(size(rotang));    % flat
            yy=repmat(mean(yy),size(yy));    % Put them all on the same y-value.
         else
            rotang = zeros(size(rotang))+90;  % upright
         end 
     end    
 end    
 
 
 if strcmp(gbox,'fancy')
    if gtickdir(1)=='i'
      fancybox(lg,MAP_VAR_LIST.longs,'xgrid','bottom',dpatch,gticklen,gtickstyle); 
      drawticks=0;
    else    
      fancybox2(lg,MAP_VAR_LIST.longs,'xgrid','bottom',dpatch,gticklen,gtickstyle); 
    end
 end    
 if drawticks
   [n,m]=size(ltx);
   line(reshape([ltx;NaN+ones(1,m)],(n+1)*m,1),reshape([lty;NaN+ones(1,m)],(n+1)*m,1),...
        'linestyle','-','color',gcolor,'linewidth',glinewidth,'tag','m_grid_xticks-lower','clipping','off');
   line(reshape([utx;NaN+ones(1,m)],(n+1)*m,1),reshape([uty;NaN+ones(1,m)],(n+1)*m,1),...
        'linestyle','-','color',gcolor,'linewidth',glinewidth,'tag','m_grid_xticks-upper','clipping','off');
 end

 % Add the labels! (whew)

 ik=1:size(X,2);

 for k=ik
   [rotang(k), horizk, vertk] = upright(rotang(k), horiz, vert);
   text(xx(k),yy(k),labs{k},'horizontalalignment',horizk,'verticalalignment',vertk, ...
        'rotation',rotang(k),'fontsize',gfontsize*scl(k),'color',gcolor,...
        'tag','m_grid_xticklabel','fontname',gfontname);
 end

 if fudge_north=='y'
   MAP_VAR_LIST.lats(2)=90;
 end
 if fudge_south=='y'
   MAP_VAR_LIST.lats(1)=-90;
 end

end

if ~isempty(ytick)
 % Y-axis labels and grid

 [X,Y,lt,ltI]=feval(MAP_PROJECTION.routine,'ygrid',ytick,gyaxisloc,gtickstyle);
 [labs,scl]=m_labels('lat',lt,ylabels,gtickstyle);

 % Draw the grid
 [n,m]=size(X);
 line(reshape([X;NaN+ones(1,m)],(n+1)*m,1),reshape([Y;NaN+ones(1,m)],(n+1)*m,1),...
      'linestyle',glinestyle,'color',ggridcolor,'linewidth',0.1,'tag','m_grid_ygrid');

 % Get the tick data
 [ltx,lty,utx,uty]=maketicks(X,Y,gticklen,gtickdir);

 % Draw ticks if labels are on left or right (not if they are in the middle)
 if strcmp(gyticklabeldir,'end')
  if ltI==size(X,1) && strcmp(gyaxisloc,'right')  % Check to see if the projection supports this option.
   horiz='left';vert='middle';drawticks=1;
   xx=utx(1,:);yy=uty(1,:);rotang=atan2(diff(uty),diff(utx))*180/pi+180;
  elseif ltI==1 && strcmp(gyaxisloc,'left')
   horiz='right';vert='middle';drawticks=1;
   xx=ltx(1,:);yy=lty(1,:);rotang=atan2(diff(lty),diff(ltx))*180/pi;
  else
   horiz='center';vert='top';ltIp1=ltI+1;drawticks=0;
   xx=X(ltI,:); yy=Y(ltI,:);rotang=atan2(Y(ltIp1,:)-Y(ltI,:),X(ltIp1,:)-X(ltI,:))*180/pi;
  end
 else
  if ltI==size(X,1) && strcmp(gyaxisloc,'right')  % Check to see if the projection supports this option.
   horiz='center';vert='top';drawticks=1;
   xx=utx(1,:);yy=uty(1,:);rotang=atan2(diff(uty),diff(utx))*180/pi+270;
  elseif ltI==1 && strcmp(gyaxisloc,'left')
   horiz='center';vert='bottom';drawticks=1;
   xx=ltx(1,:);yy=lty(1,:);rotang=atan2(diff(lty),diff(ltx))*180/pi+90;
  else
   horiz='left';vert='middle';ltIp1=ltI+1;drawticks=0;
   xx=X(ltI,:); yy=Y(ltI,:);rotang=atan2(Y(ltIp1,:)-Y(ltI,:),X(ltIp1,:)-X(ltI,:))*180/pi+90;
  end
 end

 if strcmp(gbox,'fancy')
    if gtickdir(1)=='i'
      fancybox(lt,MAP_VAR_LIST.lats,'ygrid','left',dpatch,gticklen,gtickstyle); 
      drawticks=0;
    else    
      fancybox2(lt,MAP_VAR_LIST.lats,'ygrid','left',dpatch,gticklen,gtickstyle); 
    end
 end    
 if drawticks
   [n,m]=size(ltx);
   line(reshape([ltx;NaN+ones(1,m)],(n+1)*m,1),reshape([lty;NaN+ones(1,m)],(n+1)*m,1),...
        'linestyle','-','color',gcolor,'linewidth',glinewidth,'tag','m_grid_yticks-left','clipping','off');
   line(reshape([utx;NaN+ones(1,m)],(n+1)*m,1),reshape([uty;NaN+ones(1,m)],(n+1)*m,1),...
        'linestyle','-','color',gcolor,'linewidth',glinewidth,'tag','m_grid_yticks-right','clipping','off');
 end

 % Finally - the labels!
 ik=1:size(X,2);

 for k=ik
   [rotang(k), horizk, vertk] = upright(rotang(k), horiz, vert);
   text(xx(k),yy(k),labs{k},'horizontalalignment',horizk,'verticalalignment',vertk,...
        'rotation',rotang(k),'fontsize',gfontsize*scl(k),'color',gcolor,...
	'tag','m_grid_yticklabel','fontname',gfontname);
 end

end

% Draw the plot box


if strcmp(gbox,'on')
  [X,Y]=feval(MAP_PROJECTION.routine,'box');
  line(X(:),Y(:),'linestyle','-','linewidth',glinewidth,'color',gcolor,'tag','m_grid_box','clipping','off');
end

% Give a 1-1 aspect ratio and get rid of the matlab-provided axes stuff.


if isempty(strfind(version,'R2013b'))  % 27/Sept/13 - Handling for 2013b provided by CB.
set(gca,'visible','off',...
        'dataaspectratio',[1 1 1],...
        'xlim',MAP_VAR_LIST.xlims,...
        'ylim',MAP_VAR_LIST.ylims);
else
set(gca,'visible','off',...
        'dataaspectratio',[1 1 1e16],...
        'xlim',MAP_VAR_LIST.xlims,...
        'ylim',MAP_VAR_LIST.ylims);

end


set(get(gca,'title'),'visible','on');
set(get(gca,'xlabel'),'visible','on');
set(get(gca,'ylabel'),'visible','on');

% Set coordinate system back

m_coord(Currentmap.name);


%-------------------------------------------------------------
% upright simply turns tick labels right-side up while leaving
% their positions unchanged.
% Sat  98/02/21 Eric Firing
%
function   [rotang, horiz, vert] = upright(rotang, horiz, vert);
if isnan(rotang), rotang=0; end   % Added in 2014!
if rotang > 180, rotang = rotang - 360; end
if rotang < -180, rotang = rotang + 360; end
if rotang > 90
   rotang = rotang - 180;
elseif rotang < -90
   rotang = 180 + rotang;
else
   return    % no change needed.
end
switch horiz(1)
   case 'l'
      horiz = 'right';
   case 'r'
      horiz = 'left';
end
switch vert(1)
   case 't'
      vert = 'bottom';
   case 'b'
      vert = 'top';
end
  

%--------------------------------------------------------------------------
function [L,fs]=m_labels(dir,vals,uservals,tickstyle)
% M_LONLABEL creates longitude labels
%         Default values are calculated automatically when the grid is 
%         generated. However, the user may wish to specify the labels
%         as either numeric values or as strings (in the usual way
%         for axes).
%
%         If auto-labelling occurs, minutes are labelled in a different
%         (smaller) fontsize than even degrees.

% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997

% If the user has specified [] (i.e. no labels), we return blanks.

if isempty(uservals)
  L=cellstr(char(' '*ones(length(vals),1)));
  fs=1.0*ones(length(L),1);
  return;
end

% If the user has specified strings, we merely need to make
% sure that there are enough to cover all ticks.

if any(ischar(uservals)) 
  L=cellstr( uservals((rem([0:length(vals)-1],length(uservals))+1),:) );
  fs=1.0*ones(length(L),1);
  return;
end

% Otherwise we are going to have to generate labels from numeric
% data.

if length(uservals)==1 && isnan(uservals)  % use default values
  vals=vals(:)'; % make it a row.
else                                       % or ones provided
  lv=length(vals);
  vals=uservals(:)';
  while length(vals)<lv
    vals=[vals uservals(:)'];
  end
end

% longitudes and latitudes have some differences....
if findstr(dir,'lat') 
  labname=['S';'N';' '];
else
  labname=['W';'E';' '];
  vals=rem(vals+540,360)-180;
end

i=[vals<0;vals>0;vals==0];  % get the 'names' (i.e. N/S or E/W)
vals=abs(vals);             % Convert to +ve values

L=cell(length(vals),1);
fs=ones(length(vals),1);

if strcmp(tickstyle,'dm')

   % For each label we have different options:
   %  1 - even degrees are just labelled as such.
   %  2 - ticks that fall on even minutes are just labelled as even minutes
   %      in a smaller fontsize.
   %  3 - fractional minutes are labelled to 2 decimal places in the
   %      smaller fontsize.
   for k=1:length(vals)
     if rem(vals(k),1)==0
       nam=find(i(:,k));
       L{k}=sprintf([' %3.0f^o' labname(nam) ' '],vals(k));
     elseif abs(vals*60-round(vals*60))<.01
       L{k}=sprintf([' %2.0f'' '],rem(vals(k),1)*60);
       fs(k)=0.75;
     else
       L{k}=sprintf([' %2.2f'' '],rem(vals(k),1)*60);
       fs(k)=0.75;
     end
   end

   % In most cases, the map will have at least one tick with an even degree label,
   % but for very small regions (<1 degree in size) this won't happen so we
   % want to force one label to show degrees *and* minutes.

   if ~any(fs==1)  
    k=round(length(vals)/2);
    nam=find(i(:,k));
    L{k}={sprintf([' %3.0f^o' labname(nam) ' '],fix(vals(k))),...
	  sprintf([' %2.2f'' '],rem(vals(k),1)*60)};
    fs(k)=1;
   end

elseif strcmp(tickstyle,'dd')
   % For each label we have different options:
   %  1 - even degrees are just labelled as such.
   %  2 - 2 decimal place intervals use 2 decimal places
   %  3 - the rest fo to 4
   for k=1:length(vals)
     if rem(vals(k),1)==0
       nam=find(i(:,k));
       L{k}=sprintf([' %3.0f^o' labname(nam) ' '],vals(k));
     elseif abs(vals*100-round(vals*100))<0.01
       L{k}=sprintf([' %2.2f'],vals(k));
       fs(k)=0.75;
     else
       L{k}=sprintf([' %6.4f'],vals(k));
       fs(k)=0.75;
     end
   end

  
  % write code.
end


%---------------------------------------------------------
function [ltx,lty,utx,uty]=maketicks(X,Y,gticklen,gtickdir)
% MAKETICKS makes the axis ticks.
%           AXes ticks are based on making short lines at
%           the end of the grid lines X,Y.


% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997

global MAP_VAR_LIST

tlen=gticklen*max( diff(MAP_VAR_LIST.xlims),diff(MAP_VAR_LIST.ylims));

lx=sqrt((X(2,:)-X(1,:)).^2+(Y(2,:)-Y(1,:)).^2);

if strcmp(gtickdir,'out')
  ltx=[X(1,:)-tlen*(X(2,:)-X(1,:))./lx;X(1,:)];
  lty=[Y(1,:)-tlen*(Y(2,:)-Y(1,:))./lx;Y(1,:)];
else
  ltx=[X(1,:);X(1,:)+tlen*(X(2,:)-X(1,:))./lx];
  lty=[Y(1,:);Y(1,:)+tlen*(Y(2,:)-Y(1,:))./lx];
end

lx=sqrt((X(end,:)-X(end-1,:)).^2+(Y(end,:)-Y(end-1,:)).^2);

if strcmp(gtickdir,'out')
  utx=[X(end,:)-tlen*(X(end-1,:)-X(end,:))./lx;X(end,:)];
  uty=[Y(end,:)-tlen*(Y(end-1,:)-Y(end,:))./lx;Y(end,:)];
else
  utx=[X(end,:);X(end,:)+tlen*(X(end-1,:)-X(end,:))./lx];
  uty=[Y(end,:);Y(end,:)+tlen*(Y(end-1,:)-Y(end,:))./lx];
end


%---------------------------------------------------------
function fancybox(vals,lims,gridarg1,gridarg2,dpatch,gticklen,gridarg3);
%
%  FANCYBOX  - draws fancy outlines for either top/bottom or left/right sides,
%              depending on calling parameters.

global MAP_PROJECTION

% Get xlocations including endpoints
xval=sort([lims(1) vals(vals>lims(1) & vals<lims(2)) lims(2)]);
% Add all half-way points as well.
xval=sort([xval,xval(1:end-1)+diff(xval)/2]);
	
% Interpolate extra points to handle curved boundary conditions.
	
xval=(xval(1:end-1)'*ones(1,dpatch)+diff(xval)'*[0:dpatch-1]/dpatch)';
xval=[xval(:);lims(2)];
	
% Get lat/long positions for everything
	
[X2,Y2,lg2,lgI2]=feval(MAP_PROJECTION.routine,gridarg1,xval,gridarg2,gridarg3);
[l2x,l2y,u2x,u2y]=maketicks(X2,Y2,gticklen,'in');

if gridarg1(1)=='x', sig=1; else sig=-1; end

id=[1:dpatch size(l2x,2)+[-dpatch+1:0]];
dx=diff(l2x(:,id));
l2x(2,id)=l2x(2,id)+diff(l2y(:,id)).*([dpatch:-1:1 -1:-1:-dpatch]/(dpatch))*sig;
l2y(2,id)=l2y(2,id)-dx.*([dpatch:-1:1 -1:-1:-dpatch]/(dpatch))*sig;
dx=diff(u2x(:,id));
u2x(2,id)=u2x(2,id)-diff(u2y(:,id)).*([dpatch:-1:1 -1:-1:-dpatch]/(dpatch))*sig;
u2y(2,id)=u2y(2,id)+dx.*([dpatch:-1:1 -1:-1:-dpatch]/(dpatch))*sig;

% Now actually draw the patches.

% Added the z-values 16/Oct/05.
 
px=prod(size(l2x));
kk=[0:(dpatch*4):px-3]'*ones(1,dpatch*2+2);
kk=kk+ones(size(kk,1),1)*[1 2:2:(dpatch*2+2) (dpatch*2+1):-2:3];
patch(reshape(u2x(kk),size(kk,1),size(kk,2))',...
      reshape(u2y(kk),size(kk,1),size(kk,2))',...
      repmat(MAP_PROJECTION.LARGVAL  ,size(kk,2),size(kk,1)),'w','edgecolor','k','clipping','off','tag','m_grid_fancybox1');
patch(reshape(l2x(kk),size(kk,1),size(kk,2))',...
      reshape(l2y(kk),size(kk,1),size(kk,2))',...
      repmat(MAP_PROJECTION.LARGVAL-1,size(kk,2),size(kk,1)),'k','clipping','off','tag','m_grid_fancybox1');

kk=[dpatch*2:(dpatch*4):px-3]'*ones(1,dpatch*2+2);
kk=kk+ones(size(kk,1),1)*[1 2:2:(dpatch*2+2) (dpatch*2+1):-2:3];
patch(reshape(l2x(kk),size(kk,1),size(kk,2))',...
      reshape(l2y(kk),size(kk,1),size(kk,2))',...
      repmat(MAP_PROJECTION.LARGVAL  ,size(kk,2),size(kk,1)),'w','edgecolor','k','clipping','off','tag','m_grid_fancybox1');
patch(reshape(u2x(kk),size(kk,1),size(kk,2))',...
      reshape(u2y(kk),size(kk,1),size(kk,2))',...
      repmat(MAP_PROJECTION.LARGVAL-1,size(kk,2),size(kk,1)),'k','clipping','off','tag','m_grid_fancybox1');


%---------------------------------------------------------
function fancybox2(vals,lims,gridarg1,gridarg2,dpatch,gticklen,gridarg3)
%
%  FANCYBOX  - draws fancy outlines for either top/bottom or left/right sides,
%              depending on calling parameters.

global MAP_PROJECTION
 
% Get xlocations including endpoints
xval=sort([lims(1) vals(vals>lims(1) & vals<lims(2)) lims(2)]);
% Add all half-way points as well.
xval=sort([xval,xval(1:end-1)+diff(xval)/2]);
	
% Interpolate extra points to handle curved boundary conditions.
	
xval=(xval(1:end-1)'*ones(1,dpatch)+diff(xval)'*[0:dpatch-1]/dpatch)';
xval=[xval(:);lims(2)];
	
% Get lat/long positions for everything
	
[X2,Y2,lg2,lgI2]=feval(MAP_PROJECTION.routine,gridarg1,xval,gridarg2,gridarg3);
[l2x,l2y,u2x,u2y]=maketicks(X2,Y2,gticklen,'in');
	
if gridarg1(1)=='x', sig=1; else sig=-1; end

id=[1:dpatch size(l2x,2)+[-dpatch+1:0]];
dx=diff(l2x(:,id));
l2x(2,id)=l2x(2,id)+diff(l2y(:,id)).*([dpatch:-1:1 -1:-1:-dpatch]/(dpatch))*sig;
l2y(2,id)=l2y(2,id)-dx.*([dpatch:-1:1 -1:-1:-dpatch]/(dpatch))*sig;
dx=diff(u2x(:,id));
u2x(2,id)=u2x(2,id)-diff(u2y(:,id)).*([dpatch:-1:1 -1:-1:-dpatch]/(dpatch))*sig;
u2y(2,id)=u2y(2,id)+dx.*([dpatch:-1:1 -1:-1:-dpatch]/(dpatch))*sig;

% Now actually draw the patches.

% Added large z-values 16/Oct/05
 
px=prod(size(l2x));
kk=[0:(dpatch*2):px-3]'*ones(1,dpatch*2+2);
kk=kk+ones(size(kk,1),1)*[1 2:2:(dpatch*2+2) (dpatch*2+1):-2:3];
 patch(reshape(l2x(kk),size(kk,1),size(kk,2))',...
       reshape(l2y(kk),size(kk,1),size(kk,2))',...
       repmat(MAP_PROJECTION.LARGVAL-1,size(kk,2),size(kk,1)),...
       'w','edgecolor','k','clipping','off','linewidth',.2,'tag','m_grid_fancybox2');
 patch(reshape(u2x(kk),size(kk,1),size(kk,2))',...
       reshape(u2y(kk),size(kk,1),size(kk,2))',...
       repmat(MAP_PROJECTION.LARGVAL-1,size(kk,2),size(kk,1)),...
       'w','edgecolor','k','clipping','off','linewidth',.2,'tag','m_grid_fancybox2');

kk=[0:(dpatch*2):size(l2x,2)-dpatch-1]'*ones(1,dpatch+1);
kk=(kk+ones(size(kk,1),1)*[1:dpatch+1])';
[k1,k2]=size(kk);
line(reshape(mean(l2x(:,kk)),k1,k2),reshape(mean(l2y(:,kk))',k1,k2),...
     repmat(MAP_PROJECTION.LARGVAL,k1,k2),'color','k','clipping','off','tag','m_grid_fancybox2');

kk=[dpatch:(dpatch*2):size(l2x,2)-dpatch-1]'*ones(1,dpatch+1);
kk=(kk+ones(size(kk,1),1)*[1:dpatch+1])';
[k1,k2]=size(kk);
line(reshape(mean(u2x(:,kk))',k1,k2),reshape(mean(u2y(:,kk))',k1,k2),...
     repmat(MAP_PROJECTION.LARGVAL,k1,k2),'color','k','clipping','off','tag','m_grid_fancybox2');






