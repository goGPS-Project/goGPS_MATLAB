function varargout=mu_util(optn,varargin);
% MU_UTIL Various utility routines
%           This function should not be used directly; instead it is
%           is accessed by other high- and low-level functions.

% MU_UTIL is basically a driver for a number of lower-level routines collected
% together for convenience. They are specified by the first argument:
% 'clip'      - clipping routine
% 'axisticks' - generates a "nice" set of ticks, optimized for degrees/minutes on maps
%               angles of the circle.
% 'x/ygrid'   - generates the lines for axis grids
% 'xylimits'  - finds the x/y limits for a given map
% 'lllimits'  - finds the lat/long limits for a given map (usually for rectboxes)
% 'box'       - returns a line around the boundaries of the map.

% Rich Pawlowicz (rich@ocgy.ubc.ca) 4/April/97
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

% 31/Mar/04 - added a fix in m_rectgrid that caused problems when a map
% boundary coincided with a grid line.
% 26/Oct/07 - fixed the same problem when it occurred near SOUTH pole!
% 8/Sep/13 - added 'tickstyle' to grid generation

switch optn,
  case 'clip',
    [varargout{1},varargout{2}]=m_clip(varargin{:});
  case 'axisticks',
    varargout{1}=m_getinterval(varargin{:});
  case {'xgrid','ygrid'}
    [varargout{1},varargout{2},varargout{3},varargout{4}]=m_rectgrid(optn,varargin{:});
  case 'xylimits',
    m_getxylimits;
  case 'lllimits',
    m_getlllimits;
  case 'box',
    [varargout{1},varargout{2}]=m_box(varargin{:});

end;



%---------------------------------------------------------
function [Xc,Yc]=m_clip(cliptype,X,Xedge,indx,Y);
% M_CLIP performs clipping of data. Columns of points are
%        assumed to be lines; the first points outside the
%        clip area are recomputed to lie on the edge of the
%        region; others are converted to either NaN or
%        edge points depending on CLIPTYPE.
%
%        indx = 0 for points inside the clip region, 1
%               for those outside
%
%         'on'  - replaces points outside with NaN, but interpolates
%                 to provide points right on the border.
%         'patch' - replaces points outside with nearest border point
%         'point' - does no interpolation (just checks in/out)
%
% In general m_clip will be called 4 times, since it solves a 1-edge problem.

% Rich Pawlowicz (rich@ocgy.ubc.ca) 4/April/97

Xc=X;
Yc=Y;

if ~strcmp(cliptype,'point'),

  % Find regions where we suddenly come into the area
  % (indices go from 1 to 0)

  [i,j]=find(diff(indx)==-1);

  if any(i),
    I=i+(j-1)*size(X,1); % 1-d addressing

    % Linearly interpolate to boundary

    bt=(X(I+1)-X(I));
    ibt=abs(bt)<5*eps;
    if any(ibt), bt(ibt)=1*eps; end; % In these cases the delta(Y) also = 0, so we just want
                                      % to avoid /0 warnings.
    Yc(I)=Y(I)+(Xedge-X(I)).*(Y(I+1)-Y(I))./bt;
    Yc(I(ibt))=(Y(I(ibt))+Y(I(ibt)+1))/2;
    Xc(I)=Xedge;
    indx(I(isfinite(Yc(I))))=0;
  end;

  % Find regions where we suddenly come out of the area
  % (indices go from 0 to 1)

  
  [i,j]=find(diff(indx)==1);
  if any(i),

    I=i+(j-1)*size(X,1);

    bt=(X(I+1)-X(I));
    ibt=abs(bt)<5*eps;
    if any(ibt), bt(ibt)=eps;  end; % In these cases the delta(Y) also = 0, so we just want
                                      % to avoid /0 warnings.
  
    Yc(I+1)=Y(I)+(Xedge-X(I)).*(Y(I+1)-Y(I))./bt;
    Yc(I(ibt)+1)=(Y(I(ibt))+Y(I(ibt)+1))/2;
    Xc(I+1)=Xedge;
    indx(I(isfinite(Yc(I+1)))+1)=0;
  end;

end;

switch cliptype,
  case {'on','point'}
    Xc(indx)=NaN;
    Yc(indx)=NaN;
  case 'patch',
    Xc(indx)=Xedge;
end;


%--------------------------------------------------------------------------
function gval=m_getinterval(gmin,gmax,gtick,gtickstyle);
% M_GETINTERVAL picks nice spacing for grid ticks
%        This occurs when the following call is made:
%        TICKS=M_GRID('axisticks',MIN,MAX,APPROX_NUM_TICKS)
%

% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997
%
% 9/Apr/98 - changed things so that max/min limits are not automatically
%            added (this feature made map corners messy sometimes) 

% If ticks are specified, we just make sure they are within the limits
% of the map.

if length(gtick)>1,
  gval=[gtick(gtick(:)>=gmin & gtick(:)<=gmax)];

% Otherwise, we try to fit approximately gtick ticks in the interval
else

  if gtick>2,
    
    exactint=(gmax-gmin)/(gtick-1)*60; %interval in minutes

    if strcmp(gtickstyle,'dm'),
       % These are the intervals which we will allow (they are "nice" in the sense
       % that they come to various even multiples of minutes or degrees)
       niceints=[0.1 0.2 0.25 0.5 ...
        	 1 2 3 4 5 6 10 12 15 20 30 ...
        	 60*[1 2 3 4 5 6 8 9 10 12 15 18 20 25 30 40 50 60 100 120 180]];
    elseif strcmp(gtickstyle,'dd'),
       % these are decimal intervals
       niceints=60*[1/500 1/400 1/250 1/200 1/100 ...
		     1/50 1/40 1/25 1/20 1/10 1/5 1/4 1/3 1/2 ...
		     1 2 3 4 5 6 8 9 10 12 15 18 20 25 30 40 50 60 100 120 180];
    else
       error(['bad tickstyle - ''' gtickstyle '''']);
    end;
    
 
    [dun,I]=min(abs(niceints-exactint));

    gval=niceints(I)/60*[ceil(gmin*60/niceints(I)):fix(gmax*60/niceints(I))];
    
    gval=[gval(gval>=gmin & gval<=gmax) ];
  else
    gval=[gmin gmax];
  end;
end;

 
%--------------------------------------------------------------
function [X,Y,vals,labI]=m_rectgrid(direc,Xlims,Ylims,Nx,Ny,label_pos,tickstyle);
% M_RECTGRID This handles some of the computations involved in creating grids
%            for rectangular maps. Essentially we make our "first guess" using the
%            lat/long limits. Then these curves are clipped to the boundaries, after
%            which we use the clip points as "new" boundaries and recompute the lines.

% Rich Pawlowicz (rich@ocgy.ubc.ca) 4/April/97

global MAP_PROJECTION MAP_VAR_LIST

Ny1=Ny-1;
Ny2=Ny*2;
Ny21=Ny2-1;

% First try some wildly oversampled lines (not including the boundaries)

vals=mu_util('axisticks',Xlims(1),Xlims(2),Nx,tickstyle);

if strcmp(MAP_VAR_LIST.rectbox,'on') |  strcmp(MAP_VAR_LIST.rectbox,'circle'),

 % We don't want the end limits here.
 if vals(end) == Xlims(2), vals(end) = []; end
 if vals(1)   == Xlims(1), vals(1)   = []; end

 if direc(1)=='x',
  [lg,lt]=meshgrid(vals,Ylims(1)+diff(Ylims)*[0:1/Ny1:1]);
 else
  [lt,lg]=meshgrid(vals,Ylims(1)+diff(Ylims)*[0:1/Ny1:1]);
 end;

 % But sneakily we clip them in transforming, so we end up with finite values only
 % inside the axis limits

 [X,Y]=feval(MAP_PROJECTION.routine,'ll2xy',lg,lt,'clip','on');

 % Now we find the first/last unclipped values; these will be our correct starting points.
 % (Note I am converting to one-dimensional addressing).
 
 istart=sum(cumsum(isfinite(X))==0)+1+[0:size(X,2)-1]*size(X,1);
 iend=size(X,1)-sum(cumsum(isfinite(flipud(X)))==0)+[0:size(X,2)-1]*size(X,1);

 % Now, in the case where map boundaries coincide with limits it is just possible
 % that an entire column might be NaN...so in this case make up something just
 % slight non-zero. (31/Mar/04)
 
 i3=find(iend==0); 
 if any(i3), istart(i3)=1; iend(i3)=1; end; 

 % do same fix for istart (thanks Ben Raymond for finding this bug)
 i3=find(istart>prod(size(X)));
 if any(i3), istart(i3)=1; iend(i3)=1; end;
  
 % Now go back and find the lat/longs corresponding to those points; these are our new
 % starting points for the lines (Note that the linear interpolation for clipping at boundaries
 % means that they will not *quite* be the exact longitudes due to curvature, but they should
 % be very close.

 if direc(1)=='x',

  [lgs,lts]=feval(MAP_PROJECTION.routine,'xy2ll',X(istart),Y(istart),'clip','off');
  [lgs,lte]=feval(MAP_PROJECTION.routine,'xy2ll',X(iend),Y(iend),'clip','off');

  % Finally compute the lines within those limits (these *may* include some out-of-bounds points
  % depending on the geometry of the situation; these are converted to NaN as usual.

  [X,Y]=feval(MAP_PROJECTION.routine,'ll2xy',vals(ones(Ny2,1),:),...
              lts(ones(Ny2,1),:)+[0:1/Ny21:1]'*(lte-lts),'clip','on');

 else

  [lgs,lts]=feval(MAP_PROJECTION.routine,'xy2ll',X(istart),Y(istart),'clip','off');
  [lge,lts]=feval(MAP_PROJECTION.routine,'xy2ll',X(iend),Y(iend),'clip','off');

  % Longitudes should be increasing here, but we can run into wrap problems after
  % the tansformations back and forth. It is for this line that I need to make
  % long-lims just a tad less than 360 when really they should be 360 (in m_lllimits)

  lge(lge<=lgs)=lge(lge<=lgs)+360; 

  [X,Y]=feval(MAP_PROJECTION.routine,'ll2xy',lgs(ones(Ny2,1),:)+[0:1/Ny21:1]'*(lge-lgs),...
              vals(ones(Ny2,1),:),'clip','on');

 end;

else
 if direc(1)=='x',
  [lg,lt]=meshgrid(vals,Ylims(1)+diff(Ylims)*[0:1/Ny1:1]);
 else
  [lt,lg]=meshgrid(vals,Ylims(1)+diff(Ylims)*[0:1/Ny1:1]);
 end;
 [X,Y]=feval(MAP_PROJECTION.routine,'ll2xy',lg,lt,'clip','off');

end;

switch label_pos
  case {'left','bottom','west','south'}
    labI=1;
  case 'middle',
    labI=round(size(X,1)/2+1/2);
  case {'right','top','east','north'}
    labI=size(X,1);
end;



%--------------------------------------------------------------------------------
function m_getxylimits;
% M_GET_LIMITS Converts X/Y limits to lat/long limits
%              This is a chunk of code that is needed for most projections.

% Rich Pawlowicz (rich@ocgy.ubc.ca) 4/April/97

global MAP_PROJECTION MAP_VAR_LIST

% Start with the user-specified lat/longs.

MAP_VAR_LIST.lats=MAP_VAR_LIST.ulats;
MAP_VAR_LIST.longs=MAP_VAR_LIST.ulongs;

% Now, let's get the map x/ylims
bX=MAP_VAR_LIST.longs(1)+diff(MAP_VAR_LIST.longs)*[0:1/30:1];
bY=MAP_VAR_LIST.lats(1)+diff(MAP_VAR_LIST.lats)*[0:1/30:1];
bX=[bX MAP_VAR_LIST.longs(2*ones(1,31)) fliplr(bX) MAP_VAR_LIST.longs(ones(1,31)) ];
bY=[MAP_VAR_LIST.lats(ones(1,31)) bY MAP_VAR_LIST.lats(2*ones(1,31)) fliplr(bY) ];
[X,Y]=feval(MAP_PROJECTION.routine,'ll2xy',bX,bY,'clip','off');
MAP_VAR_LIST.xlims=[min(X) max(X)];
MAP_VAR_LIST.ylims=[min(Y) max(Y)];


%-------------------------------------------------------------------------------
function m_getlllimits;
% M_GET_LIMITS Converts X/Y limits to lat/long limits
%              This is a chunk of code that is needed for most projections.

% Rich Pawlowicz (rich@ocgy.ubc.ca) 4/April/97

global MAP_PROJECTION MAP_VAR_LIST

[bX,bY]=mu_util('box',31);
% Get its lat/longs.

[lg,lt]=feval(MAP_PROJECTION.routine,'xy2ll',bX,bY,'clip','off');
% Take real part because otherwise funny things might happen if the box is very large

MAP_VAR_LIST.lats=[min(real(lt)) max(real(lt))];

% Are the poles within the axis limits? (Test necessary for oblique mercator and azimuthal)
[px,py]=feval(MAP_PROJECTION.routine,'ll2xy',[0 0],[-90 90],'clip','point');
if isfinite(px(1)), MAP_VAR_LIST.lats(1)=-90; end;
if isfinite(px(2)), MAP_VAR_LIST.lats(2)= 90; end;

if any(isfinite(px)),
  MAP_VAR_LIST.longs=[-179.9 180]+exp(1); % we add a weird number (exp(1)) to get away from 
                         % anything that might conceivably be desired as a 
                         % boundary - it makes grid generation easier.
                         % Also make the limits just a little less than 180, this
                         % is necessary because I have to have the first and last points
                         % of lines just a little different in order to figure out orientation
                         % in 'm_rectgrid'
else
  MAP_VAR_LIST.longs=[min(lg) max(lg)];
  if all(isnan(px)) & diff(MAP_VAR_LIST.longs)>360*30/31,
    ii=lg<mean(MAP_VAR_LIST.longs);
    lg(ii)=lg(ii)+360;
    MAP_VAR_LIST.longs=[min(lg) max(lg)];
  end;
end;

%------------------------------------------------------------------------
function [X,Y]=m_box(npts);
% M_BOX  Computes coordinates of the map border.
%

% Rich Pawlowicz (rich@ocgy.ubc.ca) 4/April/97

global MAP_PROJECTION MAP_VAR_LIST

n1=npts-1;

switch MAP_VAR_LIST.rectbox,
  case 'on',
    X=MAP_VAR_LIST.xlims(1)+diff(MAP_VAR_LIST.xlims)*[0:1/n1:1];
    Y=MAP_VAR_LIST.ylims(1)+diff(MAP_VAR_LIST.ylims)*[0:1/n1:1];
    X=[X MAP_VAR_LIST.xlims(2*ones(1,npts)) fliplr(X) MAP_VAR_LIST.xlims(ones(1,npts))];
    Y=[MAP_VAR_LIST.ylims(ones(1,npts)) Y  MAP_VAR_LIST.ylims(2*ones(1,npts)) fliplr(Y)];
  case 'off',
    lg=MAP_VAR_LIST.longs(1)+diff(MAP_VAR_LIST.longs)*[0:1/n1:1];
    lg=[lg MAP_VAR_LIST.longs(2*ones(1,npts)) fliplr(lg) MAP_VAR_LIST.longs(ones(1,npts))]'; 
    lt=MAP_VAR_LIST.lats(1)+diff(MAP_VAR_LIST.lats)*[0:1/n1:1];
    lt=[MAP_VAR_LIST.lats(ones(1,npts)) lt  MAP_VAR_LIST.lats(2*ones(1,npts)) fliplr(lt)]';
    [X,Y]=feval(MAP_PROJECTION.routine,'ll2xy',lg,lt,'clip','off');
  case 'circle',
    n1=npts*3-1;
    X=MAP_VAR_LIST.rhomax*cos([0:n1]/n1*pi*2);
    Y=MAP_VAR_LIST.rhomax*sin([0:n1]/n1*pi*2);
end;
