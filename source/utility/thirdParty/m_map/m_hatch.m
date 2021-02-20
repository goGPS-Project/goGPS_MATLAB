function [xi,yi,x,y]=m_hatch(lon,lat,varargin)
% M_HATCH Draws hatched or speckled interiors to a patch
%       
%    M_HATCH(LON,LAT,STYL,ANGLE,STEP,...line parameters);
%
% INPUTS:
%     X,Y - vectors of points.
%     STYL - style of fill
%     ANGLE,STEP - parameters for style
%
%     E.g.
%                 
%      'single',45,5  - single cross-hatch, 45 degrees,  5 points apart 
%      'cross',40,6   - double cross-hatch at 40 and 90+40, 6 points apart
%      'speckle',7,1  - speckled (inside) boundary of width 7 points, density 1
%                               (density >0, .1 dense 1 OK, 5 sparse)
%      'outspeckle',7,1 - speckled (outside) boundary of width 7 points, density 1
%                               (density >0, .1 dense 1 OK, 5 sparse)
%
%     
%      H=M_HATCH(...) returns handles to hatches/speckles.
%
%      [XI,YI,X,Y]=MHATCH(...) does not draw lines - instead it returns
%      vectors XI,YI of the hatch/speckle info, and X,Y of the original
%      outline modified so the first point==last point (if necessary).
%
%     Note that inside and outside speckling are done quite differently
%     and 'outside' speckling on large coastlines can be very slow.
%
%     If you get weird results - try putting an M_LINE(LON,LAT) call *before*
%     (or otherwise set the plot axis xlim/ylim parameters - this is necessary
%     because otherwise M_HATCH can't properly determine the 'points' units).
%

%
% Hatch Algorithm originally by K. Pankratov, with a bit stolen from 
% Iram Weinsteins 'fancification'. Speckle modifications by R. Pawlowicz.
%
% R Pawlowicz (rich@ocgy.ubc.ca) 15/Dec/2005
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.
%
%% 4/DEc/11 - isstr to ischar
%  Apr/12  - handle NaN-separated coastlines

global MAP_PROJECTION


if isempty(MAP_PROJECTION)
   disp('No Map Projection initialized - call M_PROJ first!');
   return;
end


styl='speckle';
angle=7;
step=1/2;

if ~isempty(varargin) && ischar(varargin{1})
  styl=varargin{1};
  varargin(1)=[];  
end
if ~isempty(varargin) && ~ischar(varargin{1})
  angle=varargin{1};
  varargin(1)=[];  
end
if ~isempty(varargin) && ~ischar(varargin{1})
  step=varargin{1};
  varargin(1)=[];
end



% If we have a naN-separated multi-line input vector:

ii=find(isnan(lon));

if any(ii)

  % if there isn't a NaN to mark the first or last segments
  if ii(1)>1, ii=[1;ii(:)]; end
  if ii(end)<length(lon), ii=[ii(:);length(lon)]; end
  
  % Compute the info, but don't draw - otherwise we end up with
  % lots of 'childred' to the plot and it is SLOWWWWW
  
  xi=[];yi=[];x=[];y=[];
  for k=1:length(ii)-1
    [xiT,yiT,xT,yT]=m_hatch(lon(ii(k)+1:ii(k+1)-1),lat(ii(k)+1:ii(k+1)-1),styl,angle,step,varargin{:});
    xi=[xi,NaN,xiT];
    yi=[yi,NaN,yiT];
    x=[x,NaN,xT];
    y=[y,NaN,yT];
  end
  
  % OK, so now plot it all.
  if nargout<2
    switch lower(styl)
      case {'single','cross'}
         xi=line(xi,yi,varargin{:});
      case {'speckle','outspeckle'}
         xi=line(xi,yi,'marker','.','linestyle','none','markersize',2,varargin{:});
    end
  end 
  
  return;
end
      
%----------------------
% otherwise, handle a single line without NaN

   
[x,y,I]=m_ll2xy(lon,lat,'clip','patch');
 
%% plot(x,y,'color','r');%%pause;
  
if x(end)~=x(1) || y(end)~=y(1)  % & to |
  x=x([1:end 1]);
  y=y([1:end 1]);
  I=I([1:end 1]);
end

if strcmp(styl,'speckle') || strcmp(styl,'outspeckle')
  angle=angle*(1-I);
end

if size(x,1)~=1
 x=x(:)';
 angle=angle(:)';
end
if size(y,1)~=1
 y=y(:)';
end
 

% Code stolen from Weinstein hatch
oldu = get(gca,'units');
set(gca,'units','points');
sza = get(gca,'position'); sza = sza(3:4);
set(gca,'units',oldu)   % Set axes units back

xlim = get(gca,'xlim');
ylim = get(gca,'ylim');
xsc = sza(1)/(xlim(2)-xlim(1)+eps);
ysc = sza(2)/(ylim(2)-ylim(1)+eps);

switch lower(styl)
 case 'single'
  [xi,yi]=drawhatch(x,y,angle,step,xsc,ysc,0);
  if nargout<2
    xi=line(xi,yi,varargin{:});
  end 
 case 'cross'
  [xi,yi]=drawhatch(x,y,angle,step,xsc,ysc,0);
  [xi2,yi2]=drawhatch(x,y,angle+90,step,xsc,ysc,0);
  xi=[xi,xi2];
  yi=[yi,yi2];
  if nargout<2
    xi=line(xi,yi,varargin{:});
  end 
 case 'speckle'
  [xi,yi ]  =drawhatch(x,y,45,   step,xsc,ysc,angle);
  [xi2,yi2 ]=drawhatch(x,y,45+90,step,xsc,ysc,angle);
  xi=[xi,xi2];
  yi=[yi,yi2];
  if nargout<2
    if any(xi)
      xi=line(xi,yi,'marker','.','linestyle','none','markersize',2,varargin{:});
    else
      if ~MAP_PROJECTION.newgraphics
        xi=NaN;yi=NaN;
      else
        xi=gobjects(1);yi=gobjects(1);
      end		
    end    
  end 
 case 'outspeckle'
  [xi,yi ]  =drawhatch(x,y,45,   step,xsc,ysc,-angle);
  [xi2,yi2 ]=drawhatch(x,y,45+90,step,xsc,ysc,-angle);
  xi=[xi,xi2];
  yi=[yi,yi2];
%  disp([size(xi) size(yi)])
  inside=logical(inpolygon(xi,yi,x,y)); % logical needed for v6!
  xi(inside)=[];yi(inside)=[];
  if nargout<2
    if any(xi)
      xi=line(xi,yi,'marker','.','linestyle','none','markersize',2,varargin{:});
    else
      if MAP_PROJECTION.newgraphics
         xi=gobjects(1);yi=gobjects(1);
      else
         xi=NaN;yi=NaN;
      end		
    end   
  end
    
end


return

%%%%%

function [xi,yi]=drawhatch(x,y,angle,step,xsc,ysc,speckle)
%
% This is the guts. 
%
 

angle=angle*pi/180;

% Idea here appears to be to rotate everthing so lines will be
% horizontal, and scaled so we go in integer steps in 'y' with
% 'points' being the units in x.
% Center it for "good behavior".
ca = cos(angle); sa = sin(angle);
x0 = mean(x); y0 = mean(y);   
x = (x-x0)*xsc; y = (y-y0)*ysc;
yi = x*ca+y*sa;              % Rotation
y = -x*sa+y*ca;
x = yi;
y = y/step;    % Make steps equal to one

 
% Compute the coordinates of the hatch line ...............
yi = ceil(y);
yd = [diff(yi) 0]; % when diff~=0 we are crossing an integer
fnd = find(yd);    % indices of crossings
dm = max(abs(yd)); % max possible #of integers between points
 

%
% This is going to be pretty space-inefficient if the line segments
% going in have very different lengths. We have one column per line
% interval and one row per hatch line within that interval.
%
A = cumsum( repmat(sign(yd(fnd)),dm,1), 1);

% Here we interpolate points along all the line segments at the
% correct intervals.
fnd1 = find(abs(A)<=abs( repmat(yd(fnd),dm,1) ));
A  = A+repmat(yi(fnd),dm,1)-(A>0);
xy = (x(fnd+1)-x(fnd))./(y(fnd+1)-y(fnd));
xi = repmat(x(fnd),dm,1)+(A-repmat(y(fnd),dm,1) ).*repmat(xy,dm,1);
yi = A(fnd1);
xi = xi(fnd1);


 % Sorting points of the hatch line ........................
%%%yi0 = min(yi); yi1 = max(yi);
% Sort them in raster order (i.e. by x, then by y)
% Add '2' to make sure we don't have problems going from a max(xi)
% to a min(xi) on the next line (yi incremented by one)
if ~isempty(xi)
  xi0 = min(xi); xi1 = max(xi);
  ci = 2*yi*(xi1-xi0)+xi;
  [ci,num] = sort(ci);
  xi = xi(num); yi = yi(num);
else
 
  xi=NaN;yi=NaN;
  return;
end

% if this happens an error has occurred somewhere (we have an odd
% # of points), and the "fix" is not correct, but for speckling anyway
% it really doesn't make a difference.
if rem(length(xi),2)==1
  disp('mhatch warning');
  xi = [xi; xi(end)];
  yi = [yi; yi(end)];
end

 % Organize to pairs and separate by  NaN's ................
li = length(xi);
xi = reshape(xi,2,li/2);
yi = reshape(yi,2,li/2);

% The speckly part - instead of taking the line we make a point some
% random distance in.
if length(speckle)>1 || speckle(1)~=0

 if length(speckle)>1
   % Now we get the speckle parameter for each line.
   
   % First, carry over the speckle parameter for the segment
%   yd=[0 speckle(1:end-1)];
   yd=[speckle(1:end)];
   A=repmat(yd(fnd),dm,1);
   speckle=A(fnd1);
   
   % Now give it the same preconditioning as for xi/yi
   speckle=speckle(num);
   if rem(length(speckle),2)==1
     speckle = [speckle; speckle(end)];
   end
   speckle=reshape(speckle,2,li/2);

 else
   speckle=[speckle;speckle];
 end
   
 % Thin out the points in narrow parts.
 % This keeps everything when abs(dxi)>2*speckle, and then makes
 % it increasingly sparse for smaller intervals.
 oldxi=xi;oldyi=yi;
 dxi=diff(xi);
 nottoosmall=sum(speckle,1)~=0 & rand(1,li/2)<abs(dxi)./(max(sum(speckle,1),eps));
 xi=xi(:,nottoosmall);
 yi=yi(:,nottoosmall);
 dxi=dxi(nottoosmall);
 if size(speckle,2)>1, speckle=speckle(:,nottoosmall); end
 % Now randomly scatter points (if there any left)
 li=length(dxi);
 if any(li)
   xi(1,:)=xi(1,:)+sign(dxi).*(1-rand(1,li).^0.5).*min(speckle(1,:),abs(dxi) );
   xi(2,:)=xi(2,:)-sign(dxi).*(1-rand(1,li).^0.5).*min(speckle(2,:),abs(dxi) );
   % Remove the 'zero' speckles
   if size(speckle,2)>1
    xi=xi(speckle~=0);
    yi=yi(speckle~=0);
   end
 end
  
else
 xi = [xi; ones(1,li/2)*nan];  % Separate the line segments
 yi = [yi; ones(1,li/2)*nan];
end
xi = xi(:)'; yi = yi(:)';

% Transform back to the original coordinate system
yi = yi*step;
xy = xi*ca-yi*sa;
yi = xi*sa+yi*ca;
xi = xy/xsc+x0;
yi = yi/ysc+y0;


