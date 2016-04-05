function [h,varargout]=m_ellipse(long,lat,fmaj,fmin,finc,fpha,scl,tpe,varargin)
% M_ELLIPSE Draws ellipses on a map.
%      [H]=M_ELLIPSE(LONG,LAT,MAJOR,MINOR,INC,PHA,SCALE,TYPE)
%      draws ellipses as patch objects or lines on a map created 
%      by the M_Map package.
%
%      H:  vector of handles to the patch or line object
%      LONG, LAT: coordinates for the ellipse center
%      MAJOR:     ellipse major axis
%      MINOR:     ellipse minor axis
%      INC:       ellipse inclination (degrees CCW of East)
%      PHA:       phase in degrees (useful for tidal ellipses)
%                 set to [] when not needed (will also omit phase lines)
%      SCALE:     scale factor (set to [] for default scale)
%      TYPE:      'patch' or 'line'
%      VARARGIN   options passed to the patch or line objects
%       
%      By default, drawing the ellipses as patch objects with phase
%      information will produce an ellipse using the full range
%      colours of the current colormap.  As time increases around the
%      ellipse, the colormap values change from, for example, blue to
%      red (jet) or black to white (grayscale).  If no phase info is
%      given, the patch colour can be specified with varargin (eg -
%      'facecolor','k').  If no colour is specified, the patch will be
%      blue.
%
%
% Mark Halverson (mhalvers@eos.ubc.ca) 25/03/2013

% Revisions:
% Rich Pawlowicz (rich@eos.ubc.ca)
%       28/03/2013
%            - vectorized ellipse parameter inputs
%            - coordinate transformation for maps covering large areas
%            - phase lines added to 'line' option
%            - phase line colours change with rotation direction
% Mark Halverson (mhalvers@eos.ubc.ca) 
%       28/03/2013  
%            - added crude input error handling


global MAP_PROJECTION MAP_VAR_LIST
 
if isempty(MAP_PROJECTION),
  error('No Map Projection initialized - call M_PROJ first!');
  return;
end;

if nargin<8,
    error('m_ellipse requires at least 8 input parameters')
    return;
end

if isempty(scl),
  scl=1;    
end
scl=scl*1/100/60*2;

% flag to determine how to deal with phase info 
pha_flag=1;  %draw phase lines for ellipse type 'line'

% if phase information is not given - take it to be 0

if isempty(fpha),
  fpha=zeros(size(fmaj));
  pha_flag=0;  % do not draw phase lines or phase patch colours 
end

% if any(finc<0) | any(finc>180),
%     error('Ellipse inclination must be between 0 and 180')
%     return;
% end

t=[0:1/24:1 NaN]*360;  % 0 to 360 for drawing ellipses
%t=[0:.1/24:1 NaN]*360;  % 0 to 360 for drawing ellipses                        
%t=[0:1/24:0.958333 NaN]*360;  % 0 to 360 for drawing ellipses

lt=length(t);

% ellipse centers in m_map (X,Y) coords
[X,Y]=m_ll2xy(long(:),lat(:),'clip','point');
X2=repmat(X,1,lt);
Y2=repmat(Y,1,lt);

fmaj2=repmat(fmaj(:)*scl,1,lt);
finc2=repmat(finc(:),1,lt);
fpha2=repmat(fpha(:),1,lt);
fmin2=repmat(fmin(:)*scl,1,lt);
t2=repmat(t,size(fmaj2,1),1);
 
% parametric equation for ellipse 
phang=t2-fpha2;  % Phase angle relative to greenwich
                 % A small positive greenwich phase should appear a
		 % small distance CW if the sense of rotation is CCW.

fmajc=fmaj2.*cosd(phang);
fmins=fmin2.*sind(phang);
cinc=cosd(finc2);
sinc=sind(finc2);
x=fmajc.*cinc-fmins.*sinc;
y=fmajc.*sinc+fmins.*cinc;

% Coordinate transformation
[XN ,YN ]=m_ll2xy([long(:) long(:)]',                           [lat(:) lat(:)+.001]','clip','off');
[XE ,YE ]=m_ll2xy([long(:) long(:)+(.001)./cos(lat(:)*pi/180)]',[lat(:) lat(:)]',     'clip','off');
mx=x.*repmat(diff(XE)',1,lt)*1000 + y.*repmat(diff(XN)',1,lt)*1000;
my=x.*repmat(diff(YE)',1,lt)*1000 + y.*repmat(diff(YN)',1,lt)*1000;


% coordinates for the ellipse
xe=(X2+mx)'; % center plus ellipse edge
ye=(Y2+my)'; 
 

% Draw each ellipse as a series of triangular patches, and shade each
% triangle according to phase.  The darkest colour (highest level) is
% current when phase at Greenwich = 0, earlier vectors are shaded
% gradually lighter.
%
% Each column of xp and yp define the vertices of a triangle.  Loop
% over times (triangles) instead of ellipses.

if strcmp(tpe,'patch'),

  for k=lt-2:-1:1,
  
     xp=[X';xe([k;k+1],:)];
     yp=[Y';ye([k;k+1],:)];

     h(k)=patch(xp,yp,ones(size(xp))*(lt-k).^(1/3),(k).^(1/3),varargin{:});
  end;   
  set(h,'edgecolor','none','clip','off')
  
  if pha_flag==0,  %no phase info, and no varargin given
      set(h,'facecolor','b',varargin{:})
  end
  
end

% draw ellipses as simple lines, with a line for phang(t=1) and
% phang(t=2)
if strcmp(tpe,'line'),

  h=line(xe(:),ye(:),varargin{:});
  
  if pha_flag==1,
    ii=fmin(:)>0; % Different colours for CW (same as ellipse colour) and CCW (opposite colour)
    if any(ii),
      hp(ii)=line([X(ii)';xe(1,ii)],[Y(ii)';ye(1,ii)],'linewi',2,'color',get(h,'color'));
      hp2(ii)=line([X(ii)';xe(2,ii)],[Y(ii)';ye(2,ii)],'linewi',2,'color',get(h,'color'),'linest',':');
    end;
    if any(~ii),
      hp(~ii)=line([X(~ii)';xe(1,~ii)],[Y(~ii)';ye(1,~ii)],'linewi',2,'color',[1 1 1]-get(h,'color'));
      hp2(~ii)=line([X(~ii)';xe(2,~ii)],[Y(~ii)';ye(2,~ii)],'linewi',2,'color',[1 1 1]-get(h,'color'),'linest',':');
    end;
    
    h=[h,hp,hp2];
  
  end %pha_flag
   
end %ellipse 'line'


set(h,'tag','m_ellipse')


if nargout==0,
 clear h
end;




