function h=m_range_ring(long,lat,range,varargin)
% M_RANGE_RING Creates range rings on a map
%    M_RANGE_RING(LONG,LAT,RANGE) creates a range ring of range RANGE
%    km centered at the position specified by LONG and LAT. Range rings
%    will generally appear as small (almost) circles for short ranges,
%    but will be distorted at longer ranges.
%
%    If RANGE is a vector, concentric rings at the specified ranges
%    are drawn. If LONG,LAT are vectors, rings are drawn around
%    all specified locations.
%
%    The appearance of lines can be modified using the usual
%    line properties thus:
%    M_RANGE_RING(LONG,LAT,RANGE, <line property/value pairs>)
%    
%    Sometimes you may need to adjust the number of points plotted
%    in each range ring (this can happen if the ring is at the extreme
%    edge of certian projections). THis can be done using
%    M_RANGE_RING(LONG,LAT,RANGE,NPTS, <line property/value pairs>)
%
%    NB: Earth radius is assumed to 6378.137km (WGS-84 value), and
%    calculations are for spherical geometry.

% Rich Pawlowicz (rich@ocgy.ubc.ca) 18/Dec/1998
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

% 6/Nov/00 - eliminate returned stuff if ';' neglected (thx to D Byrne)
% 7/Dec/11 - Octave 3.2.3 compatibility

global MAP_PROJECTION MAP_VAR_LIST


if isempty(MAP_PROJECTION)
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end


pi180=pi/180;
earth_radius=6378.137;
n=72;

if length(varargin)>0 && ~ischar(varargin{1})
 n=varargin{1};varargin(1)=[];
end

 
 

c=range(:)'/earth_radius;

h=[];
for k=1:length(long)
  rlat=lat(k)*pi180;
  rlong=long(k)*pi180;
  if long(k)<MAP_VAR_LIST.longs(1), rlong=rlong+2*pi; end
  if long(k)>MAP_VAR_LIST.longs(2), rlong=rlong-2*pi; end

  x=sin([0:n-1]'/(n-1)*2*pi)*c;
  y=cos([0:n-1]'/(n-1)*2*pi)*c;
  on=ones(n,1);

  Y=(asin(on*cos(c)*sin(rlat) + (on*cos(rlat)*(sin(c)./c)).*y))/pi180;
  switch lat(k)
    case 90
      X=(rlong+atan2(x,-y))/pi180;
    case -90
      X=(rlong+atan2(x,y))/pi180;
    otherwise
      X=(rlong+atan2(x.*(on*sin(c)),on*(cos(rlat)*cos(c).*c) - (on*sin(rlat)*sin(c)).*y ) )/pi180;
  end
 
  nz=zeros(1,length(range(:)));
  X=X+cumsum([nz;diff(X)<-300]-[nz;diff(X)>300])*360;

  kk=find(X(1,:)~=X(end,:));
  X2=X(:,kk)+360;X2(X2>MAP_VAR_LIST.longs(2))=NaN;
  X3=X(:,kk)-360;X3(X3<MAP_VAR_LIST.longs(1))=NaN;
  
  [XX,YY]=m_ll2xy([X,X2,X3],[Y,Y(:,kk),Y(:,kk)],'clip','on');
  
  % Get rid of 2-point lines (these are probably clipped lines spanning the window)
  fk=isfinite(XX(:));        
  st=find(diff(fk)==1)+1;
  ed=find(diff(fk)==-1);
  if length(st)<length(ed), st=[1;st]; end
  if length(ed)<length(st), ed=[ed;length(fk)]; end
  k2=find((ed-st)==1);
  XX(st(k2))=NaN;
 
  if MAP_PROJECTION.IsOctave
     for k2=1:size(XX,2)
         h=[h;line(XX(:,k2),YY(:,k2),varargin{:},'tag','m_range_ring')];
     end
  else   
     h=[h;line(XX,YY,varargin{:},'tag','m_range_ring')];
  end
  
end

if nargout==0
 clear h
end
