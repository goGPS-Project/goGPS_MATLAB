function [dist,lons,lats] = m_lldist(long,lat,N)
% M_LLDIST Spherical earth distance between points in long/lat coordinates. 
%   RANGE=M_LLDIST(LONG,LAT) gives the distance in kilometers between
%   successive points in the vectors LONG and LAT, computed
%   using the Haversine formula on a spherical earth of radius
%   6378.137km. Distances are probably good to better than 1% of the
%   "true" distance on the ellipsoidal earth
%
%   [RANGE,LONGS,LATS]=M_LLDIST(LONG,LAT,N) computes the N-point geodesics
%   between successive points. Each geodesic is returned on its
%   own row of length N+1.
%
%   See also M_XYDIST

% Rich Pawlowicz (rich@ocgy.ubc.ca) 6/Nov/00
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.
%
% 30/Dec/2005 - added n-point geodesic computations, based on an algorithm
%               coded by Jeff Barton at Johns Hopkins APL in an m-file
%               I looked at at mathworks.com.


pi180=pi/180;
earth_radius=6378.137;

m=length(long)-1;

long1=reshape(long(1:end-1),m,1)*pi180;
long2=reshape(long(2:end)  ,m,1)*pi180;
lat1= reshape(lat(1:end-1) ,m,1)*pi180;
lat2= reshape(lat(2:end)   ,m,1)*pi180;

dlon = long2 - long1; 
dlat = lat2 - lat1; 
a = (sin(dlat/2)).^2 + cos(lat1) .* cos(lat2) .* (sin(dlon/2)).^2;
angles = 2 * atan2( sqrt(a), sqrt(1-a) );
dist = earth_radius * angles;


if nargin==3 && nargout>1   % Compute geodesics.

  % Cartesian unit vectors in rows of v1,v2
  v1=[cos(long1).*cos(lat1)   sin(long1).*cos(lat1)   sin(lat1) ];
  v2=[cos(long2).*cos(lat2)   sin(long2).*cos(lat2)   sin(lat2) ];

  % We want to get a unit vector tangent to the great circle.
  n1=cross(v1,v2,2); 
  t1=cross(n1,v1,2);
  t1=t1./repmat(sqrt(sum(t1.^2,2)),1,3);

  lons=zeros(m,N+1);
  lats=zeros(m,N+1);
  for k=1:m

   % Radials for all points
   p1=v1(k,:)'*cos(angles(k)*[0:N]/N) + t1(k,:)'*sin(angles(k)*[0:N]/N);

   lons(k,:)=atan2(p1(2,:),p1(1,:))/pi180;
   lats(k,:)=asin(p1(3,:))/pi180;

  end

end

