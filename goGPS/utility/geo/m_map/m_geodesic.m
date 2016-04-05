function [s,lon,lat] = m_geodesic(lon1,lat1,lon2,lat2,Npoints,varargin)
% M_GEODESIC - Compute points along geodesics of an ellipsoidal earth.
%
%       [S,LON,LAT] = m_geodesic(lon1,lat1,lon2,lat2,Npoints,spheroid)
%
% lat1 = GEODETIC latitude of first point (degrees)
% lon1 = longitude of first point (degrees)
% lat2, lon2 = second point (degrees)
% Npoints = number of points along geodesic
% spheroid = (Optional) spheroid, defaults to 'wgs84'
% S = distance in meters  
% LON,LAT = points on geodesics.
%
% Note that inputs can be the same size, or a mixture of scalars and matrices.
% However, output curves will be on columns of a matrix (i.e. form of
% inputs will not be preserved).
%
% See M_FDIST, M_IDIST, M_LLDIST

% R. pawlowicz (rich@ocgy.ubc.ca) 10/Jan/2005
%
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.



% Do the equivalent of meshgrid-type calls to make all
% matrix sizes match. To do this we can have no more than
% two different numbers in any particular dimension, and one
% of these has to be a '1'.
allsize=[size(lon1);size(lat1);size(lon2);size(lat2)];
i1=ones(1,size(allsize,2));
for k=1:size(allsize,2),
 rs=unique(allsize(:,k));
 if length(rs)==2 & rs(1)==1,
   j1=i1;j1(k)=rs(2);
   if allsize(1,k)==1,lon1=repmat(lon1,j1); end;
   if allsize(2,k)==1,lat1=repmat(lat1,j1); end;
   if allsize(3,k)==1,lon2=repmat(lon2,j1); end;
   if allsize(4,k)==1,lat2=repmat(lat2,j1); end;
 elseif length(rs)>2,
  error('incompatible array sizes!');  
 end;
end;

% Get the distances and bearings

[S,A12,A21]=m_idist(lon1,lat1,lon2,lat2,varargin{:});

s=[0:Npoints-1]'/(Npoints-1)*S(:)';

[lon,lat]=m_fdist(lon1(:)',lat1(:)',A12(:)',s,varargin{:});



