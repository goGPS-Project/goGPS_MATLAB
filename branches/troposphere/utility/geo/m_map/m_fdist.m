function [lon2,lat2,a21] = m_fdist(lon1,lat1,a12,s,spheroid)
% M_FDIST - On an ellipsoidal earth, compute location of
%         a point at a given bearing/distance, all using a vectorized
%         version of Vincenty's algorithm.
%
% [lon2,lat2,a21] = m_fdist(lon1,lat1,a12,s,spheroid)
%
% lat1 = GEODETIC latitude of first point (degrees)
% lon1 = longitude of first point (degrees)
% a12 = azimuth in degrees from first point to second point (forward)
% s = distance in meters  
% spheroid = (Optional) spheroid, defaults to 'wgs84'
%
% lat2, lon2 = second point (degrees)
% a21 = azimuth in degrees from second point to first point (backward)
%       (Azimuths are in degrees clockwise from north.)
%
% Inputs can be all the same size, or a mix of scalars and matrices.
% If a mixture, inputs are expanded as required.
%
%  Original algorithm source:
%  T. Vincenty, "Direct and Inverse Solutions of Geodesics on the Ellipsoid
%  with Application of Nested Equations", Survey Review, vol. 23, no. 176,
%  April 1975, pp 88-93.
%  Available at: http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
%
% See M_IDIST, M_GEODESIC

% Notes: (1) lat1,lon1,a12,s can be any (identical) size/shape. Outputs
%            will have the same size and shape.
%        (2) Vincenty describes his distance algorithm as precise to within
%            0.01 millimeters, subject to the ellipsoidal model.
%        (3) code written by heavily modifying M. Kleder's 'vdist' in line
%            with original Vincenty paper. m_fdist and m_idist are inverses to within
%            1e-10 degrees.
%       (4) No warranties; use at your own risk.

% R. Pawlowicz (rich@ocgy.ubc.ca) 9/Jan/2005
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.


pi180=pi/180;

MAP_ELLIP = struct ( 'normal', [1.0, 0], ...
    'sphere', [6370997.0, 0], ...
    'grs80' , [6378137.0, 1/298.257], ...
    'grs67' , [6378160.0, 1/247.247], ...
    'wgs84' , [6378137.0, 1/298.257223563], ... 
    'wgs72' , [6378135.0, 1/298.260], ...
    'wgs66' , [6378145.0, 1/298.250], ...
    'wgs60' , [6378165.0, 1/298.300], ...
    'clrk66', [6378206.4, 1/294.980], ...
    'clrk80', [6378249.1, 1/293.466], ...
    'intl24', [6378388.0, 1/297.000], ...
    'intl67', [6378157.5, 1/298.250]);

if nargin<5,
  spheroid='wgs84';
end;
ellip=getfield(MAP_ELLIP,spheroid); 
if length(ellip)~=2,
 disp(MAP_ELLIP);
 error('Spheroid not chosen from above list');
end;

% Do the equivalent of meshgrid-type calls to make all
% matrix sizes match. To do this we can have no more than
% two different numbers in any particular dimension, and one
% of these has to be a '1'.
allsize=[size(lon1);size(lat1);size(a12);size(s)];
i1=ones(1,size(allsize,2));
for k=1:size(allsize,2),
 rs=unique(allsize(:,k));
 if length(rs)==2 & rs(1)==1,
   j1=i1;j1(k)=rs(2);
   if allsize(1,k)==1,lon1=repmat(lon1,j1); end;
   if allsize(2,k)==1,lat1=repmat(lat1,j1); end;
   if allsize(3,k)==1,a12=repmat(a12,j1); end;
   if allsize(4,k)==1,s=repmat(s,j1); end;
 elseif length(rs)>2,
  error('incompatible array sizes!');  
 end;
end;

% reshape inputs
keepsize = size(lat1);
lat1=lat1(:);
lon1=lon1(:);
a12=a12(:);
s=s(:);
% Input check:
if any(abs(lat1)>90  )
    error('Input latitudes must be between -90 and 90 degrees, inclusive.')
end


 % correct for errors at exact poles by adjusting 0.6 millimeters:
kidx = abs(90-abs(lat1)) < 1e-10;
if any(kidx);
    lat1(kidx) = sign(lat1(kidx))*(90-(1e-10));
end

% Algorithm begins here

a=ellip(1);
b=a*(1-ellip(2));
f = (a-b)/a;
U1 = atan((1-f)*tan(lat1*pi180));
 
a12=a12*pi180;
sigma1=atan2( tan(U1),cos(a12));
alpha = asin( cos(U1).*sin(a12) );

u2 = cos(alpha).^2.*(a^2-b^2)/b^2;

A = 1+u2./16384.*(4096+u2.*(-768+u2.*(320-175.*u2)));
B = u2./1024.*(256+u2.*(-128+u2.*(74-47.*u2)));

sigmainit=s./(b*A);
sigma=sigmainit;
sigmaold=sigmainit;

itercount = 0;
k = logical(1+0*lat1);
deltasigma =0*lat1;
sigmam2=0*lat1;
cos2sigmam=0*lat1;
while any(k)  % force at least one execution
    itercount = itercount+1;
    if itercount > 50
        warning('More than 50 iterations');
        break
    end
    sigmam2(k)=2*sigma1(k)+sigma(k);
    cos2sigmam(k)=cos(sigmam2(k));
    deltasigma(k) = B(k).*sin(sigma(k)).*(cos2sigmam(k)+B(k)./4.*(cos(sigma(k)).*(-1+ ...
                   2.*cos2sigmam(k).^2)-B(k)./6.*cos2sigmam(k).*(-3+ ...
		   4.*sin(sigma(k)).^2).*(-3+4*cos2sigmam(k).^2)));
    sigmaold(k)=sigma(k);		 
    sigma(k)=sigmainit(k)+deltasigma(k);
    k = abs(sigma-sigmaold) > 1e-12;
end


numer=sin(sigma).*sin(a12);
denom=cos(U1).*cos(sigma) - sin(U1).*sin(sigma).*cos(a12);
lambda=atan2(numer,denom);
C = f/16*cos(alpha).^2.*(4+f*(4-3*cos(alpha).^2));
L=lambda-(1-C).*f.*sin(alpha).*(sigma+C.*sin(sigma).*(cos2sigmam+ ...
                           C.*cos(sigma).*(-1+2.*cos2sigmam.^2)));

lon2=mod(reshape((lon1+L/pi180),keepsize),360);

if nargout>1,
  numer=sin(U1).*cos(sigma)+cos(U1).*sin(sigma).*cos(a12);
  denom=(1-f)*sqrt(sin(alpha).^2 +  ...
                 (sin(U1).*sin(sigma) - cos(U1).*cos(sigma).*cos(a12)).^2 );
  lat2=reshape(atan2(numer,denom)/pi180,keepsize);
end;

if nargout>2,
  a21=mod(reshape( atan2( -sin(alpha),  ...
                       sin(U1).*sin(sigma)-cos(U1).*cos(sigma).*cos(a12) )/pi180, ...
                keepsize),360);
end;


