function [s,a12,a21] = m_idist(lon1,lat1,lon2,lat2,spheroid)
% M_IDIST - On an ellipsoidal earth, compute the distance between
%         two points within a few millimeters of accuracy, compute forward
%         azimuth, and compute backward azimuth, all using a vectorized
%         version of Vincenty's algorithm.
%
%       [s,a12,a21] = m_idist(lon1,lat1,lon2,lat2,spheroid)
%
% lat1 = GEODETIC latitude of first point (degrees)
% lon1 = longitude of first point (degrees)
% lat2, lon2 = second point (degrees)
% spheroid = (Optional) spheroid, defaults to 'wgs84'
% s = distance in meters 
% a12 = azimuth in degrees from first point to second point (forward)
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
%  Code slightly modified version of vdist.m (M. Kleder) available
%  at user-contributed area of www.mathworks.com
%
% See also M_FDIST, M_GEODESIC, M_LLDIST

% (Original notes)
% Notes: (1) lat1,lon1,lat2,lon2 can be any (identical) size/shape. Outputs
%            will have the same size and shape.
%        (2) Error correcting code, convergence failure traps, antipodal
%            corrections, polar error corrections, WGS84 ellipsoid
%            parameters, testing, and comments: Michael Kleder, 2004.
%        (3) Azimuth implementation (including quadrant abiguity
%            resolution) and code vectorization, Michael Kleder, Sep 2005.
%        (4) Vectorization is convergence sensitive; that is, quantities
%            which have already converged to within tolerance are not
%            recomputed during subsequent iterations (while other
%            quantities are still converging).
%        (5) Vincenty describes his distance algorithm as precise to within
%            0.01 millimeters, subject to the ellipsoidal model.
%        (6) For distance calculations, essentially antipodal points are
%            treated as exactly antipodal, potentially reducing accuracy
%            slightly.
%        (7) Distance failures for points exactly at the poles are
%            eliminated by moving the points by 0.6 millimeters.
%        (8) The Vincenty distance algorithm was transcribed verbatim by
%            Peter Cederholm, August 12, 2003. It was modified and
%            translated to English by Michael Kleder.
%            Mr. Cederholm's website is http://www.plan.aau.dk/~pce/
%        (9) Distances agree with the Mapping Toolbox, version 2.2 (R14SP3)
%            with a max relative difference of about 5e-9, except when the
%            two points are nearly antipodal, and except when one point is
%            near the equator and the two longitudes are nearly 180 degrees
%            apart. This function (vdist) is more accurate in such cases.
%            For example, note this difference (as of this writing):
%            >>vdist(0.2,305,15,125)
%            18322827.0131551
%            >>distance(0.2,305,15,125,[6378137 0.08181919])
%            0
%       (10) Azimuths FROM the north pole (either forward starting at the
%            north pole or backward when ending at the north pole) are set
%            to 180 degrees by convention. Azimuths FROM the south pole are
%            set to 0 degrees by convention.
%       (11) Azimuths agree with the Mapping Toolbox, version 2.2 (R14SP3)
%            to within about a hundred-thousandth of a degree, except when
%            traversing to or from a pole, where the convention for this
%            function is described in (10), and except in the cases noted
%            above in (9).
%       (12) No warranties; use at your own risk.
%
%  R. Pawlowicz (rich@ocgy.ubc.ca)
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.
%
%   9/jan/2005  - changed name, altered inputs to m_map style (lon first),
%                added ellipses, many minor simplifications to tighten code.
%   9/apr/2009  - distances was NaN if start and end points are the same,
%                 and I have fixed this.
%   12/May/2014 - my fix didn't work for arrays. Mike Casey suggested a fix for
%                 this fix.

pi180=pi/180;

MAP_ELLIP = struct ( 'normal', [1.0, 0], ...
    'sphere', [6370997.0, 0], ...
    'grs80' , [6378137.0, 1/298.257], ...
    'grs67' , [6378160.0, 1/247.247], ...
    'wgs84' , [6378137.0, 1/298.257223563],  ...
    'wgs72' , [6378135.0, 1/298.260], ...
    'wgs66' , [6378145.0, 1/298.250], ...
    'wgs60' , [6378165.0, 1/298.300], ...
    'clrk66', [6378206.4, 1/294.980], ...
    'clrk80', [6378249.1, 1/293.466], ...
    'intl24', [6378388.0, 1/297.000], ...
    'intl67', [6378157.5, 1/298.250]);


if nargin<5
  spheroid='wgs84';
end
ellip=getfield(MAP_ELLIP,spheroid); 
if length(ellip)~=2
 disp(MAP_ELLIP);
 error('Spheroid not chosen from above list');
end


% Do the equivalent of meshgrid-type calls to make all
% matrix sizes match. To do this we can have no more than
% two different numbers in any particular dimension, and one
% of these has to be a '1'.
allsize=[size(lon1);size(lat1);size(lon2);size(lat2)];
i1=ones(1,size(allsize,2));
for k=1:size(allsize,2)
 rs=unique(allsize(:,k));
 if length(rs)==2 && rs(1)==1
   j1=i1;j1(k)=rs(2);
   if allsize(1,k)==1,lon1=repmat(lon1,j1); end
   if allsize(2,k)==1,lat1=repmat(lat1,j1); end
   if allsize(3,k)==1,lon2=repmat(lon2,j1); end
   if allsize(4,k)==1,lat2=repmat(lat2,j1); end
 elseif length(rs)>2
  error('incompatible array sizes!');  
 end
end
 
% reshape inputs
keepsize = size(lat1);
lat1=lat1(:);
lon1=lon1(:);
lat2=lat2(:);
lon2=lon2(:);
% Input check:
if any(abs(lat1)>90 | abs(lat2)>90)
    error('Input latitudes must be between -90 and 90 degrees, inclusive.')
end




% correct for errors at exact poles by adjusting 0.6 millimeters:
kidx = abs(90-abs(lat1)) < 1e-10;
if any(kidx)
    lat1(kidx) = sign(lat1(kidx))*(90-(1e-10));
end
kidx = abs(90-abs(lat2)) < 1e-10;
if any(kidx)
    lat2(kidx) = sign(lat2(kidx))*(90-(1e-10));
end

% Algorithm begins here

a=ellip(1);
b=a*(1-ellip(2));
f = (a-b)/a;
U1 = atan((1-f)*tan(lat1*pi180));
U2 = atan((1-f)*tan(lat2*pi180));

% Get longitude difference (short way around)
L = abs(mod(lon2,360)-mod(lon1,360))*pi180;
kidx = L > pi;
if any(kidx)
    L(kidx) = 2*pi - L(kidx);
end

lambda = L;
lambdaold = zeros(size(lat1));
alpha = lambdaold;
sinsigma=lambdaold;
cossigma=lambdaold;
sigma = lambdaold;
cos2sigmam = lambdaold;
C = lambdaold;
k = true(size(lat1));
itercount = 0;
warninggiven = false;
while any(k)  % force at least one execution
    itercount = itercount+1;
    if itercount > 50
        if ~warninggiven
            warning(['Essentially antipodal points encountered. ' ...
                'Precision may be reduced slightly.']);
        end
        lambda(k) = pi;
        break
    end
    % eliminate rare imaginary portions at limit of numerical precision:
    sinsigma(k) = real(  sqrt((cos(U2(k)).*sin(lambda(k))).^2+ ...
                       (cos(U1(k)).*sin(U2(k))-sin(U1(k)).*cos(U2(k)).*cos(lambda(k))).^2) );
    cossigma(k) = real(  sin(U1(k)).*sin(U2(k))+cos(U1(k)).*cos(U2(k)).*cos(lambda(k))  );
    sigma(k) = atan2(sinsigma(k),cossigma(k));
    alpha(k) = asin(cos(U1(k)).*cos(U2(k)).*sin(lambda(k))./sinsigma(k));
%%    if isnan(alpha(k)), alpha(k)=0; end; % this occurs when points are colocated (RP 9/apr/09)
    alpha(isnan(alpha(k)))=0;     % Fix for above line - Thanks to M. Casey
    cos2sigmam(k) = cossigma(k)-2*sin(U1(k)).*sin(U2(k))./cos(alpha(k)).^2;
    C(k) = f/16*cos(alpha(k)).^2.*(4+f*(4-3*cos(alpha(k)).^2));
    lambdaold(k) = lambda(k);
    lambda(k) = L(k)+(1-C(k)).*f.*sin(alpha(k)).*(sigma(k)+ ...
                C(k).*sinsigma(k).*(cos2sigmam(k)+ ...
		C(k).*cossigma(k).*(-1+2.*cos2sigmam(k).^2)));
 
    % correct for convergence failure in the case of essentially antipodal
    % points
    if any(lambda(k) > pi)
        warning(['Essentially antipodal points encountered. ' ...
            'Precision may be reduced slightly.']);
        warninggiven = true;
        lambdaold(lambda>pi) = pi;
        lambda(lambda>pi) = pi;
    end
    k = abs(lambda-lambdaold) > 1e-12;
end

u2 = cos(alpha).^2.*(a^2-b^2)/b^2;
A = 1+u2./16384.*(4096+u2.*(-768+u2.*(320-175.*u2)));
B = u2./1024.*(256+u2.*(-128+u2.*(74-47.*u2)));
deltasigma = B.*sin(sigma).*(cos2sigmam+B./4.*(cos(sigma).*(-1+2.*cos2sigmam.^2) ...
                -B./6.*cos2sigmam.*(-3+4.*sin(sigma).^2).*(-3+4*cos2sigmam.^2)));
s = reshape(b.*A.*(sigma-deltasigma),keepsize);

if nargout > 1
    % From point #1 to point #2
    % correct sign of lambda for azimuth calcs:
    lambda = abs(lambda);
    kidx=sign(sin((lon2-lon1)*pi180)) .* sign(sin(lambda)) < 0;
    lambda(kidx) = -lambda(kidx);
    numer = cos(U2).*sin(lambda);
    denom = cos(U1).*sin(U2)-sin(U1).*cos(U2).*cos(lambda);
    a12 = reshape( mod( atan2(numer,denom)/pi180 , 360 ), keepsize);
    % from poles:
    a12(lat1 <= -90) = 0;
    a12(lat1 >= 90 ) = 180;
end

if nargout > 2
    % From point #2 to point #1
    % correct sign of lambda for azimuth calcs:
    lambda = abs(lambda);
    kidx=sign(sin((lon1-lon2)*pi180)) .* sign(sin(lambda)) < 0;
    lambda(kidx)=-lambda(kidx);
    numer = cos(U1).*sin(lambda);
    denom = sin(U1).*cos(U2)-cos(U1).*sin(U2).*cos(lambda);
    a21 = reshape( mod(atan2(numer,denom)/pi180, 360 ), keepsize);
    % backwards from poles:
    a21(lat2 >= 90) = pi;
    a21(lat2 <= -90) = 0;
end
return

