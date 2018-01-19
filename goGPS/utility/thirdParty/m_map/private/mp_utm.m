function  [X,Y,vals,labI]=mp_utm(optn,varargin)
% MP_UTM   Universal Transverse Mercator projection
%           This function should not be used directly; instead it is
%           is accessed by various high-level functions named M_*.

% mp_utm.m, Peter Lemmond (peter@whoi.edu)

% created mp_utm.m 13Aug98 from mp_tmerc.m, v1.2d distribution, by:
%
% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.
%
% Mathematical formulas for the projections and their inverses are taken from
%
%      Snyder, John P., Map Projections used by the US Geological Survey, 
%      Geol. Surv. Bull. 1532, 2nd Edition, USGPO, Washington D.C., 1983.
%

% 10/Dec/98 - PL added various ellipsoids.
% 17/May/12 - clarified hemisphere setting

global MAP_PROJECTION MAP_VAR_LIST

% define a structure of various ellipsoids. each has a name, and
% a vector consisting of equatorial radius and flattening. the first
% two are somewhat special cases.

MAP_ELLIP = struct ( 'normal', [1.0, 0], ...
    'sphere', [6370997.0, 0], ...
    'grs80' , [6378137.0, 1/298.257], ...
    'grs67' , [6378160.0, 1/247.247], ...
    'wgs84' , [6378137.0, 1/298.257], ...
    'wgs72' , [6378135.0, 1/298.260], ...
    'wgs66' , [6378145.0, 1/298.250], ...
    'wgs60' , [6378165.0, 1/298.300], ...
    'clrk66', [6378206.4, 1/294.980], ...
    'clrk80', [6378249.1, 1/293.466], ...
    'intl24', [6378388.0, 1/297.000], ...
    'intl67', [6378157.5, 1/298.250]);


name={'UTM'};

switch optn

  case 'name'

    X=name;

  case {'usage','set'}
    
    m_names=fieldnames(MAP_ELLIP);

    X=char({['     ''' varargin{1} ''''],...
	'     <,''lon<gitude>'',[min max]>',...
	'     <,''lat<itude>'',[min max]>',...
	'     <,''zon<e>'',value>',...
	'     <,''hem<isphere>'',[1|0] (0 for N)>',...
	'     <,''ell<ipsoid>'', one of',...
    reshape(sprintf('         %6s',m_names{:}),15,length(m_names))',...
        '                               >',...
	'     <,''rec<tbox>'', ( ''on'' | ''off'' )>'});

  case 'get'

    X=char([' Projection: ' MAP_PROJECTION.name '  (function: ' ...
	  MAP_PROJECTION.routine ')'],...
	[' longitudes: ' num2str(MAP_VAR_LIST.ulongs) ],...
	[' latitudes: ' num2str(MAP_VAR_LIST.ulats) ],...
	[' zone: ' num2str(MAP_VAR_LIST.zone) ],...
	[' hemisphere: ' num2str(MAP_VAR_LIST.hemisphere) ],...
	[' ellipsoid: ' MAP_VAR_LIST.ellipsoid ], ...
	[' Rectangular border: ' MAP_VAR_LIST.rectbox ]);


  case 'initialize'

    MAP_VAR_LIST=[];
    MAP_PROJECTION.name=varargin{1};
    MAP_VAR_LIST.ulongs = [-72 -68];
    MAP_VAR_LIST.ulats = [40 44];
    MAP_VAR_LIST.zone = 0;		% will be computed if not there
    MAP_VAR_LIST.hemisphere = -1;
    MAP_VAR_LIST.ellipsoid = 'normal';
    MAP_VAR_LIST.rectbox='off';
    k=2;

    while k<length(varargin) 
      switch varargin{k}(1:3)
	case 'lon'
	  MAP_VAR_LIST.ulongs=varargin{k+1};
	case 'lat'
	  MAP_VAR_LIST.ulats=varargin{k+1};
	case 'zon'
	  MAP_VAR_LIST.zone=varargin{k+1};
	case 'hem'
	  MAP_VAR_LIST.hemisphere=varargin{k+1};
        case 'ell'
	  MAP_VAR_LIST.ellipsoid=varargin{k+1};
        case 'rec'
          MAP_VAR_LIST.rectbox=varargin{k+1};
	otherwise
	  disp(['Unknown option: ' varargin{k}]);
      end
      k=k+2;
    end

    % set the zone and hemisphere if not specified
    
    if (~MAP_VAR_LIST.zone)
      MAP_VAR_LIST.zone = 1 + fix((mod(mean(MAP_VAR_LIST.ulongs)+180,360))/6);
    end
    
    if MAP_VAR_LIST.hemisphere==-1
      MAP_VAR_LIST.hemisphere=mean(MAP_VAR_LIST.ulats)<0;
    end

    % check for a valid ellipsoid. if not, use the normalized sphere
    
    if ~isfield(MAP_ELLIP,MAP_VAR_LIST.ellipsoid)
       MAP_VAR_LIST.ellipsoid = 'normal';
    end
    
    % Get X/Y and (if rectboxs are desired) recompute lat/long limits.

    mu_util('xylimits');
    if strcmp(MAP_VAR_LIST.rectbox,'on'), mu_util('lllimits'); end

  case 'll2xy'

    long=varargin{1};
    lat=varargin{2};
    vals=zeros(size(long));

    % Clip out-of-range values (lat/long boxes)
    
    if  strcmp(MAP_VAR_LIST.rectbox,'off') && ~strcmp(varargin{4},'off')
        vals=vals | long<=MAP_VAR_LIST.longs(1)+eps*10 | long>=MAP_VAR_LIST.longs(2)-eps*10 | ...
	              lat<=MAP_VAR_LIST.lats(1)+eps*10 |   lat>=MAP_VAR_LIST.lats(2)-eps*10;
        [long,lat]=mu_util('clip',varargin{4},long,MAP_VAR_LIST.longs(1),long<MAP_VAR_LIST.longs(1),lat);
        [long,lat]=mu_util('clip',varargin{4},long,MAP_VAR_LIST.longs(2),long>MAP_VAR_LIST.longs(2),lat);
        [lat,long]=mu_util('clip',varargin{4},lat,MAP_VAR_LIST.lats(1),lat<MAP_VAR_LIST.lats(1),long);
        [lat,long]=mu_util('clip',varargin{4},lat,MAP_VAR_LIST.lats(2),lat>MAP_VAR_LIST.lats(2),long);
    end

    % do the forward transformation

    [X,Y] = mu_ll2utm(lat,long,MAP_VAR_LIST.zone,MAP_VAR_LIST.hemisphere, ...
	getfield(MAP_ELLIP,MAP_VAR_LIST.ellipsoid));
    
    % Clip out-of-range values (rectboxes)

    if strcmp(MAP_VAR_LIST.rectbox,'on')  && ~strcmp(varargin{4},'off')
        vals= vals | X<=MAP_VAR_LIST.xlims(1)+eps*10 | X>=MAP_VAR_LIST.xlims(2)-eps*10 | ...
                     Y<=MAP_VAR_LIST.ylims(1)+eps*10 | Y>=MAP_VAR_LIST.ylims(2)-eps*10;
        [X,Y]=mu_util('clip',varargin{4},X,MAP_VAR_LIST.xlims(1),X<MAP_VAR_LIST.xlims(1),Y);
        [X,Y]=mu_util('clip',varargin{4},X,MAP_VAR_LIST.xlims(2),X>MAP_VAR_LIST.xlims(2),Y);
        [Y,X]=mu_util('clip',varargin{4},Y,MAP_VAR_LIST.ylims(1),Y<MAP_VAR_LIST.ylims(1),X);
        [Y,X]=mu_util('clip',varargin{4},Y,MAP_VAR_LIST.ylims(2),Y>MAP_VAR_LIST.ylims(2),X);
    end

  case 'xy2ll'

    [Y,X] = mu_utm2ll(varargin{1}, varargin{2}, MAP_VAR_LIST.zone, ...
	MAP_VAR_LIST.hemisphere, getfield(MAP_ELLIP,MAP_VAR_LIST.ellipsoid));
        
  case 'xgrid'
   
    [X,Y,vals,labI]=mu_util('xgrid',MAP_VAR_LIST.longs,MAP_VAR_LIST.lats,varargin{1},31,varargin{2:3});

  case 'ygrid'
   
    [X,Y,vals,labI]=mu_util('ygrid',MAP_VAR_LIST.lats,MAP_VAR_LIST.longs,varargin{1},31,varargin{2:3});

  case 'box'

    [X,Y]=mu_util('box',31);

end


%-------------------------------------------------------------------

function [x,y] = mu_ll2utm (lat,lon, zone, hemisphere,ellipsoid)
%mu_ll2utm		Convert geodetic lat,lon to X/Y UTM coordinates
%
%	[x,y] = mu_ll2utm (lat, lon, zone, hemisphere,ellipsoid)
%
%	input is latitude and longitude vectors, zone number, 
%		hemisphere(N=0,S=1), ellipsoid info [eq-rad, flat]
%	output is X/Y vectors
%
%	see also	mu_utm2ll, utmzone


% some general constants

DEG2RADS    = 0.01745329252;
RADIUS      = ellipsoid(1);
FLAT        = ellipsoid(2);
K_NOT       = 0.9996;
FALSE_EAST  = 500000;
FALSE_NORTH = 10000000;

% check for valid numbers

if (max(abs(lat)) > 90)
  error('latitude values exceed 89 degree');
  return;
end

if ((zone < 1) || (zone > 60))
  error ('utm zones only valid from 1 to 60');
  return;
end

% compute some geodetic parameters

lambda_not  = ((-180 + zone*6) - 3) * DEG2RADS;

e2  = 2*FLAT - FLAT*FLAT;
e4  = e2 * e2;
e6  = e4 * e2;
ep2 = e2/(1-e2);

% some other constants, vectors

lat = lat * DEG2RADS;
lon = lon * DEG2RADS;

sinL = sin(lat);
tanL = tan(lat);
cosL = cos(lat);

T = tanL.*tanL;
C = ep2 * (cosL.*cosL);
A = (lon - lambda_not).*cosL;
A2 = A.*A;
A4 = A2.*A2;
S = sinL.*sinL;

% solve for N

N = RADIUS ./ (sqrt (1-e2*S));

% solve for M

M0 = 1 - e2*0.25 - e4*0.046875 - e6*0.01953125;
M1 = e2*0.375 + e4*0.09375 + e6*0.043945313;
M2 = e4*0.05859375 + e6*0.043945313;
M3 = e6*0.011393229;
M = RADIUS.*(M0.*lat - M1.*sin(2*lat) + M2.*sin(4*lat) - M3.*sin(6*lat));

% solve for x

X0 = A4.*A/120;
X1 = 5 - 18*T + T.*T + 72*C - 58*ep2;
X2 = A2.*A/6;
X3 = 1 - T + C;
x = N.*(A + X3.*X2 + X1.* X0);

% solve for y

Y0 = 61 - 58*T + T.*T + 600*C - 330*ep2;
Y1 = 5 - T + 9*C + 4*C.*C;

y = M + N.*tanL.*(A2/2 + Y1.*A4/24 + Y0.*A4.*A2/720);


% finally, do the scaling and false thing. if using a unit-normal radius,
% we don't bother.

x = x*K_NOT + (RADIUS>1) * FALSE_EAST;

y = y*K_NOT;
if (hemisphere)
  y = y + (RADIUS>1) * FALSE_NORTH;
end

return



%-------------------------------------------------------------------

function [lat,lon] = mu_utm2ll (x,y, zone, hemisphere,ellipsoid)
%mu_utm2ll		Convert X/Y UTM coordinates to geodetic lat,lon 
%
%	[lat,lon] = mu_utm2ll (x,y, zone, hemisphere,ellipsoid)
%
%	input is X/Y vectors, zone number, hemisphere(N=0,S=1),
%		ellipsoid info [eq-rad, flat]
%	output is lat/lon vectors
%
%	see also	mu_ll2utm, utmzone


% some general constants

DEG2RADS    = 0.01745329252;
RADIUS      = ellipsoid(1);
FLAT        = ellipsoid(2);
K_NOT       = 0.9996;
FALSE_EAST  = 500000;
FALSE_NORTH = 10000000;

if ((zone < 1) || (zone > 60))
  error ('utm zones only valid from 1 to 60');
  return;
end

% compute some geodetic parameters

e2  = 2*FLAT - FLAT*FLAT;
e4  = e2 * e2;
e6  = e4 * e2;
eps = e2 / (1-e2);
em1 = sqrt(1-e2);
e1  = (1-em1)/(1+em1);
e12 = e1*e1;

lambda_not  = ((-180 + zone*6) - 3) * DEG2RADS;

% remove the false things

x = x - (RADIUS>1)*FALSE_EAST;
if (hemisphere)
  y = y - (RADIUS>1)*FALSE_NORTH;
end

% compute the footpoint latitude

M = y/K_NOT;
mu = M/(RADIUS * (1 - 0.25*e2 - 0.046875*e4 - 0.01953125*e6));
foot = mu + (1.5*e1 - 0.84375*e12*e1)*sin(2*mu) ...
    + (1.3125*e12 - 1.71875*e12*e12)*sin(4*mu) ...
    + (1.57291666667*e12*e1)*sin(6*mu) ...
    + (2.142578125*e12*e12)*sin(8*mu);

% some other terms

sinF = sin(foot);
cosF = cos(foot);
tanF = tan(foot);

N = RADIUS ./ sqrt(1-e2*(sinF.*sinF));
T = tanF.*tanF;
T2 = T.*T;
C = eps * cosF.*cosF;
C2 = C.*C;
denom = sqrt(1-e2*(sinF.*sinF));
R = RADIUS * em1*em1 ./ (denom.*denom.*denom);
D = x./(N*K_NOT);
D2 = D.*D;
D4 = D2.*D2;

% can now compute the lat and lon

lat = foot - (N.*tanF./R) .* (0.5*D2 - (5 + 3*T + 10*C - 4*C2 - 9*eps).*D4/24 ...
    + (61 + 90*T + 298*C + 45*T2 - 252*eps - 3*C2) .* D4 .* D2/720);

lon = lambda_not + (D - (1 + 2*T +C).*D2.*D/6 + ...
    (5 - 2*C + 28*T - 3*C2 + 8*eps + 24*T2).*D4.*D./120)./cosF;


% convert back to degrees;

lat=lat/DEG2RADS;
lon=lon/DEG2RADS;

return
