function [dpsi, deps] = nut2000a (date1, date2)

% nutation based on iau 2000a theory

% input

%  date1, date2 = tt julian date 
%  (julian date = date1 + date2)

% output

%  dpsi = nutation in longitude in radians

%  deps = nutation in obliquity in radians

% reference

%  Nutation Series Evaluation in NOVAS 3.0
%  USNO Circular No. 181, December 15, 2009

% ported from NOVAS 3.0 Fortran subroutine

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global inutate nals napl icpl cls

if (inutate == 1)
    
    % read data files
    
    nals = csvread('nals.csv');
    
    napl = csvread('napl.csv');
    
    icpl = csvread('icpl.csv');
    
    cls = csvread('cls.csv');
    
    % transpose matrices
    
    nals = nals';
    
    napl = napl';
    
    icpl = icpl';
    
    cls = cls';
    
    % reset flag
    
    inutate = 0;
    
end

% arc seconds to radians

das2r = 4.848136811095359935899141d-6;

% arc seconds in a full circle

turnas = 1296000.0d0;

% 2 * pi

d2pi = 6.283185307179586476925287d0;

% units of 0.1 microarcsecond to radians

u2r = das2r / 1.0d7;

% reference epoch (j2000)

dj0 = 2451545.0d0;

% days per julian century

djc = 36525.0d0;

% -------------------------
% luni-solar nutation model
% -------------------------

%  ---------------
%  planetary terms
%  ---------------

% interval between fundamental epoch j2000.0 and given date (jc)

t = ((date1 - dj0) + date2) / djc;

% -------------------
% luni-solar nutation
% -------------------

% fundamental (delaunay) arguments from simon et al. (1994)

% mean anomaly of the moon (radians)

el = mod (485868.249036d0 + t * (1717915923.2178d0 ...
    + t * (31.8792d0 + t * (0.051635d0 ...
    + t * (-0.00024470d0)))), turnas) * das2r;

% mean anomaly of the sun (radians)

elp = mod (1287104.79305d0 + t * (129596581.0481d0 ...
    + t * (-0.5532d0 + t * (0.000136d0 ...
    + t * (-0.00001149d0)))), turnas) * das2r;

% mean argument of the latitude of the moon (radians)

f = mod (335779.526232d0 + t * (1739527262.8478d0 ...
    + t * (-12.7512d0 + t * (-0.001037d0 ...
    + t * (0.00000417d0)))), turnas) * das2r;

% mean elongation of the moon from the sun (radians)

d = mod (1072260.70369d0 + t * (1602961601.2090d0 ...
    + t * (-6.3706d0 + t * (0.006593d0 ...
    + t * (-0.00003169d0)))), turnas) * das2r;

% mean longitude of the ascending node of the moon (radians)

om = mod (450160.398036d0 + t * (-6962890.5431d0 ...
    + t * (7.4722d0 + t * (0.007702d0 ...
    + t * (-0.00005939d0)))), turnas) * das2r;

% initialize the nutation values

dp = 0.0d0;

de = 0.0d0;

% summation of luni-solar nutation series (in reverse order)

for i = 678: -1: 1
    
    arg = mod ((nals(1, i)) * el  + (nals(2, i)) * elp + (nals(3, i)) * f ...
        + (nals(4, i)) * d  + (nals(5, i)) * om, d2pi);
    
    sarg = sin(arg);
    
    carg = cos(arg);
    
    dp = dp + (cls(1,i) + cls(2,i) * t) * sarg + cls(3,i) * carg;
    
    de = de + (cls(4,i) + cls(5,i) * t) * carg + cls(6,i) * sarg;
    
end

% convert from 0.1 microarcsec units to radians

dpsils = dp * u2r;

depsls = de * u2r;

% ------------------
% planetary nutation
% ------------------

% mean anomaly of the moon (radians)

al = mod (2.35555598d0 + 8328.6914269554d0 * t, d2pi);

% mean anomaly of the sun (radians)

alsu = mod (6.24006013d0 + 628.301955d0 * t, d2pi);

% mean argument of the latitude of the moon (radians)

af = mod (1.627905234d0 + 8433.466158131d0 * t, d2pi);

% mean elongation of the moon from the sun (radians)

ad = mod (5.198466741d0 + 7771.3771468121d0 * t, d2pi);

% mean longitude of the ascending node of the moon (radians)

aom = mod (2.18243920d0 - 33.757045d0 * t, d2pi);

% general accumulated precession in longitude (radians)

apa = (0.02438175d0 + 0.00000538691d0 * t) * t;

% planetary longitudes, mercury through neptune (souchay et al. 1999)

alme = mod (4.402608842d0 + 2608.7903141574d0 * t, d2pi);

alve = mod (3.176146697d0 + 1021.3285546211d0 * t, d2pi);

alea = mod (1.753470314d0 +  628.3075849991d0 * t, d2pi);

alma = mod (6.203480913d0 +  334.0612426700d0 * t, d2pi);

alju = mod (0.599546497d0 +   52.9690962641d0 * t, d2pi);

alsa = mod (0.874016757d0 +   21.3299104960d0 * t, d2pi);

alur = mod (5.481293871d0 +    7.4781598567d0 * t, d2pi);

alne = mod (5.321159000d0 +    3.8127774000d0 * t, d2pi);

% initialize the nutation values

dp = 0.0d0;

de = 0.0d0;

% summation of planetary nutation series (in reverse order)

for i = 687: -1: 1
    
    arg = mod ((napl(1, i)) * al + (napl(2, i)) * alsu + (napl(3, i)) * af ...
        + (napl(4, i)) * ad + (napl(5, i)) * aom  + (napl(6, i)) * alme ...
        + (napl(7, i)) * alve + (napl(8, i)) * alea + (napl(9, i)) * alma ...
        + (napl(10, i)) * alju + (napl(11, i)) * alsa + (napl(12, i)) * alur ...
        + (napl(13, i)) * alne + (napl(14, i)) * apa, d2pi);
    
    sarg = sin(arg);
    
    carg = cos(arg);
    
    dp = dp + (icpl(1,i)) * sarg + (icpl(2,i)) * carg;
    
    de = de + (icpl(3,i)) * sarg + (icpl(4,i)) * carg;
    
end

% convert from 0.1 microarcsec units to radians

dpsipl = dp * u2r;

depspl = de * u2r;

% add planetary and luni-solar components

dpsi = dpsipl + dpsils;

deps = depspl + depsls;

end
