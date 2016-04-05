function [ra2, dec2, pmra2, pmdec2, parx2, rv2] = catran ...
    (it, date1, ra1, dec1, pmra1, pmdec1, parx1, rv1, date2)

% this function transforms a star's catalog quantities for
% a change of epoch and/or equator and equinox.  it can also be
% used to rotate catalog quantities on the dynamical equator and
% equinox of j2000.0 to the icrs or vice versa.

% input

%  it     = transformation option 

%           set it=1 to change epoch (same equator and equinox)
%           set it=2 to change equator and equinox (same epoch)
%           set it=3 to change equator and equinox and epoch
%           set it=4 to change equator and equinox of j2000.0 to icrs
%           set it=5 to change icrs to equator and equinox of j2000.0

%  date1  = tt julian date, or year, of original catalog
%           data (the following six arguments)

%  ra1    = original mean right ascension in hours

%  dec1   = original mean declination in degrees

%  pmra1  = original proper motion in ra in milliarcseconds/year

%  pmdec1 = original proper motion in dec in milliarcseconds/year

%  parx1  = original parallax in milliarcseconds 

%  rv1    = original radial velocity in kilometers/second
%               
%  date2  = tt julian date, or year, for transformed
%           output data (the following six arguments) 

% output

%  ra2    = transformed mean right ascension in hours 

%  dec2   = transformed mean declination in degrees 

%  pmra2  = transformed proper motion in ra in milliarcseconds/year

%  pmdec2 = transformed proper motion in dec in milliarcseconds/year 

%  parx2  = transformed parallax in milliarcseconds 

%  rv2    = transformed radial velocity in kilometers/second               

% note 1:  date1 and date2 may be specified either as a julian
% date (e.g., 2433282.5d0) or a julian year and fraction
% (e.g., 1950.0d0).  values less than 10000 are assumed to
% be years.  for it=2 or it=3, either date1 or date2 must be
% 2451545.0 or 2000.0 (j2000.0).  for it=4 and it=5, date1 and
% date2 are ignored.

% note 2:  it=1 updates the star's data to account for
% the star's space motion between the first and second dates,
% within a fixed reference system.  it=2 applies a rotation
% of the reference system corresponding to precession between
% the first and second dates, but leaves the star fixed in space.
% it=3 provides both transformations.  it=4 and it=5 provide a
% a fixed rotation about very small angles (<0.1 arcsecond) to
% take data from the dynamical system of j2000.0 to the icrs (it=4)
% or vice versa (it=5).

% note 3:  for it=1, input data can be in any fixed reference
% system. for it=2 or it=3, this subroutine assumes the input data
% is in the dynamical system and produces output in the dynamical
% system.  for it=4, the input data must be on the dynamical equator
% and equinox of j2000.0.  for it=5, the input data must be in the
% icrs.

% note 4:  this subroutine cannot be properly used to bring data
% from old star catalogs into the modern system, because
% old catalogs were compiled using a set of constants that are
% incompatible with modern values.  in particular, it should not
% be used for catalogs whose positions and proper motions were
% derived by assuming a precession constant significantly different
% from the value implicit in subroutine preces.

% ported from NOVAS 3.0

%%%%%%%%%%%%%%%%%%%%%%%

seccon = 180.0d0 * 3600.0d0 / pi;

% speed of light in kilometers/second

c = 1.0d-3 * 2997924580.0d0;

% au in kilometers

aukm = 1.0d-3 * 499.0047838061d0 * c;

% --- if necessary, compute julian dates

% subroutine uses tdb julian dates internally, but no
% distinction between tdb and tt is necessary

if (date1 < 10000.0d0)
    
    tjd1 = 2451545.0d0 + (date1 - 2000.0d0) * 365.25d0;
    
else
    tjd1 = date1;
end

if (date2 < 10000.0d0)
    
    tjd2 = 2451545.0d0 + (date2 - 2000.0d0) * 365.25d0;
    
else
    tjd2 = date2;
end

% --- convert input angular components to vectors

% if parallax is unknown, undetermined, or zero, set it to 1e-6
% milliarcsecond, corresponding to a distance of 1 gigaparsec

paralx = parx1;

if (paralx <= 0.0d0)
    
    paralx = 1.0d-6;
    
end

% convert right ascension, declination, and parallax to position
% vector in equatorial system with units of au

dist = 1.0d0 / sin(paralx * 1.0d-3 / seccon);

r = ra1 * 54000.0d0 / seccon;

d = dec1 * 3600.0d0 / seccon;

cra = cos(r);

sra = sin(r);

cdc = cos(d);

sdc = sin(d);

pos1(1) = dist * cdc * cra;

pos1(2) = dist * cdc * sra;

pos1(3) = dist * sdc;

% compute doppler factor, which accounts for change in
% light travel time to star

k = 1.0d0 / (1.0d0 - rv1 / c);

% convert proper motion and radial velocity to orthogonal
% components of motion, in spherical polar system at star's
% original position, with units of au/day

pmr = pmra1  / (paralx * 365.25d0) * k;

pmd = pmdec1 / (paralx * 365.25d0) * k;

rvl = rv1 * 86400.0d0 / aukm * k;

% transform motion vector to equatorial system

vel1(1) = - pmr * sra - pmd * sdc * cra + rvl * cdc * cra;

vel1(2) =   pmr * cra - pmd * sdc * sra + rvl * cdc * sra;

vel1(3) =               pmd * cdc       + rvl * sdc;

% --- update star's position vector for space motion

% (only if it = 1 or it = 3)

if (it == 1 || it == 3)

    for j = 1:1:3
        pos2(j) = pos1(j) + vel1(j) * (tjd2 - tjd1);
        
        vel2(j) = vel1(j);
    end
    
else

    for j = 1:1:3
        pos2(j) = pos1(j);
        
        vel2(j) = vel1(j);
    end
    
end

% --- precess position and velocity vectors (only if it = 2 or it = 3)

if (it == 2 || it == 3)

    for j = 1:1:3
        pos1(j) = pos2(j);
        
        vel1(j) = vel2(j);
    end

    pos2 = preces (tjd1, pos1, tjd2);
    
    vel2 = preces (tjd1, vel1, tjd2);
end

% --- rotate dynamical j2000.0 position and velocity vectors to icrs
% (only if it = 4)

if (it == 4)
    pos2 = frame (pos1, -1);
    
    vel2 = frame (vel1, -1);
end

% --- rotate icrs position and velocity vectors to dynamical j2000.0
% (only if it = 5)

if (it == 5)
    pos2 = frame (pos1, 1);
    
    vel2 = frame (vel1, 1);
end

% --- convert vectors back to angular components for output

% from updated position vector, obtain star's new position
% expressed as angular quantities

xyproj = sqrt(pos2(1)^2 + pos2(2)^2);

r = 0.0d0;

if (xyproj > 0.0d0)
    r = atan2 (pos2(2), pos2(1));
end

ra2 = r * seccon / 54000.0d0;

if (ra2 <  0.0d0)
    ra2 = ra2 + 24.0d0;
end

if (ra2 >= 24.0d0)
    ra2 = ra2 - 24.0d0;
end

d = atan2 (pos2(3), xyproj);

dec2 = d * seccon / 3600.0d0;

dist = sqrt (pos2(1)^2 + pos2(2)^2 + pos2(3)^2);

paralx = asin (1.0d0 / dist) * seccon * 1.0d3;

parx2 = paralx;

% transform motion vector back to spherical polar system at star's
% new position

cra = cos(r);

sra = sin(r);

cdc = cos(d);

sdc = sin(d);

pmr = - vel2(1) * sra       + vel2(2) * cra;

pmd = - vel2(1) * cra * sdc - vel2(2) * sra * sdc + vel2(3) * cdc;

rvl =   vel2(1) * cra * cdc + vel2(2) * sra * cdc + vel2(3) * sdc;

% convert components of motion from au/day to normal catalog units

pmra2  = pmr * paralx * 365.25d0 / k;

pmdec2 = pmd * paralx * 365.25d0 / k;

rv2    = rvl * (aukm / 86400.d0) / k;

% take care of zero-parallax case

if (parx2 <= 1.01d-6)

    parx2 = 0.0d0;

    rv2 = rv1;
end

