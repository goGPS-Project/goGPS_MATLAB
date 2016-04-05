function skypos = place (tjd, object, locatn, icoord, star, observ)

% this function computes the apparent direction of a star or solar
% system body at a specified time and in a specified coordinate
% system.  based on kaplan, et al. (1989), astronomical journal 97,
% 1197-1210, with some enhancements from klioner (2003),
% astronomical journal 125, 1580-1597.

%  tjd    = tt julian date for place (in)

%  object = character string identifying object of interest (in)

%           for solar system
%           body,             specify the name using all upper-
%                             case letters ('sun', 'moon',
%                             'jupiter', etc.),
%                             - or -
%                             specify the body id number
%                             in a 4-character string of the
%                             form '=nnn', where nnn is the
%                             body id number
%           for star,         provide a blank string, the word
%                             'star', or any string beginning
%                             with '*'

%  locatn = integer code specifying location of observer (in)

%           set locatn=0 for observer at geocenter
%           set locatn=1 for observer on surface of earth
%           set locatn=2 for observer on near-earth spacecraft

%  icoord = integer code specifying coordinate system of output
%           position (in)

%           set icoord=0 for gcrs (or 'local gcrs')
%           set icoord=1 for true equator and equinox of date
%           set icoord=2 for true equator and cio of date
%           set icoord=3 for astrometric coordinates, i.e.,
%                        without light deflection or aberration

%  star   = array of catalog data for star (in)

%           (not used if solar system body requested)

%           star(1) = icrs right ascension in hours
%           star(2) = icrs declination in degrees
%           star(3) = icrs proper motion in ra in
%                     milliarcseconds/year
%           star(4) = icrs proper motion in dec in
%                     milliarcseconds/year
%           star(5) = parallax in milliarcseconds
%           star(6) = radial velocity in kilometers/second
%           further star array elements are not used here
%           but are reserved for future use

%  observ = array of data specifying location of observer (in)

%           (not used if locatn=0)

%           for locatn = 1,

%           observ(1) = geodetic longitude (wgs-84) of observer
%                       (east +) in degrees
%           observ(2) = geodetic latitude (wgs-84) of observer
%                       (north +) in degrees
%           observ(3) = height of observer above ellipsoid
%                       in meters
%           observ(4) = value of delta-t in seconds
%                       (delta-t=tt-ut1)
%           observ(5) = (not used, reserved for future use)
%           observ(6) = (not used, reserved for future use)

%           for locatn = 2,

%           observ(1) = geocentric x in kilometers
%           observ(2) = geocentric y in kilometers
%           observ(3) = geocentric z in kilometers
%           observ(4) = geocentric x-dot in kilometers/second
%           observ(5) = geocentric y-dot in kilometers/second
%           observ(6) = geocentric z-dot in kilometers/second
%           with respect to true equator and equinox of date

%  skypos = array of output data specifying object's place
%           on the sky at time tjd, with respect to the
%           specified output coordinate system (out)

%           skypos(1) = x, dimensionless      unit vector
%           skypos(2) = y, dimensionless      toward object
%           skypos(3) = z, dimensionless
%           skypos(4) = apparent, topocentric, or astrometric
%                       right ascension in hours
%           skypos(5) = apparent, topocentric, or astrometric
%                       declination in degrees
%           skypos(6) = true (geometric, euclidian) distance
%                       to solar system body in au at time tjd,
%                       or 0.0d0 for star
%           skypos(7) = radial velocity in kilometers/second

%           further skypos array elements are not used here
%           but are reserved for future use

% note 1: values of locatn and icoord for various standard kinds of place:

%   locatn = 0 and icoord = 1 apparent place
%   locatn = 1 and icoord = 1 topocentric place
%   locatn = 0 and icoord = 0 virtual place
%   locatn = 1 and icoord = 0 local place
%   locatn = 0 and icoord = 3 astrometric place
%   locatn = 1 and icoord = 3 topocentric astrometric place

% note 2: arrays star and skypos may be expanded in the future, and
% this can be allowed for in the calling code by dimensioning
% these arrays with 20 and 10 elements, respectively, even though
% elements beyond star(6) and skypos(7) are not now referred to in
% this subroutine.

% note 3: if locatn = 1 and observ(4) = 0.0d0, the value of delta-t will
% be obtained from getdt, which provides the last value of delta-t
% defined by the user via call to setdt.

% note 4: skypos(7), the radial velocity, is the predicted
% radial velocity measure (z) times the speed of light, an
% inherently spectroscopic measure. for a star, it
% includes all effects, such as gravitational red shift,
% contained in the catalog barycentric radial velocity measure,
% which is assumed given in star(6). for a solar system
% body, it applies to a fictitious emitter at the center of the
% observed object, assumed massless (no gravitational red shift),
% and does not in general apply to reflected light.

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

% t0 = tdb julian date of epoch j2000.0 (tt)

t0 = 2451545.0d0;

% get body ids

iearth = idss_novas ('earth');

isun = idss_novas ('sun');

% get c, the speed of light in au/day

c = astcon ('c(au/day)', 1.0d0);

% check on earth as an observed object

if (strcmp(object, '=3') == 1 && locatn ~= 2)
    
    fprintf ('\n place: will not process earth as observed object except when locatn = 2');
    
    return
    
end

% compute tdbjd, the tdb julian date corresponding to ttjd

ttjd = tjd;

tdbjd = tjd;

[x, secdif] = novas_times (tdbjd);

tdbjd = ttjd + secdif / 86400.0d0;

% get position and velocity of the earth wrt barycenter of solar system, in icrs

% fprintf('\n\nearth state vector wrt barycenter\n');

[peb, veb, ierr] = solsys (tdbjd, iearth, 0);

if (ierr ~= 0)
    
    fprintf ('\nplace: cannot obtain coordinates of earth at jd %16.8f', tjd);
    
    return
    
end

% get position and velocity of the sun wrt barycenter of solar system, in icrs

% fprintf('\n\nsun state vector wrt barycenter\n');

[psb, vsb, ierr] = solsys (tdbjd, isun, 0);

if (ierr ~= 0)
    
    fprintf ('\nplace: cannot obtain coordinates of sun at jd %16.8f', tjd);
    
    return
    
end

% get position and velocity of observer

if (locatn == 1 || locatn == 2)
    
    % for topocentric place, get geocentric position and velocity
    % vectors of observer
    
    % fprintf('\n\nobserver state vector wrt geocenter\n');
    
    [pog, vog] = geopos (ttjd, locatn, observ);
    
    loc = 1;
    
else
    
    % for geocentric place, there is nothing to do
    
    for j = 1:3
        
        pog(j) = 0.0d0;
        
        vog(j) = 0.0d0;
        
    end
    
    loc = 0;
    
end

% compute position and velocity of observer wrt barycenter of
% solar system (galilean transformation from gcrs to bcrs)

for j = 1:3
    
    pob(j) = peb(j) + pog(j);
    
    vob(j) = veb(j) + vog(j);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find geometric position of observed object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (strcmp(object, 'star') == 1 || strcmp(object, ' ') == 1 || strcmp(object(1:1), '*'))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % observed object is star
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    idbody = -9999;
    
    % get position of star updated for its space motion
    
    [pos1, vel1] = vectrs (star(1), star(2), star(3), star(4), star(5), star(6));
    
    dt = dlight (pos1, pob);
    
    pos2 = propmo (t0, pos1, vel1, tdbjd + dt);
    
    % get position of star wrt observer (corrected for parallax)
    
    [pos3, tlight] = geocen (pos2, pob);
    
    dis = 0.0d0;
    
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % observed object is solar system body
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get id number of body
    
    if (strcmp(object(1:1), '=') == 1)
        
        idbody = object(2:length(object));
        
    else
        
        idbody = idss_novas (object);
        
        if (idbody == -9999)
            
            fprintf ('\nplace: cannot obtain coordinates of object at jd %16.8f', tjd);
            
            return
            
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get position of body wrt barycenter of solar system
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % fprintf('\n\nobject state vector wrt barycenter\n');
    
    [pos1, vel1, ierr] = solsys (tdbjd, str2num(idbody), 0);
    
    if (ierr ~= 0)
        
        fprintf ('\nplace: cannot obtain coordinates of object at jd %16.8f', tjd);
        
        return
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get position of body wrt observer, and true (euclidian) distance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % fprintf('\n\nobject state vector wrt observer\n');
    
    [pos2, tlight] = geocen (pos1, pob);
    
    dis = tlight * c;
    
    % get position of body wrt observer, antedated for light-time
    
    % fprintf('\n\nobject state vector wrt observer - littim\n');
    
    [pos3, tlight] = littim (tdbjd, idbody, pob, 0.0);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply gravitational deflection of light and aberration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (icoord == 3)
    
    % these calculations are skipped for astrometric place
    
    for j = 1:3
        
        pos5(j) = pos3(j);
        
    end
    
else
    
    % variable loc determines whether earth deflection is included
    
    if (loc == 1)
        
        [x, frlimb] = limang (pos3, pog);
        
        if (frlimb < 0.8d0)
            
            loc = 0;
            
        end
        
    end
    
    % compute gravitational deflection and aberration
    
    pos4 = grvdef (tdbjd, loc, pos3, pob);
    
    pos5 = aberat (pos4, vob, tlight);
    
    % position vector is now in gcrs
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transform, if necessary, to output coordinate system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (icoord == 1)
    
    % transform to equator and equinox of date
    
    pos6 = frame (pos5, 1);
    
    pos7 = preces (t0, pos6, tdbjd);
    
    pos8 = nutate (tdbjd, pos7);
    
elseif (icoord == 2)
    
    % transform to equator and cio of date
    
    % obtain the basis vectors, in the gcrs, of the celestial
    % intermediate system
    
    kcio = cioloc (tdbjd, rcio);
    
    [px, py, pz] = ciobas (tdbjd, rcio, kcio);
    
    % transform position vector to celestial intermediate system
    
    pos8(1) = px(1) * pos5(1) + px(2) * pos5(2) + px(3) * pos5(3);
    
    pos8(2) = py(1) * pos5(1) + py(2) * pos5(2) + py(3) * pos5(3);
    
    pos8(3) = pz(1) * pos5(1) + pz(2) * pos5(2) + pz(3) * pos5(3);
    
else
    
    % no transformation -- keep coordinates in gcrs (or icrs for astrometric place)
    
    for j = 1:3
        
        pos8(j) = pos5(j);
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up star data, if applicable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (idbody == -9999)
    
    rvs(1) = star(1);
    
    rvs(2) = star(2);
    
    rvs(3) = star(6);
    
    if (star(5) <= 0.0d0)
        
        vel1(1) = 0.0d0;
        
        vel1(2) = 0.0d0;
        
        vel1(3) = 0.0d0;
        
    end
    
else
    
    rvs(1) = 0.0d0;
    
    rvs(2) = 0.0d0;
    
    rvs(3) = 0.0d0;
    
end

% compute distances: observer-geocenter, observer-sun, object-sun

rvd(1) = sqrt((pob(1) - peb(1))^2 + (pob(2) - peb(2))^2 + (pob(3) - peb(3))^2);

rvd(2) = sqrt((pob(1) - psb(1))^2 + (pob(2) - psb(2))^2 + (pob(3) - psb(3))^2);

rvd(3) = sqrt((pos1(1) - psb(1))^2 + (pos1(2) - psb(2))^2 + (pos1(3) - psb(3))^2);

rv = radvl (pos3, vel1, vob, rvs, rvd);

% finish up

[ra, dec] = angles (pos8);

x = sqrt (pos8(1)^2 + pos8(2)^2 + pos8(3)^2);

for j = 1:3
    
    skypos(j) = pos8(j) / x;
    
end

skypos(4) = ra;

skypos(5) = dec;

skypos(6) = dis;

skypos(7) = rv;


