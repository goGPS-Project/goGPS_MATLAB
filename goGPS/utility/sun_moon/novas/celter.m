function vec2 = celter (tjdh, tjdl, xp, yp, vec1)

% this function rotates a vector from the celestial to the
% terrestrial system.  specifically, it transforms a vector in the
% gcrs (a local space-fixed system) to the itrs (a rotating
% earth-fixed system) by applying rotations for the gcrs-to-
% dynamical frame tie, precession, nutation, earth rotation,
% and polar motion.

%      tjdh   = ut1 julian date, high-order part (in)

%      tjdl   = ut1 julian date, low-order part (in)
%               the julian date may be split at any point, but
%               for highest precision, set tjdh to be the integral
%               part of the julian date, and set tjdl to be the
%               fractional part

%      xp     = conventionally-defined x coordinate of celestial
%               intermediate pole with respect to itrs pole,
%               in arcseconds (in)

%      yp     = conventionally-defined y coordinate of celestial
%               intermediate pole with respect to itrs pole,
%               in arcseconds (in)

%      vec1   = position vector, geocentric equatorial rectangular
%               coordinates, referred to gcrs axes (celestial
%               system) (in)

%      vec2   = position vector, geocentric equatorial rectangular
%               coordinates, referred to itrs axes (terrestrial
%               system) (out)

% note 1:  set xp = yp = 0.0d0 to eliminate polar motion rotation.

% note 2:  see also subroutine setdt to set the value of delta-t
%          (delta-t = tt - ut1) to be used here.

% note 3:  both tjdh and tjdl should be non-negative for normal use
%          (tjdl=0.0d0 is ok). a negative value of tjdh is used to invoke a
%          special option where the input vector is assumed to be with
%          respect to the equator and equinox of date, and the date for which
%          the transformation applies is taken from tjdl only. this option
%          works only in 'equinox' mode.

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

% t0 = tdb julian date of epoch j2000.0 (tt)

t0 = 2451545.0d0;

deltat = getdt;

if (tjdh >= 0.0d0)
    
    utjdh = tjdh;
    
    utjdl = tjdl;
    
else
    
    utjdh = tjdl;
    
    utjdl = 0.0d0;
    
end

utjd = utjdh + utjdl;

% time argument for precession and nutation is tdb

ttjd = utjd + deltat;

tdbjd = utc2tdb(ttjd);

% tdbjd = ttjd;
% 
% [t, secdif] = novas_times (tdbjd);
% 
% tdbjd = ttjd + secdif / 86400.0d0;

% get method/accuracy mode

mode = getmod;

if (mode >= 2)
    
    % 'equinox' mode
    
    % special option skips initial transformations
    
    if (tjdh < 0.0d0)
        
        v3(1) = vec1(1);
        
        v3(2) = vec1(2);
        
        v3(3) = vec1(3);
        
    else
        
        % apply frame-tie matrix
        
        v1 = frame (vec1, 1);
        
        % apply precession and nutation
        
        v2 = preces (t0, v1, tdbjd);
        
        v3 = nutate (tdbjd, v2);
        
    end
    
    % apply earth rotation
    
    gast = sidtim (utjdh, utjdl, 1);
    
    v4 = spin (gast * 15.0d0, v3);
    
    % apply polar motion
    
    if (xp == 0.0d0 && yp == 0.0d0)
        
        vec2(1) = v4(1);
        
        vec2(2) = v4(2);
        
        vec2(3) = v4(3);
        
    else
        
        vec2 = wobble (-tdbjd, xp, yp, v4);
        
    end
    
else
    
    % 'cio-tio-theta' mode
    % see g. kaplan (2003), 'another look at non-rotating origins',
    % proceedings of iau xxv joint discussion 16 (preprint),
    % eq. (3) and (4).
    
    % obtain the basis vectors, in the gcrs, of the celestial
    % intermediate system
    
    [rcio, kcio] = cioloc (tdbjd);
    
    [x, y, z] = ciobas (tdbjd, rcio, kcio);
    
    % transform the vector from the gcrs to the
    % celestial intermediate system
    
    v1(1) = x(1) * vec1(1) + x(2) * vec1(2) + x(3) * vec1(3);
    
    v1(2) = y(1) * vec1(1) + y(2) * vec1(2) + y(3) * vec1(3);
    
    v1(3) = z(1) * vec1(1) + z(2) * vec1(2) + z(3) * vec1(3);
    
    % compute and apply the earth rotation angle theta, transforming
    % the vector to the terrestrial intermediate system
    
    theta = erot (utjdh, utjdl);
    
    v2 = spin (theta, v1);
    
    % apply polar motion, transforming the vector to the itrs
    
    if (xp == 0.0d0 && yp == 0.0d0)
        
        vec2(1) = v2(1);
        
        vec2(2) = v2(2);
        
        vec2(3) = v2(3);
        
    else
        
        vec2 = wobble (-tdbjd, xp, yp, v2);
        
    end
    
end

