function vec2 = tercel (tjdh, tjdl, xp, yp, vec1)

% this function rotates a vector from the terrestrial to the
% celestial system. specifically, it transforms a vector in the
% itrs (a rotating earth-fixed system) to the gcrs (a local
% space-fixed system) by applying rotations for polar motion,
% earth rotation, nutation, precession, and the dynamical-to-gcrs
% frame tie.

%  tjdh   = ut1 julian date, high-order part (in)

%  tjdl   = ut1 julian date, low-order part (in)
%           the julian date may be split at any point, but
%           for highest precision, set tjdh to be the integral
%           part of the julian date, and set tjdl to be the
%           fractional part

%  xp     = conventionally-defined x coordinate of celestial
%           intermediate pole with respect to itrs pole, in arcseconds (in)

%  yp     = conventionally-defined y coordinate of celestial
%           intermediate pole with respect to itrs pole, in arcseconds (in)

%  vec1   = position vector, geocentric equatorial rectangular
%           coordinates, referred to itrs axes (terrestrial system) (in)

%  vec2   = position vector, geocentric equatorial rectangular
%           coordinates, referred to gcrs axes (celestial system) (out)

% note 1:  set xp = yp = 0.0d0 to eliminate polar motion rotation.

% note 2:  see also function setdt to set the value of delta-t
%          (delta-t = tt - ut1) to be used here.

% note 3:  both tjdh and tjdl should be non-negative for normal use
%          (tjdl=0.0d0 is ok). a negative value of tjdh is used to invoke a
%          special option where the output vector is produced with respect
%          to the equator and equinox of date, and the date for which the
%          transformation applies is taken from tjdl only. this option
%          works only in 'equinox' mode.

% note 4: input parameters xp, yp were xpole, ypole in novas f3.0.
%         the names were changed for consistancy throughout novas and with
%         iers conventions.

% t0 = tdb julian date of epoch j2000.0 (tt)

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

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

tdbjd = ttjd;

[ttjd, secdif] = novas_times (tdbjd);

tdbjd = ttjd + secdif / 86400.0d0;

% get method/accuracy mode

mode = getmod;

if (mode >= 2)
    
    %%%%%%%%%%%%%%%%
    % 'equinox' mode
    %%%%%%%%%%%%%%%%
    
    % apply polar motion
    
    if (xp == 0.0d0 && yp == 0.0d0)
        
        v1(1) = vec1(1);
        
        v1(2) = vec1(2);
        
        v1(3) = vec1(3);
        
    else
        
        v1 = wobble (tdbjd, xp, yp, vec1);
        
    end
    
    % apply earth rotation
    
    gast = sidtim (utjdh, utjdl, 1);
    
    v2 = spin (-gast * 15.0d0, v1);
    
    if (tjdh <= 0.0d0)
        
        % special option skips remaining transformations
        
        vec2(1) = v2(1);
        
        vec2(2) = v2(2);
        
        vec2(3) = v2(3);
        
    else
        
        % apply nutation and precession
        
        v3 = nutate (-tdbjd, v2);
        
        v4 = preces (tdbjd, v3, t0);
        
        % apply frame-tie matrix
        
        vec2 = frame (v4, -1);
        
    end
    
else
    
    %%%%%%%%%%%%%%%%%%%%%%
    % 'cio-tio-theta' mode
    %%%%%%%%%%%%%%%%%%%%%%
    
    % see g. kaplan (2003), 'another look at non-rotating origins',
    % proceedings of iau xxv joint discussion 16 (preprint), eq. (3) and (4).
    
    % apply polar motion, transforming the vector to the terrestrial
    % intermediate system
    
    if (xp == 0.0d0 && yp == 0.0d0)
        
        v1(1) = vec1(1);
        
        v1(2) = vec1(2);
        
        v1(3) = vec1(3);
        
    else
        
        v1 = wobble (tdbjd, xp, yp, vec1);
        
    end
    
    % obtain the basis vectors, in the gcrs, of the celestial
    % intermediate system
    
    [rcio, kcio] = cioloc (tdbjd);
    
    [x, y, z] = ciobas (tdbjd, rcio, kcio);
    
    % compute and apply the earth rotation angle theta, transforming
    % the vector to the celestial intermediate system
    
    theta = erot (utjdh, utjdl);
    
    v2 = spin (-theta, v1);
    
    % transform the vector from the celestial intermediate system
    % to the gcrs
    
    vec2(1) = x(1) * v2(1) + y(1) * v2(2) + z(1) * v2(3);
    
    vec2(2) = x(2) * v2(1) + y(2) * v2(2) + z(2) * v2(3);
    
    vec2(3) = x(3) * v2(1) + y(3) * v2(2) + z(3) * v2(3);
    
end
