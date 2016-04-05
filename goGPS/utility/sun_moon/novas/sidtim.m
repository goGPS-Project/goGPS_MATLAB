function gst = sidtim (tjdh, tjdl, k)

% this function computes the greenwich sidereal time
% (either mean or apparent) at julian date tjdh + tjdl.

%      tjdh   = ut1 julian date, high-order part (in)
%      tjdl   = ut1 julian date, low-order part (in)
%               the julian date may be split at any point, but
%               for highest precision, set tjdh to be the integral
%               part of the julian date, and set tjdl to be the
%               fractional part
%      k      = time selection code (in)
%               set k=0 for greenwich mean sidereal time
%               set k=1 for greenwich apparent sidereal time
%      gst    = greenwich (mean or apparent) sidereal time
%               in hours (out)

% note:  see also subroutine setdt to set the value of delta-t
% (delta-t = tt - ut1) to be used here.

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

degcon = 180.0d0 / pi;

% t0 = tdb julian date of epoch j2000.0 (tt)

t0 = 2451545.0d0;

unitx(1) = 1.0d0;
unitx(2) = 0.0d0;
unitx(3) = 0.0d0;

deltat = getdt;

% time argument for precession and nutation components of sidereal
% time is tdb

utjd = tjdh + tjdl;

ttjd = utjd + deltat;

tdbjd = ttjd;

[a, secdif] = novas_times (tdbjd);

tdbjd = ttjd + secdif / 86400.0d0;

% get method/accuracy mode

mode = getmod;

if (mode == 2)
    
    % equinox-based mode
    % see usno circular 179, section 2.6.2
    
    % get -1 times the mean or true right ascension of the cio
    
    rcio = eqxra (tdbjd, k);
    
    % get earth rotation angle
    
    theta = erot (tjdh, tjdl);
    
    % combine to obtain sidereal time
    
    gst = mod (theta / 15.0d0 - rcio, 24.0d0);
    
    if (gst < 0.0d0)
        
        gst = gst + 24.0d0;
        
    end
    
else
    
    % cio-based mode
    % see usno circular 179, section 6.5.4
    
    % get earth rotation angle
    
    theta = erot (tjdh, tjdl);
    
    % obtain the basis vectors, in the gcrs, of the celestial
    % intermediate system
    
    [rcio, kcio] = cioloc (tdbjd);
    
    if (rcio == 99.0d0)
        
        pause
        
    end
    
    [x, y, z] = ciobas (tdbjd, rcio, kcio);
    
    % compute the direction of the true equinox in the gcrs
    
    w1 = nutate (-tdbjd, unitx);
    
    w2 = preces (tdbjd, w1, t0);
    
    eq = frame (w2, -1);
    
    % compute the hour angle of the equinox wrt the tio meridian
    % (near greenwich, but passes through the cip and tio)
    
    haeq = theta - atan2 (eq(1) * y(1) + eq(2) * y(2) + eq(3) * y(3), ...
        eq(1) * x(1) + eq(2) * x(2) + eq(3) * x(3)) * degcon;
    
    % for mean sidereal time, obtain the equation of the equinoxes
    % and subtract it
    
    if (k == 0)
        
        [a, a, ee, a, a] = etilt (tdbjd);
        
        haeq = haeq - ee / 240.0d0;
        
    end
    
    haeq = mod (haeq, 360.0d0) / 15.0d0;
    
    if (haeq < 0.0d0)
        
        haeq = haeq + 24.0d0;
        
    end
    
    gst = haeq;
    
end

