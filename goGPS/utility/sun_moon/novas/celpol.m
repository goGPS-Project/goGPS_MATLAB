function [psicor, epscor] = celpol (tjd, itype, dpole1, dpole2)

% this function allows for the specification of celestial pole
% offsets for high-precision applications. each set of offsets is
% a correction to the modeled position of the pole for a specific
% date, derived from observations and published by the iers.
% this entry, if used, should be called before any other routines
% for a given date. values of the pole offsets specified via a call
% to this entry will be used until explicitly changed.

%      tjd    = tdb or tt julian date for pole offsets (in)

%      itype  = type of pole offset (in)

%               set itype = 1 for corrections to angular coordinates
%                             of modeled pole referred to mean
%                             ecliptic of date, that is,
%                             delta-delta-psi and delta-delta-epsilon

%               set itype = 2 for corrections to components of
%                             modeled pole unit vector with referred
%                             to gcrs axes, that is, dx and dy

%      dpole1 = value of celestial pole offset in first coordinate,
%               (delta-delta-psi or dx) in milliarcseconds (in)

%      dpole2 = value of celestial pole offset in second coordinate,
%               (delta-delta-epsilon or dy) in milliarcseconds (in)

% note 1: tjd is used only for itype = 2, to transform dx and dy to
%         the equivalent delta-delta-psi and delta-delta-epsilon values.

% note 2: for itype = 2, dx and dy are unit vector component
%         corrections, but are expressed in milliarcseconds simply by
%         multiplying by 206264806, the number of milliarcseconds in one
%         radian.

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

seccon = 180.d0 * 3600.d0 / pi;

t0 = 2451545.0d0;

if (itype == 1)
    
    psicor = dpole1 * 1.0d-3;
    
    epscor = dpole2 * 1.0d-3;
    
else
    
    dx = dpole1;
    
    dy = dpole2;
    
    t = (tjd - t0) / 36525.0d0;
    
    % compute sine of mean obliquity of date
    
    sine = sin (obliq(t) / seccon);
    
    % the following algorithm, to transform dx and dy to
    % delta-delta-psi and delta-delta-epsilon, is from g. kaplan
    % (2003), usno/aa technical note 2003-03, eqs. (7)-(9).
    
    % trivial model of pole trajectory in gcrs allows computation of dz
    
    x = (2004.19d0 * t) / seccon;
    
    dz = - (x + 0.5d0 * x^3) * dx;
    
    % form pole offset vector (observed - modeled) in gcrs
    
    dp1(1) = dx * 1.0d-3 / seccon;
    
    dp1(2) = dy * 1.0d-3 / seccon;
    
    dp1(3) = dz * 1.0d-3 / seccon;
    
    % precess pole offset vector to mean equator and equinox of date
    
    dp2 = frame (dp1, 1);
    
    dp3 = preces (t0, dp2, tjd);
    
    % compute delta-delta-psi and delta-delta-epsilon in arcseconds
    
    psicor = (dp3(1) / sine) * seccon;
    
    epscor = (dp3(2)) * seccon;
    
end

