function deltat = getdt

% this entry returns the current value of delta-t
% (delta-t = tt - ut1), previously set by the user.  the value
% returned is to be used in the calculation of sidereal time and
% the terrestrial-to-celestial transformation.  it allows these
% calculations to be performed, correctly, using ut1 as the time
% argument for the earth rotation angle and tdb as the time argument
% for the precession and nutation components.

%  deltat = value of delta-t (tt-ut1) in days (out)

global dt

deltat = dt;
