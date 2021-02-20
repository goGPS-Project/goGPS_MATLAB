function setdt (delt)

% this function allows for the specification of the delta-t value
% (delta-t = tt - ut1) to be used in the calculation of sidereal
% time and the terrestrial-to-celestial transformation.  it allows
% these calculations to be performed, correctly, using ut1 as the
% time argument for the earth rotation angle and tdb as the time
% argument for the precession and nutation components.  this
% function, if used, should be called before any function
% related to earth rotation (e.g., sidtim or tercel) for a given
% date.  the value of delta-t specified here will be used until
% explicitly changed.

%  delt = value of delta-t (tt-ut1) in seconds (in)

% note 1:  the computed value of sidereal time, and the equivalent
% earth orientation angles, are relatively insensitive to the value
% of delta-t: up to only ~3 microarcseconds per second of delta-t.
% therefore, for many applications, this function either need not
% be called at all, or can be called just once for a wide range of
% dates (e.g., a year).  if this call is not used, a default
% value of delta-t of 64 seconds is used, which is appropriate to
% 2000.0.

% note 2:  the input time arguments to sidtim and tercel (tjdh and
% tjdl) are expressed in ut1 regardless of whether this call is
% used.

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

global dt

dt = delt / 86400.0d0;

end


