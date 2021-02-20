function [racio, k] = cioloc (tjd)

% this function returns the location of the celestial
% intermediate origin (cio) for a given julian date, as a
% right ascension with respect to either the gcrs (geocentric icrs)
% origin or the true equinox of date. the cio is always located on
% the true equator (= intermediate equator) of date.

%      tjd    = tdb julian date (in)

%      racio  = right ascension of the cio, in hours (out)

%      k      = reference system in which right ascension is
%               given (out)

%               k = 1 means gcrs

%               k = 2 means true equator and equinox of date

% note: if an external file of cio right ascensions is available,
% it will be used and k will be set to 1. otherwise an internal
% computation will be used and k will be set to 2.

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

% default to internal calculation

k = 2;

% get equation of the origins

eqor = eqxra (tjd, 1);

racio = -eqor;

end
