function gst = sidtim1 (tjdh, tjdl, k)

% this function computes the greenwich sidereal time
% (either mean or apparent) at julian date tjdh + tjdl

% uses low-precision nutation algorithm

% input

%  tjdh = ut1 julian date, high-order part

%  tjdl = ut1 julian date, low-order part

%         the julian date may be split at any point, but
%         for highest precision, set tjdh to be the integral
%         part of the julian date, and set tjdl to be the
%         fractional part

%  k = time selection code

%      set k=0 for greenwich mean sidereal time
%      set k=1 for greenwich apparent sidereal time

% output

%  gst = greenwich (mean or apparent) sidereal time in hours

% ported from NOVAS 3.0

%%%%%%%%%%%%%%%%%%%%%%%

% note: value of delta-t is passed via global

global deltat

% time argument for precession and nutation components
% of sidereal time is tdb

utjd = tjdh + tjdl;

ttjd = utjd + deltat;

tdbjd = ttjd;

[ttjd, secdif] = times_novas (tdbjd);

tdbjd = ttjd + secdif / 86400.0d0;

% ------------------
% equinox-based mode
% ------------------

% see usno circular 179, section 2.6.2

% get -1 times the mean or true right ascension of the cio

rcio = eqxra (tdbjd, k);

% get earth rotation angle

theta = erot (tjdh, tjdl);

% combine to obtain sidereal time

gst = mod(theta / 15.0d0 - rcio, 24.0d0);

if (gst < 0.0d0)
    gst = gst + 24.0d0;
end


