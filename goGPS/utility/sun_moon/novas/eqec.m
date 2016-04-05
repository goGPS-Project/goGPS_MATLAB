function pos2 = eqec (tjd, icoord, pos1)

% this function converts an equatorial position vector to
% an ecliptic position vector.

% input

%  tjd    = tt julian date of equator, equinox, and ecliptic
%           used for coordinates

%  icoord = coordinate system selection

%           set icoord = 0 for mean equator and equinox of date
%           set icoord = 1 for true equator and equinox of date

%           (ecliptic is always the mean plane)

%  pos1   = position vector, referred to specified
%           equator and equinox of date

% output

%  pos2 = position vector, referred to specified
%         ecliptic and equinox of date

% note: to convert icrs vectors to ecliptic vectors (mean ecliptic
% and equinox of j2000.0 only), set tjd = 0.d0 and icoord = 0.
% except for the input to this case, all vectors are assumed to
% be with respect to a dynamical system.

% important: initialize tlast and ob2000 in main script

% ported from NOVAS 3.0

%%%%%%%%%%%%%%%%%%%%%%%

global tlast ob2000

radcon = pi / 180.0d0;

% t0 = tdb julian date of epoch j2000.0 (tt)

t0 = 2451545.0d0;

% t1 is the tdb julian date

[t, secdif] = times_novas (tjd);

t1 = tjd + secdif / 86400.0d0;

if (tjd == 0.0d0)

    % case where input vector is in icrs system

    pos2 = frame (pos1, 1);

    % get mean obliquity at j2000.0 if necessary

    if (ob2000 == 0.0d0)

        [ob2000, x, x, x, x] = etilt (t0);

    end

    obl = ob2000 * radcon;

else

    % case where input vector is in equator of date system

    pos0(1) = pos1(1);

    pos0(2) = pos1(2);

    pos0(3) = pos1(3);

    % get mean and true obliquity

    if (abs(tjd - tlast) > 1.0d-8)

        [oblm, oblt, x, x, x] = etilt (t1);

        tlast = tjd;
    end

    % select mean or true obliquity

    obl = oblm * radcon;

    if (icoord == 1)

        obl = oblt * radcon;

    end

end

% rotate equatorial position vector to ecliptic system

pos2(1) =  pos0(1);

pos2(2) =  pos0(2) * cos(obl) + pos0(3) * sin(obl);

pos2(3) = -pos0(2) * sin(obl) + pos0(3) * cos(obl);


