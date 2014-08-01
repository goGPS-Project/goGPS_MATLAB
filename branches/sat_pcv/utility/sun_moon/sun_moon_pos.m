function [sun_ECEF, moon_ECEF] = sun_moon_pos(year, month, day)

global iephem km ephname inutate psicor epscor ob2000

sun_id = 11; moon_id = 10; earth_id = 3;

readleap; iephem = 1; ephname = 'de421.bin'; km = 1; inutate = 1; ob2000 = 0.0d0;

tmatrix = j2000_icrs(1);

setmod(2);
setdt(3020092e-7);
xp = 171209e-6; yp = 414328e-6;

%if the binary JPL ephemeris file is not available, generate it
if (~exist('./utility/sun_moon/jpl_ephem/de421.bin','file'))
    fprintf('Warning: file "de421.bin" not found in directory ./utility/sun_moon/jpl_ephem/ ... generating a new "de421.bin" file\n')
    fprintf('         (this procedure should be done only once on each machine):\n')
    fprintf('-------------------------------------------------------------------\n\n')
    asc2eph(421, {'ascp1900.421', 'ascp2050.421'}, './utility/sun_moon/jpl_ephem/de421.bin');
    fprintf('-------------------------------------------------------------------\n\n')
end

jdutc = julian(month, day, year);
jdtdb = utc2tdb(jdutc);

%precise celestial pole (disabled)
[psicor, epscor] = celpol(jdtdb, 1, 0.0d0, 0.0d0);

%compute the Sun position (ICRS coordinates)
rrd = jplephem(jdtdb, sun_id, earth_id);
sun_ECI = rrd(1:3);
sun_ECI = tmatrix*sun_ECI;

%Sun ICRS coordinates to ITRS coordinates
deltat = getdt;
jdut1 = jdutc - deltat;
tjdh = floor(jdut1); tjdl = jdut1 - tjdh;
sun_ECEF = celter(tjdh, tjdl, xp, yp, sun_ECI);

if (nargout > 1)
    %compute the Moon position (ICRS coordinates)
    rrd = jplephem(jdtdb, moon_id, earth_id);
    moon_ECI = rrd(1:3);
    moon_ECI = tmatrix*moon_ECI;
    
    %Moon ICRS coordinates to ITRS coordinates
    deltat = getdt;
    jdut1 = jdutc - deltat;
    tjdh = floor(jdut1); tjdl = jdut1 - tjdh;
    moon_ECEF = celter(tjdh, tjdl, xp, yp, moon_ECI);
end

