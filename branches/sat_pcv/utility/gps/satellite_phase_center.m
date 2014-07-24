function [XS_pc, Az, El, D, time] = satellite_phase_center(XS_cm)

global iephem km ephname inutate psicor epscor ob2000

% compute the Sun position (ICRS)
sun_id = 11; earth_id = 3;
readleap; iephem = 1; ephname = 'de421.bin'; km = 1; inutate = 1; ob2000 = 0.0d0;

%if the binary JPL ephemeris file is not available, generate it
if (~exist('./utility/sun_moon/jpl_ephem/de421.bin','file'))
    fprintf('Warning: file "de421.bin" not found in directory ./utility/sun_moon/jpl_ephem/ ... generating a new "de421.bin" file\n')
    fprintf('         (this procedure should be done only once on each machine):\n')
    fprintf('-------------------------------------------------------------------\n\n')
    asc2eph(421, {'ascp1900.421', 'ascp2050.421'}, './utility/sun_moon/jpl_ephem/de421.bin');
    fprintf('-------------------------------------------------------------------\n\n')
end

year = 2014;
month = 7;

i = 1;
for day = 1 : 0.01 : 2
    jdutc = julian(month, day, year);
    jdtdb = utc2tdb(jdutc);
    rrd = jplephem(jdtdb, sun_id, earth_id);
    
    %ICRS coordinates (ecliptic)
    sun_ECEF = rrd(1:3);
    
    %ICRS ecliptic to equatorial coordinates
    sun_ECEF = eceq(0.0d0, 0, sun_ECEF);
    
    %precise celestial pole (now disabled)
    [psicor, epscor] = celpol(jdtdb, 1, 0.0d0, 0.0d0);
    
    %ICRS equatorial coordinates to ITRS equatorial coordinates
    xp = 0.0d0; yp = 0.0d0;
    tjdh = floor(jdutc); tjdl = jdutc - tjdh;
    setdt(64);
    TRF(:,i) = celter(tjdh, tjdl, xp, yp, sun_ECEF);
    
    %topocentric coordinates
    [Az(i), El(i), D(i)] = topocent([4.398972273573140e+06; 7.006644334305141e+05; 4.549920877913681e+06], TRF(:,i)');
    
    [gps_week, gps_sow] = jd2gps(jdutc);
    gps_date = gps2date(gps_week, gps_sow);
    time(i) = datenum(gps_date);
    
    i = i + 1;
end

XS_pc = TRF;
