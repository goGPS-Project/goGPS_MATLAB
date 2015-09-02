function [XS_pc, Az, El, D, time] = satellite_phase_center(XS_cm)

global iephem km ephname inutate psicor epscor ob2000

% compute the Sun position
sun_id = 11; earth_id = 3;
readleap; iephem = 1; ephname = 'de421.bin'; km = 1; inutate = 1; ob2000 = 0.0d0;
setmod(2);
setdt(3020092e-7);

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
for day = 1 : 0.001 : 2
    jdutc = julian(month, day, year);
    jdtdb = utc2tdb(jdutc);
    rrd = jplephem(jdtdb, sun_id, earth_id);
    
    %ICRS coordinates
    sun_ECI = rrd(1:3);
    tmatrix = j2000_icrs(1);
    sun_ECI = tmatrix*sun_ECI;
    
    %precise celestial pole (disabled)
    [psicor, epscor] = celpol(jdtdb, 1, 0.0d0, 0.0d0);
    
    %ICRS coordinates to ITRS coordinates
    xp = 171209e-6; yp = 414328e-6;
    deltat = getdt;
    jdut1 = jdutc - deltat;
    tjdh = floor(jdut1); tjdl = jdut1 - tjdh;
    sun_ECEF(:,i) = celter(tjdh, tjdl, xp, yp, sun_ECI);
    
    %topocentric coordinates
    [Az(i), El(i), D(i)] = topocent([4.398834576291538e+06; 7.006425011286493e+05; 4.549777495792151e+06], sun_ECEF(:,i)'*1e3);
    
    [gps_week, gps_sow] = jd2gps(jdutc);
    gps_date = gps2date(gps_week, gps_sow);
    time(i) = datenum(gps_date);
    
    i = i + 1;
end

XS_pc = sun_ECEF;
