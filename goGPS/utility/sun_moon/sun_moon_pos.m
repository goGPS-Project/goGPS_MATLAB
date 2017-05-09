function [sun_ECEF, moon_ECEF] = sun_moon_pos(epoch)

global iephem km ephname inutate psicor epscor ob2000

year = epoch(:,1);
month = epoch(:,2);
day = epoch(:,3) + epoch(:,4)/24 + epoch(:,5)/1440 + epoch(:,6)/86400;

sun_id = 11; moon_id = 10; earth_id = 3;

readleap; iephem = 1; ephname = 'de436.bin'; km = 1; inutate = 1; ob2000 = 0.0d0;

tmatrix = j2000_icrs(1);

setmod(2);
% setdt(3020092e-7);
setdt(5.877122033683494);
xp = 171209e-6; yp = 414328e-6;

gs = Go_State.getInstance();
go_dir = gs.getLocalStorageDir();

%if the binary JPL ephemeris file is not available, generate it
if (exist(fullfile(go_dir, 'de436.bin'),'file') ~= 2)
    fprintf('Warning: file "de436.bin" not found in at %s\n         ... generating a new "de436.bin" file\n',fullfile(go_dir, 'de436.bin'));
    fprintf('         (this procedure may take a while, but it will be done only once on each installation):\n')
    fprintf('-------------------------------------------------------------------\n\n')
    asc2eph(436, {'ascp01950.436', 'ascp02050.436'}, fullfile(go_dir, 'de436.bin'));
    fprintf('-------------------------------------------------------------------\n\n')
end

sun_ECEF = zeros(3,length(year));
moon_ECEF = zeros(3,length(year));

for e = 1 : length(year)

    %UTC to TDB
    jdutc = julian(month(e,1), day(e,1), year(e,1));
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
    sun_ECEF(:,e) = celter(tjdh, tjdl, xp, yp, sun_ECI);

    if (nargout > 1)
        %compute the Moon position (ICRS coordinates)
        rrd = jplephem(jdtdb, moon_id, earth_id);
        moon_ECI = rrd(1:3);
        moon_ECI = tmatrix*moon_ECI;

        %Moon ICRS coordinates to ITRS coordinates
        deltat = getdt;
        jdut1 = jdutc - deltat;
        tjdh = floor(jdut1); tjdl = jdut1 - tjdh;
        moon_ECEF(:,e) = celter(tjdh, tjdl, xp, yp, moon_ECI);
    end
end

sun_ECEF  = sun_ECEF*1e3;
moon_ECEF = moon_ECEF*1e3;
