function rrd = jplephem (et, ntarg, ncent)

% reads the jpl planetary ephemeris and gives
% the position and velocity of the point 'ntarg'
% with respect to point 'ncent'

% input

%   et    = julian ephemeris date at which interpolation is wanted

%   ntarg = integer number of 'target' point

%   ncent = integer number of center point

%   the numbering convention for 'ntarg' and 'ncent' is:

%        1 = mercury           8 = neptune
%        2 = venus             9 = pluto
%        3 = earth            10 = moon
%        4 = mars             11 = sun
%        5 = jupiter          12 = solar-system barycenter
%        6 = saturn           13 = earth-moon barycenter
%        7 = uranus           14 = nutations (longitude and obliq)
%                             15 = librations, if on ephemeris file

%        if nutations are wanted, set ntarg = 14.
%        for librations, set ntarg = 15. set ncent = 0.

% output

%   rrd = output 6-word array containing position and velocity
%         of point 'ntarg' relative to 'ncent'. the units are au and
%         au/day. for librations the units are radians and radians
%         per day. in the case of nutations the first four words of
%         rrd will be set to nutations and rates, having units of
%         radians and radians/day.

%         the option is available to have the units in km and km/sec.
%         for this, set km = 1 via global in the calling program.

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global jplephem_cval jplephem_ss jplephem_au jplephem_emrat jplephem_ncon jplephem_ipt jplephem_np jplephem_nv ...
    jplephem_twot jplephem_pc jplephem_vc

global jplephem_iephem jplephem_ephname jplephem_bary jplephem_pvsun jplephem_nrl jplephem_fid jplephem_lpt

rrd = zeros(6, 1);

list = zeros(12, 1);

et2 = zeros(2, 1);

% load time array

et2(1) = et;

et2(2) = 0;

% first entry?

if (jplephem_iephem == 1)

    jplephem_pvsun = zeros(6, 1);

    % read header file data

    go_dir = Core.getLocalStorageDir();
    jplephem_fid = fopen(fullfile(go_dir, jplephem_ephname), 'r');

    ttl = fread(jplephem_fid, 252);

    cnam = fread(jplephem_fid, 2400);

    jplephem_ss = fread(jplephem_fid, 3, 'double');

    jplephem_ncon = fread(jplephem_fid, 1, 'int');

    % astronomical unit

    jplephem_au = fread(jplephem_fid, 1, 'double');

    % earth-moon ratio

    jplephem_emrat = fread(jplephem_fid, 1, 'double');

    jplephem_ipt = fread(jplephem_fid, [3 12], 'int');

    numde = fread(jplephem_fid, 1, 'int');

    jplephem_lpt = fread(jplephem_fid, 3, 'int');

    % move to next record

    status = fseek(jplephem_fid, 8144, 'bof');

    % read "constant" values

    jplephem_cval = fread(jplephem_fid, 400, 'double');

    % initialization

    jplephem_nrl = 0;

    jplephem_bary = 0;

    jplephem_pc(1) = 1;
    jplephem_pc(2) = 0;

    jplephem_vc(2) = 1;

    jplephem_np = 2;
    jplephem_nv = 3;

    jplephem_twot = 0;

    jplephem_iephem = 0;
end

if (ntarg == ncent)

    return;

end

%%%%%%%%%%%%%%%%%%%%%%
% nutations requested?
%%%%%%%%%%%%%%%%%%%%%%

if (ntarg == 14)

    if (jplephem_ipt(2, 12) > 0)

        list(11) = 2;

        [pv, rrd] = state(et2, list);

        list(11) = 0;

        return;

    else

        fprintf('\n\njplephem - no nutations on this ephemeris file \n');

        return;

    end
end

%%%%%%%%%%%%%%%%%%%%%%%
% librations requested?
%%%%%%%%%%%%%%%%%%%%%%%

if (ntarg == 15)

    if (jplephem_lpt(2) > 0)

        list(12) = 2;

        [pv, rrd] = state(et2, list);

        list(12) = 0;

        for i = 1:1:6
            rrd(i) = pv(i, 11);
        end

        return

    else

        fprintf('\n\n no librations on this ephemeris file \n');

        return;

    end

end

% force barycentric output by function 'state'

bsave = jplephem_bary;

jplephem_bary = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up proper entries in 'list' array for state call
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:1:2

    k = ntarg;

    if (i == 2)
        k = ncent;
    end

    if (k <= 10)
        list(k) = 2;
    end

    if (k == 10)
        list(3) = 2;
    end

    if (k == 3)
        list(10) = 2;
    end

    if (k == 13)
        list(3) = 2;
    end

end

%%%%%%%%%%%%%%%%%%%%
% make call to state
%%%%%%%%%%%%%%%%%%%%

[pv, rrd] = state(et2, list);

if (ntarg == 11 || ncent == 11)
    for i = 1:1:6
        pv(i, 11) = jplephem_pvsun(i);
    end
end

if (ntarg == 12 || ncent == 12)
    for i = 1:1:6
        pv(i, 12) = 0;
    end
end

if (ntarg == 13 || ncent == 13)
    for i = 1:1:6
        pv(i, 13) = pv(i, 3);
    end
end

if (ntarg * ncent == 30 && ntarg + ncent == 13)
    for i = 1:1:6
        pv(i, 3) = 0;
    end
else
    if (list(3) == 2)
        for i = 1:1:6
            pv(i, 3) = pv(i, 3) - pv(i, 10) / (1 + jplephem_emrat);
        end
    end

    if (list(10) == 2)
        for i = 1:1:6
            pv(i, 10) = pv(i, 3) + pv(i, 10);
        end
    end
end

for i = 1:1:6
    rrd(i) = pv(i, ntarg) - pv(i, ncent);
end

jplephem_bary = bsave;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pv, nut] = state(et2, list)

% reads and interpolates the jpl planetary ephemeris file

% input

%   et2    2-word julian ephemeris epoch at which interpolation
%          is wanted.  any combination of et2(1) + et2(2) which falls
%          within the time span on the file is a permissible epoch.

%          a. for ease in programming, the user may put the
%             entire epoch in et2(1) and set et2(2) = 0.

%          b. for maximum interpolation accuracy, set et2(1) equal
%             to the most recent midnight at or before interpolation
%             epoch and set et2(2) equal to fractional part of a day
%             elapsed between et2(1) and epoch.

%          c. as an alternative, it may prove convenient to set
%             et2(1) = some fixed epoch, such as start of integration,
%             and et2(2) = elapsed interval between and epoch.

%   list   12-word integer array specifying what interpolation
%          is wanted for each of the bodies on the file.
%
%          list(i) = 0 => no interpolation for body i
%                  = 1 => position only
%                  = 2 => position and velocity

%          the designation of the astronomical bodies by i is:

%              i =  1 => mercury
%                =  2 => venus
%                =  3 => earth-moon barycenter
%                =  4 => mars
%                =  5 => jupiter
%                =  6 => saturn
%                =  7 => uranus
%                =  8 => neptune
%                =  9 => pluto
%                = 10 => geocentric moon
%                = 11 => nutations in longitude and obliquity
%                = 12 => lunar librations (if on file)

% output

%   pv   6 x 11 array that will contain requested interpolated
%        quantities. the body specified by list(i) will have its
%        state in the array starting at pv(1, i). (on any given
%        call, only those words in 'pv' which are affected by the
%        first 10 'list' entries (and by list(12) if librations are
%        on the file) are set. the rest of the 'pv' array
%        is untouched).  the order of components starting in
%        pv(1, i) is x, y, z, dx, dy, dz.

%        all output vectors are referenced to the earth mean
%        equator and equinox of j2000 if the de number is 200 or
%        greater; of b1950 if the de number is less than 200.

%        the moon state is always geocentric; the other nine states
%        are either heliocentric or solar-system barycentric,
%        depending on the setting of common flags (see below).

%        lunar librations, if on file, are put into pv(k, 11) if
%        list(12) is 1 or 2.
%
%   nut  4-word array that will contain nutations and rates,
%        depending on the setting of list(11). the order of
%        quantities in nut is:

%        d psi  (nutation in longitude)
%        d epsilon (nutation in obliquity)
%        d psi dot
%        d epsilon dot

% global

%   km    logical flag defining physical units of the output states
%         = 1 => kilometers and kilometers/second
%         = 0 => au and au/day

%         default value = 0 (km determines time unit
%         for nutations and librations. angle unit is always radians.)

%   jplephem_bary  logical flag defining output center.
%         only the 9 planets are affected.
%         jplephem_bary = 1 => center is solar-system barycenter
%              = 0 => center is sun
%         default value = 0

%   jplephem_pvsun 6-word array containing the jplephem_barycentric position and
%         velocity of the sun

global jplephem_ss jplephem_au jplephem_ipt jplephem_lpt

global jplephem_nrl jplephem_fid jplephem_km jplephem_bary jplephem_pvsun jplephem_coef

nut = zeros(4, 1);

pv = zeros(6, 11);

if (et2(1) == 0)
    return
end

s = et2(1) - 0.5;

tmp = splitDouble(s);

pjd(1) = tmp(1);

pjd(2) = tmp(2);

tmp = splitDouble(et2(2));

pjd(3) = tmp(1);

pjd(4) = tmp(2);

pjd(1) = pjd(1) + pjd(3) + 0.5;

pjd(2) = pjd(2) + pjd(4);

tmp = splitDouble(pjd(2));

pjd(3) = tmp(1);

pjd(4) = tmp(2);

pjd(1) = pjd(1) + pjd(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% error return for epoch out of range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (pjd(1) + pjd(4) < jplephem_ss(1) || pjd(1) + pjd(4) > jplephem_ss(2))
    fprintf('\n\n error in state - epoch out of range error in jplephem\n');
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate record number and relative time in interval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nr = fix((pjd(1) - jplephem_ss(1)) / jplephem_ss(3)) + 2;

if (pjd(1) == jplephem_ss(2))
    nr = nr - 1;
end

t(1) = ((pjd(1) - ((nr-2) * jplephem_ss(3) + jplephem_ss(1))) + pjd(4)) / jplephem_ss(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read correct record if not in core
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nr ~= jplephem_nrl)
    jplephem_nrl = nr;

    status = fseek(jplephem_fid, nr * 8144, 'bof');

    jplephem_coef = fread(jplephem_fid, 1018, 'double');
end

if (jplephem_km == 1)
    t(2) = jplephem_ss(3) * 86400;
    aufac = 1;
else
    t(2) = jplephem_ss(3);
    aufac = 1 / jplephem_au;
end

% interpolate barycentric state vector of sun

tmpv = zeros(3, 2);

ibuf = jplephem_ipt(1, 11);
ncf = jplephem_ipt(2, 11);
na = jplephem_ipt(3, 11);

tmpv = interp(ibuf, t, ncf, 3, na, 2);

k = 0;

for j = 1:1:2
    for i = 1:1:3
        k = k + 1;
        jplephem_pvsun(k) = tmpv(i, j) * aufac;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check and interpolate bodies requested
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmpv = zeros(3, 2);

tempv = zeros(6, 1);

for i = 1:1:10
    if (list(i) ~= 0)

        ibuf = jplephem_ipt(1, i);
        ncf = jplephem_ipt(2, i);
        na = jplephem_ipt(3, i);

        tmpv = interp(ibuf, t, ncf, 3, na, list(i));

        k = 0;

        for j = 1:1:2
            for m = 1:1:3
                k = k + 1;
                tempv(k) = tmpv(m, j);
            end
        end

        for j = 1:1:6
            if (i <= 9 && jplephem_bary == 0)
                pv(j, i) = tempv(j) * aufac - jplephem_pvsun(j);
            else
                pv(j, i) = tempv(j) * aufac;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do nutations if requested (and if on file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (list(11) > 0 && jplephem_ipt(2, 12) > 0)

    ibuf = jplephem_ipt(1, 12);
    ncf = jplephem_ipt(2, 12);
    na = jplephem_ipt(3, 12);

    nut = interp(ibuf, t, ncf, 2, na, list(11));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get librations if requested (and if on file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (list(12) > 0 && jplephem_lpt(2) > 0)
    tmpv = interp(jplephem_lpt(1), t, jplephem_lpt(2), 3, jplephem_lpt(3), list(12));

    for i = 1:1:3
        pv(i, 11) = tmpv(i, 1);

        pv(i + 3, 11) = tmpv(i, 2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pv = interp(ibuf, t, ncf, ncm, na, ifl)

% this function differentiates and interpolates a set
% of chebyshev coefficients to give position and velocity

% input

%   buf   1st location of array of chebyshev coefficients of position

%   t     t(1) is fractional time in interval covered by
%         coefficients at which interpolation is wanted
%         (0 <= t(1) <= 1). t(2) is length of whole
%         interval in input time units.

%   ncf   number of coefficients per component

%   ncm   number of components per set of coefficients

%   na    number of sets of coefficients in full array
%         (i.e., number of sub-intervals in full interval)

%   ifl   integer flag
%         = 1 for positions only
%         = 2 for pos and vel

% output

%   pv    interpolated quantities requested.  dimension
%         expected is pv(ncm, ifl)

global jplephem_coef jplephem_np jplephem_nv jplephem_twot jplephem_pc jplephem_vc

pv = zeros(6, 12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get correct sub-interval number for this set of
% coefficients and get normalized chebyshev time
% within that subinterval.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dna = na;

dt1 = fix(t(1));

temp = dna * t(1);

ll = fix(temp - dt1) + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tc is the normalized chebyshev time (-1 <= tc <= 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tc = 2 * (mod(temp, 1) + dt1) - 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check to see whether chebyshev time has changed,
% and compute new polynomial values if it has.
% (the element jplephem_pc(2) is the value of t1(tc) and hence
% contains the value of tc on the previous call)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (tc ~= jplephem_pc(2))
    jplephem_np = 2;
    jplephem_nv = 3;
    jplephem_pc(2) = tc;
    jplephem_twot = tc + tc;
end

% be sure that at least 'ncf' polynomials have been evaluated
% and are stored in the array 'jplephem_pc'.

if (jplephem_np < ncf)
    for i = jplephem_np + 1:1:ncf
        jplephem_pc(i) = jplephem_twot * jplephem_pc(i - 1) - jplephem_pc(i - 2);
    end

    jplephem_np = ncf;
end

bcoef = ncf * na * ncm;

cbody = zeros(bcoef, 1);

n = ibuf;

for m = 1:1:bcoef
    cbody(m) = jplephem_coef(n);
    n = n + 1;
end

cbuf = zeros(ncf, ncm, na);

n = 0;

for l = 1:1:na
    for i = 1:1:ncm
        for j = 1:1:ncf
            n = n + 1;
            cbuf(j, i, l) = cbody(n);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolate to get position for each component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:1:ncm
    pv(i, 1) = 0;

    for j = ncf:-1:1
        pv(i, 1) = pv(i, 1) + jplephem_pc(j) * cbuf(j, i, ll);
    end
end

if (ifl <= 1)
    return
end

% if velocity interpolation is wanted, be sure enough
% derivative polynomials have been generated and stored.

vfac = (dna + dna) / t(2);

jplephem_vc(3) = jplephem_twot + jplephem_twot;

if (jplephem_nv < ncf)
    for i = jplephem_nv + 1:1:ncf
        jplephem_vc(i) = jplephem_twot * jplephem_vc(i - 1) + jplephem_pc(i - 1) + jplephem_pc(i - 1) - jplephem_vc(i - 2);
    end

    jplephem_nv = ncf;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolate to get velocity for each component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:1:ncm
    pv(i, 2) = 0;

    for j = ncf:-1:2
        pv(i, 2) = pv(i, 2) + jplephem_vc(j) * cbuf(j, i, ll);
    end

    pv(i, 2) = pv(i, 2) * vfac;
end

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

function fr = splitDouble(tt)

% this function breaks a number into a integer
% and a fractional part.

% input

%   tt = input number

% output

%   fr = 2-word output array
%        fr(1) contains integer part
%        fr(2) contains fractional part

%        for negative input numbers, fr(1) contains the next
%        more negative integer; fr(2) contains a positive fraction.

fr = zeros(2, 1);

fr(1) = fix(tt);

fr(2) = tt - fr(1);

if (tt >= 0 || fr(2) == 0)
    return
end

% make adjustments for negative input number

fr(1) = fr(1) - 1;

fr(2) = fr(2) + 1;


