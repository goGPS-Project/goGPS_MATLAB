%   CLASS Solid_Earth_Tides
% =========================================================================
%
% DESCRIPTION:
%
% This class is a direct port of the International Earth Rotation and
% Reference Systems Service (IERS) Conventions software collection.
%
% Ported in 2024-01-06 from:
%   https://iers-conventions.obspm.fr/content/chapter7/software/dehanttideinel/
%
% The class is used to compute the station tidal displacements
% caused by lunar and solar gravitational attraction (see References). 
%
% NOTE:
%  This class is a library of Static functions, to execute dehantTideInelIERS
%   
% REFERENCES:
%
%    Groten, E., 2000, Geodesists Handbook 2000, Part 4,
%    http://www.gfy.ku.dk/~iag/HB2000/part4/groten.htm. See also
%    ''Parameters of Common Relevance of Astronomy, Geodesy, and
%    Geodynamics," J. Geod., 74, pp. 134-140
%
%    Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
%    displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
%
%    Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
%    IERS Technical Note No. 36, BKG (2010)
%    
%    Pitjeva, E. and Standish, E. M., 2009, ''Proposals for the masses
%    of the three largest asteroids, the Moon-Earth mass ratio and the
%    Astronomical Unit," Celest. Mech. Dyn. Astr., 103, pp. 365-372
%
%    Ries, J. C., Eanes, R. J., Shum, C. K. and Watkins, M. M., 1992,
%    ''Progress in the Determination of the Gravitational Coefficient
%    of the Earth," Geophys. Res. Lett., 19(6), pp. 529-531

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0
%
%--------------------------------------------------------------------------
%  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
%  Written by:        Andrea Gatti
%  Contributors:      Andrea Gatti ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

classdef Solid_Earth_Tides
    properties(Constant)
        % Leap seconds data: [Year, Month, Delta(AT)]
        leap_seconds = [
            1972, 1, 10; 1972, 7, 11; 1973, 1, 12; 1974, 1, 13; 1975, 1, 14;
            1976, 1, 15; 1977, 1, 16; 1978, 1, 17; 1979, 1, 18; 1980, 1, 19;
            1981, 7, 20; 1982, 7, 21; 1983, 7, 22; 1985, 7, 23; 1988, 1, 24;
            1990, 1, 25; 1991, 1, 26; 1992, 7, 27; 1993, 7, 28; 1994, 7, 29;
            1996, 1, 30; 1997, 7, 31; 1999, 1, 32; 2006, 1, 33; 2009, 1, 34;
            2012, 7, 35; 2015, 7, 36; 2017, 1, 37
        ];
    end
    
    methods(Static)
        function dx_tide = dehantTideInelIERS(x_sta, time, x_sun, x_moon)
            % Computes the station tidal displacement caused by lunar and solar gravitational attraction
            %
            % This function computes the station tidal displacement caused by lunar and solar gravitational
            % attraction. It includes corrections for out-of-phase parts of Love numbers for the diurnal and
            % semi-diurnal band, and corrections for the latitude dependence of Love numbers.
            %
            % INPUTS:
            %   x_sta   - Geocentric position of the IGS station [X, Y, Z] in meters
            %   time    - GPS_Time
            %   x_sun   - Geocentric position of the Sun [X, Y, Z] in meters
            %   x_mon   - Geocentric position of the Moon [X, Y, Z] in meters
            %
            % OUTPUTS:
            %   dxtide - Displacement vector [dX, dY, dZ] in meters
            %
            % Example usage:
            %   dxtide = dehantTideInelIERS(x_sta, time, x_sun, x_moon)

            dx_tide = nan(3, time.length());
            date6 = time.get6ColDate;
            yr = date6(:,1);
            month = date6(:,2);
            day = date6(:,3);
            fhr = rem(time.getMatlabTime,1)*24;
            for t = 1:time.length()                
                dx_tide(:,t) = Solid_Earth_Tides.dehantTideInel(x_sta(t,:)', yr(t), month(t), day(t), fhr(t), x_sun(t,:)', x_moon(t,:)');
            end
        end

        function dxtide = dehantTideInel(xsta, yr, month, day, fhr, xsun, xmon)
            % Computes the station tidal displacement caused by lunar and solar gravitational attraction
            %
            % This function computes the station tidal displacement caused by lunar and solar gravitational
            % attraction. It includes corrections for out-of-phase parts of Love numbers for the diurnal and
            % semi-diurnal band, and corrections for the latitude dependence of Love numbers.
            %
            % INPUTS:
            %   xsta   - Geocentric position of the IGS station [X, Y, Z] in meters
            %   yr     - Year
            %   month  - Month
            %   day    - Day of the month
            %   fhr    - Hour in the day
            %   xsun   - Geocentric position of the Sun [X, Y, Z] in meters
            %   xmon   - Geocentric position of the Moon [X, Y, Z] in meters
            %
            % OUTPUTS:
            %   dxtide - Displacement vector [dX, dY, dZ] in meters
            %
            % Example usage:
            %   dxtide = dehantTideInel([4075578.385, 931852.890, 4801570.154], 2009, 4, 13, 0, [137859926952.015, 54228127881.4350, 23509422341.6960], [-179996231.920342, -312468450.131567, -169288918.592160])

            % Constants
            h20 = 0.6078; l20 = 0.0847; h3 = 0.292; l3 = 0.015;

            % Scalar product of station vector with Sun/Moon vector
            [scs, rsta, rsun] = Solid_Earth_Tides.sprod(xsta, xsun);
            [scm, ~, rmon] = Solid_Earth_Tides.sprod(xsta, xmon);
            scsun = scs / rsta / rsun;
            scmon = scm / rsta / rmon;

            % Computation of new h2 and l2
            cosphi = sqrt(xsta(1)^2 + xsta(2)^2) / rsta;
            h2 = h20 - 0.0006 * (1 - 3/2 * cosphi^2);
            l2 = l20 + 0.0002 * (1 - 3/2 * cosphi^2);

            % P2 and P3 terms
            p2sun = 3 * (h2/2 - l2) * scsun^2 - h2/2;
            p2mon = 3 * (h2/2 - l2) * scmon^2 - h2/2;
            p3sun = 5/2 * (h3 - 3 * l3) * scsun^3 + 3/2 * (l3 - h3) * scsun;
            p3mon = 5/2 * (h3 - 3 * l3) * scmon^3 + 3/2 * (l3 - h3) * scmon;

            % Term in direction of Sun/Moon vector
            x2sun = 3 * l2 * scsun;
            x2mon = 3 * l2 * scmon;
            x3sun = 3 * l3 / 2 * (5 * scsun^2 - 1);
            x3mon = 3 * l3 / 2 * (5 * scmon^2 - 1);

            % Factors for Sun/Moon
            mass_ratio_sun = 332946.0482;
            mass_ratio_moon = 0.0123000371;
            re = 6378136.6;
            fac2sun = mass_ratio_sun * re * (re / rsun)^3;
            fac2mon = mass_ratio_moon * re * (re / rmon)^3;
            fac3sun = fac2sun * (re / rsun);
            fac3mon = fac2mon * (re / rmon);

            % Total displacement
            dxtide = zeros(3,1);
            for i = 1:3
                dxtide(i) = fac2sun * (x2sun * xsun(i) / rsun + p2sun * xsta(i) / rsta) + ...
                    fac2mon * (x2mon * xmon(i) / rmon + p2mon * xsta(i) / rsta) + ...
                    fac3sun * (x3sun * xsun(i) / rsun + p3sun * xsta(i) / rsta) + ...
                    fac3mon * (x3mon * xmon(i) / rmon + p3mon * xsta(i) / rsta);
            end

            % Corrections for the out-of-phase part of Love numbers
            xcorsta = Solid_Earth_Tides.st1idiu(xsta, xsun, xmon, fac2sun, fac2mon);
            dxtide = dxtide + xcorsta;

            xcorsta = Solid_Earth_Tides.st1isem(xsta, xsun, xmon, fac2sun, fac2mon);
            dxtide = dxtide + xcorsta;

            xcorsta = Solid_Earth_Tides.st1l1(xsta, xsun, xmon, fac2sun, fac2mon);
            dxtide = dxtide + xcorsta;

            % Step 2 corrections
            [jjm0, jjm1, ~] = Solid_Earth_Tides.cal2jd(yr, month, day);
            fhrd = fhr / 24;
            t = ((jjm0 - 2451545) + jjm1 + fhrd) / 36525;

            [dtt, ~] = Solid_Earth_Tides.dat(yr, month, day, fhr);
            dtt = dtt + 32.184; % Conversion of T in TT time
            t = t + dtt / (3600 * 24 * 36525);

            % The biggest correction effect is here:
            xcorsta = Solid_Earth_Tides.step2diu(xsta, fhr, t);
            dxtide = dxtide + xcorsta;

            xcorsta = Solid_Earth_Tides.step2lon(xsta, t);
            dxtide = dxtide + xcorsta;

        end
       
        function [djm0, djm, j] = cal2jd(iy, im, id)
            % Converts Gregorian calendar date to Modified Julian Date (MJD).
            % Given:
            %   iy, im, id - Year, Month, Day in the Gregorian calendar
            % Returned:
            %   djm0 - MJD zero-point: always 2400000.5
            %   djm - Modified Julian Date for 0 hrs
            %   j - status: 0 = OK, -1 = bad year, -2 = bad month, -3 = bad day

            % Earliest year allowed (4800BC)
            iymin = -4799;

            % Preset status
            j = 0;

            % Validate year
            if iy < iymin
                j = -1;
            else
                % Validate month
                if im >= 1 && im <= 12
                    % Month lengths in days
                    mtab = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

                    % Days in current month
                    ndays = mtab(im);

                    % Allow for leap year
                    if im == 2
                        if mod(iy, 4) == 0
                            ndays = 29;
                        end
                        if mod(iy, 100) == 0 && mod(iy, 400) ~= 0
                            ndays = 28;
                        end
                    end

                    % Validate day
                    if id < 1 || id > ndays
                        j = -3;
                    end

                    % Result
                    my = floor((im - 14) / 12);
                    iypmy = iy + my;
                    djm0 = 2400000.5;

                    %julian day
                    jd = floor(365.25*(iy+4716)) + floor(30.6001*(im+1)) + id - 1537.5;

                    %modified julian day
                    djm = jd - djm0;
                else
                    % Bad month
                    j = -2;
                end
            end
        end

        function [deltat, j] = dat(iy, im, id, fd)
            % Initialize the result to zero and the status to OK
            deltat = 0;
            j = 0;
            
            % If invalid fraction of a day, set error status and give up.
            if fd < 0 || fd > 1
                j = -4;
                return;
            end
            
            % Convert the date into an MJD using the cal2jd method.
            [djm0, djm, j] = Solid_Earth_Tides.cal2jd(iy, im, id);
            
            % If invalid date, give up.
            if j < 0
                return;
            end
            
            % Combine year and month for comparison.
            m = 12*iy + im;
            
            % Find the most recent table entry.
            is = 0;
            for n = size(Solid_Earth_Tides.leap_seconds, 1):-1:1
                leap_year_month = 12*Solid_Earth_Tides.leap_seconds(n, 1) + Solid_Earth_Tides.leap_seconds(n, 2);
                if m >= leap_year_month
                    is = n;
                    break;
                end
            end
            
            % Get the Delta(AT) if an entry was found.
            if is > 0
                deltat = Solid_Earth_Tides.leap_seconds(is, 3);
            else
                % No applicable leap second found; use default or handle as needed.
                j = -5;
            end
        end

        function norm8 = norm8(a)
            % Normalizes a given three-dimensional vector A
            % Input:
            %   a - A 3-element vector
            % Output:
            %   norm8 - The norm (magnitude) of the vector A
            
            % Validate input vector
            if length(a) ~= 3
                error('Input must be a 3-element vector.');
            end
            
            % Calculate the norm of the vector
            norm8 = sqrt(a(1)^2 + a(2)^2 + a(3)^2);
        end

        function [scal, r1, r2] = sprod(x, y)
            % Computes the scalar product of two vectors and their norms
            % Input:
            %   x - Components of vector x (3-element vector)
            %   y - Components of vector y (3-element vector)
            % Output:
            %   scal - Scalar product of vector x and vector y
            %   r1 - Length (norm) of vector x
            %   r2 - Length (norm) of vector y
            
            % Validate input vectors
            if length(x) ~= 3 || length(y) ~= 3
                error('Both input vectors must be 3-element vectors.');
            end
            
            % Compute the norms of the vectors
            r1 = sqrt(x(1)^2 + x(2)^2 + x(3)^2);
            r2 = sqrt(y(1)^2 + y(2)^2 + y(3)^2);
            
            % Compute the scalar product
            scal = x(1)*y(1) + x(2)*y(2) + x(3)*y(3);
        end

        function xcorsta = st1idiu(xsta, xsun, xmon, fac2sun, fac2mon)
            % This function gives the out-of-phase corrections induced by
            % mantle anelasticity in the diurnal band.
            %
            % Input:
            %   xsta - Geocentric position of the IGS station (3-element vector)
            %   xsun - Geocentric position of the Sun (3-element vector)
            %   xmon - Geocentric position of the Moon (3-element vector)
            %   fac2sun - Degree 2 TGP factor for the Sun
            %   fac2mon - Degree 2 TGP factor for the Moon
            %
            % Output:
            %   xcorsta - Out of phase station corrections for diurnal band (3-element vector)
            
            % Constants
            dhi = -0.0025;
            dli = -0.0007;
            
            % Compute the normalized position vector of the IGS station.
            rsta = Solid_Earth_Tides.norm8(xsta);
            sinphi = xsta(3) / rsta;
            cosphi = sqrt(xsta(1)^2 + xsta(2)^2) / rsta;
            cos2phi = cosphi^2 - sinphi^2;
            sinla = xsta(2) / (cosphi * rsta);
            cosla = xsta(1) / (cosphi * rsta);
            
            % Compute the normalized position vector of the Moon and the Sun.
            rmon = Solid_Earth_Tides.norm8(xmon);
            rsun = Solid_Earth_Tides.norm8(xsun);
            
            % Compute corrections.
            drsun = -3 * dhi * sinphi * cosphi * fac2sun * xsun(3) * (xsun(1) * sinla - xsun(2) * cosla) / rsun^2;
            drmon = -3 * dhi * sinphi * cosphi * fac2mon * xmon(3) * (xmon(1) * sinla - xmon(2) * cosla) / rmon^2;
            dnsun = -3 * dli * cos2phi * fac2sun * xsun(3) * (xsun(1) * sinla - xsun(2) * cosla) / rsun^2;
            dnmon = -3 * dli * cos2phi * fac2mon * xmon(3) * (xmon(1) * sinla - xmon(2) * cosla) / rmon^2;
            desun = -3 * dli * sinphi * fac2sun * xsun(3) * (xsun(1) * cosla + xsun(2) * sinla) / rsun^2;
            demon = -3 * dli * sinphi * fac2mon * xmon(3) * (xmon(1) * cosla + xmon(2) * sinla) / rmon^2;
            
            dr = drsun + drmon;
            dn = dnsun + dnmon;
            de = desun + demon;
            
            % Compute the corrections for the station.
            xcorsta = zeros(3,1);
            xcorsta(1) = dr * cosla * cosphi - de * sinla - dn * sinphi * cosla;
            xcorsta(2) = dr * sinla * cosphi + de * cosla - dn * sinphi * sinla;
            xcorsta(3) = dr * sinphi + dn * cosphi;
        end        

        function xcorsta = st1isem(xsta, xsun, xmon, fac2sun, fac2mon)
            % This function gives the out-of-phase corrections induced by
            % mantle anelasticity in the semi-diurnal band.
            %
            % Input:
            %   xsta - Geocentric position of the IGS station (3-element vector)
            %   xsun - Geocentric position of the Sun (3-element vector)
            %   xmon - Geocentric position of the Moon (3-element vector)
            %   fac2sun - Degree 2 TGP factor for the Sun
            %   fac2mon - Degree 2 TGP factor for the Moon
            %
            % Output:
            %   xcorsta - Out of phase station corrections for semi-diurnal band (3-element vector)

            % Constants
            dhi = -0.0022;
            dli = -0.0007;

            % Compute the normalized position vector of the IGS station.
            rsta = Solid_Earth_Tides.norm8(xsta);
            sinphi = xsta(3) / rsta;
            cosphi = sqrt(xsta(1)^2 + xsta(2)^2) / rsta;
            sinla = xsta(2) / (cosphi * rsta);
            cosla = xsta(1) / (cosphi * rsta);
            costwola = cosla^2 - sinla^2;
            sintwola = 2 * cosla * sinla;

            % Compute the normalized position vector of the Moon and the Sun.
            rmon = Solid_Earth_Tides.norm8(xmon);
            rsun = Solid_Earth_Tides.norm8(xsun);

            % Compute corrections.
            drsun = -3/4 * dhi * cosphi^2 * fac2sun * ((xsun(1)^2 - xsun(2)^2) * ...
                sintwola - 2 * xsun(1) * xsun(2) * costwola) / rsun^2;
            drmon = -3/4 * dhi * cosphi^2 * fac2mon * ((xmon(1)^2 - xmon(2)^2) * ...
                sintwola - 2 * xmon(1) * xmon(2) * costwola) / rmon^2;

            dnsun = 3/2 * dli * sinphi * cosphi * fac2sun * ((xsun(1)^2 - xsun(2)^2) * ...
                sintwola - 2 * xsun(1) * xsun(2) * costwola) / rsun^2;
            dnmon = 3/2 * dli * sinphi * cosphi * fac2mon * ((xmon(1)^2 - xmon(2)^2) * ...
                sintwola - 2 * xmon(1) * xmon(2) * costwola) / rmon^2;

            desun = -3/2 * dli * cosphi * fac2sun * ((xsun(1)^2 - xsun(2)^2) * ...
                costwola + 2 * xsun(1) * xsun(2) * sintwola) / rsun^2;
            demon = -3/2 * dli * cosphi * fac2mon * ((xmon(1)^2 - xmon(2)^2) * ...
                costwola + 2 * xmon(1) * xmon(2) * sintwola) / rmon^2;

            dr = drsun + drmon;
            dn = dnsun + dnmon;
            de = desun + demon;

            % Compute the corrections for the station.
            xcorsta = zeros(3,1);
            xcorsta(1) = dr * cosla * cosphi - de * sinla - dn * sinphi * cosla;
            xcorsta(2) = dr * sinla * cosphi + de * cosla - dn * sinphi * sinla;
            xcorsta(3) = dr * sinphi + dn * cosphi;
        end

        function xcorsta = st1l1(xsta, xsun, xmon, fac2sun, fac2mon)
            % This function computes corrections induced by the latitude
            % dependence given by L^1 in the diurnal and semi-diurnal band.
            %
            % Input:
            %   xsta - Geocentric position of the IGS station (3-element vector)
            %   xsun - Geocentric position of the Sun (3-element vector)
            %   xmon - Geocentric position of the Moon (3-element vector)
            %   fac2sun - Degree 2 TGP factor for the Sun
            %   fac2mon - Degree 2 TGP factor for the Moon
            %
            % Output:
            %   xcorsta - Station corrections (3-element vector)

            % Constants
            l1d = 0.0012;
            l1sd = 0.0024;

            % Normalized position vectors and related calculations
            rsta = Solid_Earth_Tides.norm8(xsta);
            sinphi = xsta(3) / rsta;
            cosphi = sqrt(xsta(1)^2 + xsta(2)^2) / rsta;
            sinla = xsta(2) / (cosphi * rsta);
            cosla = xsta(1) / (cosphi * rsta);

            rmon = Solid_Earth_Tides.norm8(xmon);
            rsun = Solid_Earth_Tides.norm8(xsun);

            % Corrections for diurnal band
            l1 = l1d;
            dnsun = -l1 * sinphi^2 * fac2sun * xsun(3) * (xsun(1) * cosla + xsun(2) * sinla) / rsun^2;
            dnmon = -l1 * sinphi^2 * fac2mon * xmon(3) * (xmon(1) * cosla + xmon(2) * sinla) / rmon^2;
            desun = l1 * sinphi * (cosphi^2 - sinphi^2) * fac2sun * xsun(3) * ...
                    (xsun(1) * sinla - xsun(2) * cosla) / rsun^2;
            demon = l1 * sinphi * (cosphi^2 - sinphi^2) * fac2mon * xmon(3) * ...
                    (xmon(1) * sinla - xmon(2) * cosla) / rmon^2;

            de = 3 * (desun + demon);
            dn = 3 * (dnsun + dnmon);

            xcorsta = [-de * sinla - dn * sinphi * cosla; 
                       de * cosla - dn * sinphi * sinla;
                       dn * cosphi];

            % Corrections for semi-diurnal band
            l1 = l1sd;
            costwola = cosla^2 - sinla^2;
            sintwola = 2 * cosla * sinla;

            dnsun = -l1 / 2 * sinphi * cosphi * fac2sun * ...
                    ((xsun(1)^2 - xsun(2)^2) * costwola + 2 * xsun(1) * xsun(2) * sintwola) / rsun^2;
            dnmon = -l1 / 2 * sinphi * cosphi * fac2mon * ...
                    ((xmon(1)^2 - xmon(2)^2) * costwola + 2 * xmon(1) * xmon(2) * sintwola) / rmon^2;
            desun = -l1 / 2 * sinphi^2 * cosphi * fac2sun * ...
                    ((xsun(1)^2 - xsun(2)^2) * sintwola - 2 * xsun(1) * xsun(2) * costwola) / rsun^2;
            demon = -l1 / 2 * sinphi^2 * cosphi * fac2mon * ...
                    ((xmon(1)^2 - xmon(2)^2) * sintwola - 2 * xmon(1) * xmon(2) * costwola) / rmon^2;

            de = 3 * (desun + demon);
            dn = 3 * (dnsun + dnmon);

            % Combine corrections
            xcorsta = xcorsta + [-de * sinla - dn * sinphi * cosla; 
                                 de * cosla - dn * sinphi * sinla;
                                 dn * cosphi];
        end

        function xcorsta = step2diu(xsta, fhr, t)
            % Constants
            d2pi = 6.283185307179586476925287;
            deg2rad = d2pi / 360.0;
            
            % DATDI array initialization
            datdi = [
                -3, 0, 2, 0, 0, -0.01, 0, 0, 0;
                -3, 2, 0, 0, 0, -0.01, 0, 0, 0;
                -2, 0, 1, -1, 0, -0.02, 0, 0, 0;
                -2, 0, 1, 0, 0, -0.08, 0, -0.01, 0.01;
                -2, 2, -1, 0, 0, -0.02, 0, 0, 0;
                -1, 0, 0, -1, 0, -0.10, 0, 0, 0;
                -1, 0, 0, 0, 0, -0.51, 0, -0.02, 0.03;
                -1, 2, 0, 0, 0, 0.01, 0, 0, 0;
                0, -2, 1, 0, 0, 0.01, 0, 0, 0;
                0, 0, -1, 0, 0, 0.02, 0, 0, 0;
                0, 0, 1, 0, 0, 0.06, 0, 0, 0;
                0, 0, 1, 1, 0, 0.01, 0, 0, 0;
                0, 2, -1, 0, 0, 0.01, 0, 0, 0;
                1, -3, 0, 0, 1, -0.06, 0, 0, 0;
                1, -2, 0, -1, 0, 0.01, 0, 0, 0;
                1, -2, 0, 0, 0, -1.23, -0.07, 0.06, 0.01;
                1, -1, 0, 0, -1, 0.02, 0, 0, 0;
                1, -1, 0, 0, 1, 0.04, 0, 0, 0;
                1, 0, 0, -1, 0, -0.22, 0.01, 0.01, 0;
                1, 0, 0, 0, 0, 12.00, -0.80, -0.67, -0.03;
                1, 0, 0, 1, 0, 1.73, -0.12, -0.10, 0;
                1, 0, 0, 2, 0, -0.04, 0, 0, 0;
                1, 1, 0, 0, -1, -0.50, -0.01, 0.03, 0;
                1, 1, 0, 0, 1, 0.01, 0, 0, 0;
                0, 1, 0, 1, -1, -0.01, 0, 0, 0;
                1, 2, -2, 0, 0, -0.01, 0, 0, 0;
                1, 2, 0, 0, 0, -0.11, 0.01, 0.01, 0;
                2, -2, 1, 0, 0, -0.01, 0, 0, 0;
                2, 0, -1, 0, 0, -0.02, 0, 0, 0;
                3, 0, 0, 0, 0, 0, 0, 0, 0;
                3, 0, 0, 1, 0, 0, 0, 0, 0
                ];
            % Phase angles in degrees
            s = 218.31664563 + (481267.88194 + (-0.0014663889 + (0.00000185139) * t) * t) * t;
            tau = fhr * 15.0 + 280.4606184 + (36000.7700536 + (0.00038793 + (-0.0000000258) * t) * t) * t - s;
            pr = (1.396971278 + (0.000308889 + (0.000000021 + (0.000000007) * t) * t) * t) * t;
            s = s + pr;

            h = 280.46645 + (36000.7697489 + (0.00030322222 + (0.000000020 - 0.00000000654 * t) * t) * t) * t;
            p = 83.35324312 + (4069.01363525 + (-0.01032172222 - 0.0000124991 + (0.00000005263) * t) * t) * t;
            zns = 234.95544499 + (1934.13626197 + (-0.00207561111 - 0.00000213944 + (0.00000001650) * t) * t) * t;
            ps = 282.93734098 + (1.71945766667 + (0.00045688889 - 0.00000001778 - 0.00000000334 * t) * t) * t;
            
            rsta = sqrt(sum(xsta.^2));
            sinphi = xsta(3) / rsta;
            cosphi = sqrt(xsta(1)^2 + xsta(2)^2) / rsta;
            cosla = xsta(1) / (cosphi * rsta);
            sinla = xsta(2) / (cosphi * rsta);
            zla = atan2(xsta(2), xsta(1));

            % Initialize correction vector
            xcorsta = zeros(3, 1);

            % Loop through the DATDI matrix rows to calculate corrections
            for j = 1:size(datdi, 1)
                thetaf = (tau + datdi(j, 1)*s + datdi(j, 2)*h + datdi(j, 3)*p + datdi(j, 4)*zns + datdi(j, 5)*ps) * deg2rad;

                dr = datdi(j, 6) * 2.0 * sinphi * cosphi * sin(thetaf + zla) + datdi(j, 7) * 2.0 * sinphi * cosphi * cos(thetaf + zla);
                dn = datdi(j, 8) * (cosphi^2 - sinphi^2) * sin(thetaf + zla) + datdi(j, 9) * (cosphi^2 - sinphi^2) * cos(thetaf + zla);
                de = datdi(j, 8) * sinphi * cos(thetaf + zla) - datdi(j, 9) * sinphi * sin(thetaf + zla);

                xcorsta(1) = xcorsta(1) + dr * cosla * cosphi - de * sinla - dn * sinphi * cosla;
                xcorsta(2) = xcorsta(2) + dr * sinla * cosphi + de * cosla - dn * sinphi * sinla;
                xcorsta(3) = xcorsta(3) + dr * sinphi + dn * cosphi;
            end

            % Convert from mm to meters
            xcorsta = xcorsta / 1000.0;
        end

        function xcorsta = step2lon(xsta, t)
            % Constants
            d2pi = 6.283185307179586476925287;
            deg2rad = d2pi / 360.0;
            
            % DATDI array
            datdi = [0, 0, 0, 1, 0, 0.47, 0.23, 0.16, 0.07;
                     0, 2, 0, 0, 0, -0.20, -0.12, -0.11, -0.05;
                     1, 0, -1, 0, 0, -0.11, -0.08, -0.09, -0.04;
                     2, 0, 0, 0, 0, -0.13, -0.11, -0.15, -0.07;
                     2, 0, 0, 1, 0, -0.05, -0.05, -0.06, -0.03];
            
            % Compute the phase angles in degrees
            s = 218.31664563 + (481267.88194 + (-0.0014663889 + (0.00000185139)*t)*t)*t;
            pr = (1.396971278 + (0.000308889 + (0.000000021 + (0.000000007)*t)*t)*t)*t;
            s = s + pr;
            h = 280.46645 + (36000.7697489 + (0.00030322222 + (0.000000020 - 0.00000000654*t)*t)*t)*t;
            p = 83.35324312 + (4069.01363525 + (-0.01032172222 - 0.0000124991 + (0.00000005263)*t)*t)*t;
            zns = 234.95544499 + (1934.13626197 + (-0.00207561111 - 0.00000213944 + (0.00000001650)*t)*t)*t;
            ps = 282.93734098 + (1.71945766667 + (0.00045688889 - 0.00000001778 - 0.00000000334*t)*t)*t;
            
            % Reduce angles to between the range 0 and 360
            s = mod(s, 360.0);
            h = mod(h, 360.0);
            p = mod(p, 360.0);
            zns = mod(zns, 360.0);
            ps = mod(ps, 360.0);
            
            rsta = sqrt(sum(xsta.^2));
            sinphi = xsta(3) / rsta;
            cosphi = sqrt(xsta(1)^2 + xsta(2)^2) / rsta;
            cosla = xsta(1) / (cosphi * rsta);
            sinla = xsta(2) / (cosphi * rsta);
            
            % Initialize correction vector
            xcorsta = zeros(3, 1);
            dr_tot = 0;
            dn_tot = 0;
            
            % Loop through the DATDI matrix rows to calculate corrections
            for j = 1:5
                thetaf = (datdi(j, 1)*s + datdi(j, 2)*h + datdi(j, 3)*p + datdi(j, 4)*zns + datdi(j, 5)*ps) * deg2rad;
                
                dr = datdi(j, 6)*(3*sinphi^2 - 1)/2*cos(thetaf) + datdi(j, 8)*(3*sinphi^2 - 1)/2*sin(thetaf);
                dn = datdi(j, 7)*(cosphi*sinphi*2)*cos(thetaf) + datdi(j, 9)*(cosphi*sinphi*2)*sin(thetaf);
                
                dr_tot = dr_tot + dr;
                dn_tot = dn_tot + dn;
                
                xcorsta(1) = xcorsta(1) + dr*cosla*cosphi - dn*sinphi*cosla;
                xcorsta(2) = xcorsta(2) + dr*sinla*cosphi - dn*sinphi*sinla;
                xcorsta(3) = xcorsta(3) + dr*sinphi + dn*cosphi;
            end
            
            % Convert from mm to meters
            xcorsta = xcorsta / 1000.0;
        end

    end
    
    % Tester functions -----------------------------------------------------------------------------

    methods(Static)
        function test_dehanttideinel()
            % Define test cases
            test_cases = {
                % Test case 1
                struct('xsta', [4075578.385; 931852.890; 4801570.154], ...
                'xsun', [137859926952.015; 54228127881.4350; 23509422341.6960], ...
                'xmon', [-179996231.920342; -312468450.131567; -169288918.592160], ...
                'yr', 2009, 'month', 4, 'day', 13, 'fhr', 0.00, ...
                'expected', [0.7700420357108125891e-01; 0.6304056321824967613e-01; 0.5516568152597246810e-01]),
                % Test case 2
                struct('xsta', [1112189.660; -4842955.026; 3985352.284], ...
                'xsun', [-54537460436.2357; 130244288385.279; 56463429031.5996], ...
                'xmon', [300396716.912; 243238281.451; 120548075.939], ...
                'yr', 2012, 'month', 7, 'day', 13, 'fhr', 0.00, ...
                'expected', [-0.2036831479592075833e-01; 0.5658254776225972449e-01; -0.7597679676871742227e-01]),
                % Test case 3
                struct('xsta', [1112200.5696; -4842957.8511; 3985345.9122], ...
                'xsun', [100210282451.6279; 103055630398.3160; 56855096480.4475], ...
                'xmon', [369817604.4348; 1897917.5258; 120804980.8284], ...
                'yr', 2015, 'month', 7, 'day', 15, 'fhr', 0.00, ...
                'expected', [0.00509570869172363845; 0.0828663025983528700; -0.0636634925404189617]),
                % Test case 4
                struct('xsta', [1112152.8166; -4842857.5435; 3985496.1783], ...
                'xsun', [8382471154.1312895; 10512408445.356153; -5360583240.3763866], ...
                'xmon', [380934092.93550891; 2871428.1904491195; 79015680.553570181], ...
                'yr', 2017, 'month', 1, 'day', 15, 'fhr', 0.00, ...
                'expected', [-18.217357581922339; -23.505348376537949; 12.097611382175685]) %[0.0050957086917236384; 0.082866302598352870; -0.063663492540418962])
                % Test case 5
                struct('xsta', [1112152.8166; -4842857.5435; 3985496.1783], ...
                'xsun', [8382471154.1312895; 10512408445.356153; -5360583240.3763866], ...
                'xmon', [380934092.93550891; 2871428.1904491195; 79015680.553570181], ...
                'yr', 2017, 'month', 1, 'day', 15, 'fhr', 0.00, ...
                'expected', [0.0050957086917236384; 0.082866302598352870; -0.063663492540418962])
                };

            % Set tolerance for comparison
            tolerance = 1e-6;

            % Run test cases
            for i = 1:length(test_cases)
                tc = test_cases{i};
                dxtide = Solid_Earth_Tides.dehantTideInel(tc.xsta, tc.yr, tc.month, tc.day, tc.fhr, tc.xsun, tc.xmon);

                % Check if the result is within the tolerance
                if all(abs(dxtide - tc.expected) < tolerance)
                    fprintf('Test case %d passed.\n', i);
                else
                    fprintf('Test case %d failed.\n', i);
                    disp('Computed output:');
                    disp(dxtide);
                    disp('Expected output:');
                    disp(tc.expected);
                end
            end
        end

        function test_iauDat()
            % Test case date: December 31, 2016, just before the leap second introduction
            iy = 2016; im = 12; id = 31; fd = 0.5;  % Halfway through the day

            % Expected Delta T for this date (leap seconds up to this date - 37 seconds)
            expected_deltat = 36.0;

            % Call the function
            [deltat, j] = Solid_Earth_Tides.dat(iy, im, id, fd);

            % Check results
            tolerance = 1e-6;  % Tolerance for numerical comparison
            if abs(deltat - expected_deltat) <= tolerance && j == 0
                disp('Test passed successfully.');
                fprintf('Delta T (TAI - UTC) on %d-%02d-%02d is %.1f seconds as expected.\n', iy, im, id, deltat);
            else
                disp('Test failed.');
                fprintf('Expected Delta T: %.1f seconds, Computed Delta T: %.1f seconds, Status: %d\n', expected_deltat, deltat, j);
            end
        end

        function test_norm8()
            % Given input vector A
            a = [2; 2; 1]; % 3-element vector

            % Expected output from the test case
            expected_output = 3; % The norm (magnitude) of vector A

            % Call the norm8 function
            norm8_result = Solid_Earth_Tides.norm8(a);

            % Display the results
            fprintf('Computed norm: %f\n', norm8_result);
            fprintf('Expected norm: %f\n', expected_output);

            % Tolerance for comparison
            tolerance = 1e-12;

            % Check if the computed result is within the tolerance of the expected result
            if abs(norm8_result - expected_output) < tolerance
                disp('Test passed successfully.');
            else
                disp('Test failed.');
            end
        end

        function test_sprod()
            % Given input vectors
            x = [2; 2; 3];
            y = [1; 3; 4];

            % Expected output
            expected_scal = 20;
            expected_r1 = 4.123105625617660586;
            expected_r2 = 5.099019513592784492;

            % Tolerance for comparison
            tolerance = 1e-12;

            % Call the sprod function
            [scal, r1, r2] = Solid_Earth_Tides.sprod(x, y);

            % Check results
            scal_diff = abs(scal - expected_scal);
            r1_diff = abs(r1 - expected_r1);
            r2_diff = abs(r2 - expected_r2);

            % Display results
            if scal_diff < tolerance && r1_diff < tolerance && r2_diff < tolerance
                disp('Test passed successfully.');
            else
                disp('Test failed.');
                fprintf('Computed scalar product: %f, Expected: %f\n', scal, expected_scal);
                fprintf('Computed norm of x: %f, Expected: %f\n', r1, expected_r1);
                fprintf('Computed norm of y: %f, Expected: %f\n', r2, expected_r2);
            end
        end

        function test_st1idiu()
            % Given input parameters
            xsta = [4075578.385; 931852.890; 4801570.154]; % in meters
            xsun = [137859926952.015; 54228127881.4350; 23509422341.6960]; % in meters
            xmon = [-179996231.920342; -312468450.131567; -169288918.592160]; % in meters
            fac2sun = 0.163271964478954; % Degree 2 TGP factor for the Sun in 1/meters
            fac2mon = 0.321989090026845; % Degree 2 TGP factor for the Moon in 1/meters

            % Expected output from the test case in meters
            expected_output = [-0.2836337012840008001e-03;
                0.1125342324347507444e-03;
                -0.2471186224343683169e-03];

            % Tolerance for comparison
            tolerance = 1e-12;

            % Call the implemented function
            xcorsta = Solid_Earth_Tides.st1idiu(xsta, xsun, xmon, fac2sun, fac2mon);

            % Compute the difference between expected and computed results
            diff = abs(xcorsta - expected_output);

            % Check if the difference is within the tolerance
            if all(diff < tolerance)
                disp('Test passed successfully.');
            else
                disp('Test failed.');
                disp('Computed output:');
                disp(xcorsta);
                disp('Expected output:');
                disp(expected_output);
            end
        end

        function test_st1isem()
            % Given input parameters
            xsta = [4075578.385; 931852.890; 4801570.154]; % in meters
            xsun = [137859926952.015; 54228127881.4350; 23509422341.6960]; % in meters
            xmon = [-179996231.920342; -312468450.131567; -169288918.592160]; % in meters
            fac2sun = 0.163271964478954; % Degree 2 TGP factor for the Sun in 1/meters
            fac2mon = 0.321989090026845; % Degree 2 TGP factor for the Moon in 1/meters

            % Expected output from the test case in meters
            expected_output = [-0.2801334805106874015e-03;
                0.2939522229284325029e-04;
                -0.6051677912316721561e-04];

            % Tolerance for comparison
            tolerance = 1e-12;

            % Call the implemented function
            xcorsta = Solid_Earth_Tides.st1isem(xsta, xsun, xmon, fac2sun, fac2mon);

            % Compute the difference between expected and computed results
            diff = abs(xcorsta - expected_output);

            % Check if the difference is within the tolerance
            if all(diff < tolerance)
                disp('Test passed successfully.');
            else
                disp('Test failed.');
                disp('Computed output:');
                disp(xcorsta);
                disp('Expected output:');
                disp(expected_output);
            end
        end

        function test_st1l1()
            % Given input parameters for the test case
            xsta = [4075578.385; 931852.890; 4801570.154]; % Geocentric position of the IGS station in meters
            xsun = [137859926952.015; 54228127881.4350; 23509422341.6960]; % Geocentric position of the Sun in meters
            xmon = [-179996231.920342; -312468450.131567; -169288918.592160]; % Geocentric position of the Moon in meters
            fac2sun = 0.163271964478954; % Degree 2 TGP factor for the Sun in 1/meters
            fac2mon = 0.321989090026845; % Degree 2 TGP factor for the Moon in 1/meters
            t = 0.1059411362080767; % Julian centuries (not used in st1l1 but listed for completeness)

            % Expected output from the test case in meters
            expected_output = [0.2367189532359759044e-03;
                0.5181609907284959182e-03;
                -0.3014881422940427977e-03];

            % Tolerance for comparison
            tolerance = 1e-12;

            % Call the st1l1 function
            xcorsta = Solid_Earth_Tides.st1l1(xsta, xsun, xmon, fac2sun, fac2mon);

            % Compute the difference between expected and computed results
            diff = abs(xcorsta - expected_output);

            % Check if the difference is within the tolerance
            if all(diff < tolerance)
                disp('Test passed successfully.');
            else
                disp('Test failed.');
                disp('Computed output:');
                disp(xcorsta);
                disp('Expected output:');
                disp(expected_output);
            end
        end

        function test_step2diu()
            % Given input parameters
            xsta = [4075578.385; 931852.890; 4801570.154]; % in meters
            fhr = 0.00; % hours
            t = 0.1059411362080767; % Julian centuries

            % Expected output from test case
            expected_output = [0.4193085327321284701e-02;
                0.1456681241014607395e-02;
                0.5123366597450316508e-02]; % in meters

            % Tolerance for comparison (considering numerical precision)
            tolerance = 1e-12;

            % Call the implemented function
            xcorsta = Solid_Earth_Tides.step2diu(xsta, fhr, t);

            % Compute the difference between expected and computed results
            diff = abs(xcorsta - expected_output);

            % Check if the difference is within the tolerance
            if all(diff < tolerance)
                disp('Test passed successfully.');
            else
                disp('Test failed.');
                disp('Computed output:');
                disp(xcorsta);
                disp('Expected output:');
                disp(expected_output);
            end
        end

        function test_step2lon()
            % Given input parameters
            xsta = [4075578.385; 931852.890; 4801570.154]; % in meters
            t = 0.1059411362080767; % Julian centuries

            % Expected output from test case
            expected_output = [-0.9780962849562107762e-04;
                -0.2236349699932734273e-04;
                0.3561945821351565926e-03]; % in meters

            % Tolerance for comparison (considering numerical precision)
            tolerance = 1e-12;

            % Call the implemented function
            xcorsta = Solid_Earth_Tides.step2lon(xsta, t);

            % Compute the difference between expected and computed results
            diff = abs(xcorsta - expected_output);

            % Check if the difference is within the tolerance
            if all(diff < tolerance)
                disp('Test passed successfully.');
            else
                disp('Test failed.');
                disp('Computed output:');
                disp(xcorsta);
                disp('Expected output:');
                disp(expected_output);
            end
        end
    end
end
