function theta = erot (date1, date2)

% this function returns the value of the earth rotation angle
% (theta) for a given ut1 julian date. the expression used is
% taken from the note to iau resolution b1.8 of 2000.

% input

%  date1 = high-order part of ut1 julian date

%  date2 = low-order part of ut1 julian date

% output

%  theta = earth rotation angle in degrees

% ported from NOVAS 3.0

%%%%%%%%%%%%%%%%%%%%%%%

t0 = 2451545.0d0;

% the algorithm used below is equivalent to the cannonical
% theta = 0.7790572732640d0 + 1.00273781191135448d0 * t,
% where t is the time in ut1 days from 2000.0 (t = date1 + date2 - t0),
% but it avoids many two-pi 'wraps' that decrease precision
% (adopted from sofa routine iau_era00 by pat wallace; see also
% expression at top of page 35 of iers conventions (1996))

thet1 = 0.7790572732640d0 + 0.00273781191135448d0 * (date1 - t0);

thet2 = 0.00273781191135448d0 * date2;

thet3 = mod(date1, 1.0d0) + mod(date2, 1.0d0);

theta = mod(thet1 + thet2 + thet3, 1.0d0) * 360.0d0;

if (theta < 0.0d0)
    
    theta = theta + 360.0d0;
    
end


