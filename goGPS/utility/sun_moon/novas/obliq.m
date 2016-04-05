function y = obliq(t)

% function to compute mean obliquity of the ecliptic in arcseconds
% capitaine et al. (2003), astronomy and astrophysics 412, 567-586,
% expression from eq. (39) with obliquity at j2000.0 taken from
% eq. (37) or table 8

% input

%  t = tdb julian centuries

% output

%  y = mean obliquity of the ecliptic in arcseconds

% ported from NOVAS 3.0

%%%%%%%%%%%%%%%%%%%%%%%

y = (((( - 0.0000000434d0 * t ...
    - 0.000000576d0) * t ...
    + 0.00200340d0) * t ...
    - 0.0001831d0) * t ...
    - 46.836769d0) * t + 84381.406d0;