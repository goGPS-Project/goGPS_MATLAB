function jd = julday(y,m,d,h)

% SYNTAX:
%   jd = julday(y,m,d,h);
%
% INPUT:
%   y = year
%   m = month
%   d = day
%   h = hour
%
% OUTPUT:
%   jd = julian day
%
% DESCRIPTION:
%   Julian day computation.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.0 beta
%
% Copyright (C) Kai Borre
% Kai Borre 02-14-01
%
% Adapted by Mirko Reguzzoni, Eugenio Realini, 2009
%----------------------------------------------------------------------------------------------

if m <= 2
    y = y-1;
    m = m+12;
end

%return julian day
jd = floor(365.25*(y+4716))+floor(30.6001*(m+1))+d+h/24-1537.5;
