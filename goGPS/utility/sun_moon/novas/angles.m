function [ra, dec] = angles (pos)

% this function converts a vector to angular quantities

% input

%  pos = position vector, equatorial rectangular coordinates

% output

%  ra = right ascension in hours 

%  dec = declination in degrees 

% ported from NOVAS 3.1

%%%%%%%%%%%%%%%%%%%%%%%

seccon = 180.0d0 * 3600.0d0 / pi;

xyproj = sqrt(pos(1)^2 + pos(2)^2);

r = 0.0d0;

if (xyproj > 0.0d0)
    
    r = atan2(pos(2), pos(1));
    
end

ra = r * seccon / 54000.0d0;

if (ra < 0.0d0)
    
    ra = ra + 24.0d0;
    
end

if (ra >= 24.0d0)
    
    ra = ra - 24.0d0;
    
end

d = atan2(pos(3), xyproj);

dec = d * seccon / 3600.0d0;


