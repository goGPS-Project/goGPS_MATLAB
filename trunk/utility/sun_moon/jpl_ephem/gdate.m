function [month, day, year] = gdate (jdate)

% convert Julian date to Gregorian (calendar) date
   
% input

%  jdate = julian day

% output

%  month = calendar month [1 - 12]
%  day   = calendar day [1 - 31]
%  year  = calendar year [yyyy]

%  note: day may include fractional part

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

jd = jdate;

z = fix(jd + .5);
fday = jd + .5 - z;

if (fday < 0)
   fday = fday + 1;
   z = z - 1;
end

if (z < 2299161)
   a = z;
else
   alpha = floor((z - 1867216.25) / 36524.25);
   a = z + 1 + alpha - floor(alpha / 4);
end
 
b = a + 1524;
c = fix((b - 122.1) / 365.25);
d = fix(365.25 * c);
e = fix((b - d) / 30.6001);
day = b - d - fix(30.6001 * e) + fday;
 
if (e < 14)
   month = e - 1;
else
   month = e - 13;
end
 
if (month > 2)
   year = c - 4716;
else
   year = c - 4715;
end


