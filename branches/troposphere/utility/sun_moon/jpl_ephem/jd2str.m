function [cdstr, utstr] = jd2str(jdate)

% convert Julian date to string equivalent
% calendar date and universal time

% input

%  jdate = Julian date

% output

%  cdstr = calendar date string
%  utstr = universal time string

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[month, day, year] = gdate(jdate);

% serial date number

sdn = datenum(year, month, day);

% create calendar date string

cdstr = datestr(sdn, 1);

% create universal time string

utstr = datestr(day - fix(day), 'HH:MM:SS.FFF');
