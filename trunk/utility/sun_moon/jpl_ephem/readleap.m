function readleap

% read leap seconds data file

% output via global

%  jdateleap = array of utc julian dates
%  leapsec   = array of leap seconds

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global jdateleap leapsec

leapdata = csvread('tai-utc.dat');

% find length of data arrays

ndata = length(leapdata);

% put julian dates and leap seconds into data arrays

for i = 1:1:ndata
    
    jdateleap(i) = leapdata(i, 1);

    leapsec(i) = leapdata(i, 2);
    
end

