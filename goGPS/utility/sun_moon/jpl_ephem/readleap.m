function readleap

% read leap seconds data file

% output via global

%  jdateleap = array of utc julian dates
%  leapsec   = array of leap seconds

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global jdateleap leapsec

% I cannot load this external file when the program is deployed
% Set here manually the leap-seconds as julian date
leapdata = csvread('tai-utc.dat');

% find length of data arrays

ndata = length(leapdata);

% put julian dates and leap seconds into data arrays

jdateleap = leapdata(:, 1);
leapsec = leapdata(:, 2);
