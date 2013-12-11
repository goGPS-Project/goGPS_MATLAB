function [week, sec_of_week] = gps_time(julday)

% SYNTAX:
%   [week, sec_of_week] = gps_time(julday);
%
% INPUT:
%   julday = julian day
%
% OUTPUT:
%   week = GPS week
%   sec_of_week = GPS Seconds of Week
%
% DESCRIPTION:
%   Conversion of Julian Day number to GPS week and
%	Seconds of Week reckoned from Saturday midnight

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.1 beta
%
% Copyright (C) 2009-2012 Mirko Reguzzoni, Eugenio Realini
%----------------------------------------------------------------------------------------------

deltat = julday - 2444244.5;
week = floor(deltat/7);
sec_of_week = (deltat - week*7)*86400;

