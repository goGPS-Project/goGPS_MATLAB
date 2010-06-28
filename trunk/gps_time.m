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
%                           goGPS v0.1.1 alpha
%
% Copyright (C) Kai Borre
% Kai Borre 05-20-96
%
% Adapted by Mirko Reguzzoni, Eugenio Realini, 2009
%
%----------------------------------------------------------------------------------------------


a = floor(julday+.5);
b = a+1537;
c = floor((b-122.1)/365.25);
e = floor(365.25*c);
f = floor((b-e)/30.6001);
d = b-e-floor(30.6001*f)+rem(julday+.5,1);
day_of_week = rem(floor(julday+.5),7);
week = floor((julday-2444244.5)/7);
% We add +1 as the GPS week starts at Saturday midnight
sec_of_week = (rem(d,1)+day_of_week+1)*86400;
