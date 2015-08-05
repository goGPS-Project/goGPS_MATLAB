function [corrTime] = check_t(time)

% SYNTAX:
%   [corrTime] = check_t(time);
%
% INPUT:
%   time = GPS time
%
% OUTPUT:
%   corrTime = corrected GPS time
%
% DESCRIPTION:
%   Function accounting for beginning or end of week crossover. From the
%   Interface Specification document revision E (IS-GPS-200E), page 93.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) Kai Borre
% Kai Borre 04-01-96
%
% Adapted by Mirko Reguzzoni, Eugenio Realini, 2009
%----------------------------------------------------------------------------------------------

half_week = 302400;     % seconds

corrTime = time;

if time > half_week
    corrTime = time - 2*half_week;
elseif time < -half_week
    corrTime = time + 2*half_week;
end
