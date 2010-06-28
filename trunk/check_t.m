function [corrTime] = check_t(time)

% SYNTAX:
%   [corrTime] = check_t(time);
%
% INPUT:
%   t = GPS time
%
% OUTPUT:
%   tt = corrected GPS time
%
% DESCRIPTION:
%   Function accounting for beginning or end of week crossover.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.1 alpha
%
% Copyright (C) Kai Borre
% Kai Borre 04-01-96
%
% Adapted by Mirko Reguzzoni, Eugenio Realini, 2009
%
%----------------------------------------------------------------------------------------------

half_week = 302400;     % seconds

corrTime = time;

if time > half_week
    corrTime = time - 2*half_week;
elseif time < -half_week
    corrTime = time + 2*half_week;
end
