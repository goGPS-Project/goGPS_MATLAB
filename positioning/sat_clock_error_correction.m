function [corr] = sat_clock_error_correction(time, Eph)

% SYNTAX:
%   [corr] = sat_clock_error_correction(time, Eph);
%
% INPUT:
%   time = GPS time
%   Eph = satellite ephemeris vector
%
% OUTPUT:
%   corr = satellite clock correction term
%
% DESCRIPTION:
%   Computation of the satellite clock error correction term. From the
%   Interface Specification document revision E (IS-GPS-200E), page 86.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.0 beta
%
% Copyright (C) 2009-2012 Mirko Reguzzoni, Eugenio Realini
%
%----------------------------------------------------------------------------------------------

af2   = Eph(2);
af0   = Eph(19);
af1   = Eph(20);
tom   = Eph(21);

dt = check_t(time - tom);
corr = (af2 * dt + af1) * dt + af0;
