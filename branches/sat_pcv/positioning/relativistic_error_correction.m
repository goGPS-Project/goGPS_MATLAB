function [corr] = relativistic_error_correction(time, Eph, XS, VS)

% SYNTAX:
%   [corr] = relativistic_error_correction(time, Eph, XS, VS);
%
% INPUT:
%   time = GPS time
%   Eph = satellite ephemeris vector
%   XS = satellite position (X,Y,Z)
%   VS = satellite velocity (X,Y,Z)
%
% OUTPUT:
%   corr = relativistic correction term
%
% DESCRIPTION:
%   Computation of the relativistic correction term. From the
%   Interface Specification document revision E (IS-GPS-200E), page 86.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
%----------------------------------------------------------------------------------------------

if (Eph(4)~=0) %if not using SP3 ephemeris
    roota = Eph(4);
    ecc   = Eph(6);
    
    Ek = ecc_anomaly(time, Eph);
    corr = -4.442807633e-10 * ecc * roota * sin(Ek);
else
    corr = -2*dot(XS,VS)/(goGNSS.V_LIGHT^2);
end
