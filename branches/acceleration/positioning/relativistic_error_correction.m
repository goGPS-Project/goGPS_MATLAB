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
%                           goGPS v0.4.2 beta
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
%----------------------------------------------------------------------------------------------

if (Eph(4)~=0) %if not using SP3 ephemeris
    roota = Eph(4);
    ecc   = Eph(6);
    
    switch char(Eph(31))
        case 'G'
            %GM = goGNSS.GM_GPS;
            GM = 3.986005e14;
        case 'R'
            %GM = goGNSS.GM_GLO;
            GM = 3.9860044e14;
        case 'E'
            %GM = goGNSS.GM_GAL;
            GM = 3.986004418e14;
        case 'C'
            %GM = goGNSS.GM_BDS;
            GM = 3.986004418e14;
        case 'J'
            %GM = goGNSS.GM_QZS;
            GM = 3.986005e14;
        otherwise
            fprintf('Something went wrong in ecc_anomaly.m\nUnrecongized Satellite system!\n');
            %GM = goGNSS.GM_GPS;
            GM = 3.986005e14;
    end

    Ek = ecc_anomaly(time, Eph, GM);
    F = -2*sqrt(GM)/(goGNSS.V_LIGHT^2);
    corr = F * ecc * roota * sin(Ek);
else
    corr = -2*dot(XS,VS)/(goGNSS.V_LIGHT^2);
end
