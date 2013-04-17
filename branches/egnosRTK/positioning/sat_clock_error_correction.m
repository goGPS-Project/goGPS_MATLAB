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
%                           goGPS v0.3.1 beta
%
% Copyright (C) 2009-2012 Mirko Reguzzoni, Eugenio Realini
%
%----------------------------------------------------------------------------------------------

%if GLONASS
if (strcmp(char(Eph(31)),'R'))

    TauN   = Eph(2);
    GammaN = Eph(3);
    toe    = Eph(18);
    
    dt = check_t(time - toe);
    corr = -TauN + GammaN*dt;
    
else %if GPS/Galileo/QZSS/BeiDou

    af2 = Eph(2);
    af0 = Eph(19);
    af1 = Eph(20);
    toc = Eph(21);
    
    %consider BeiDou time (BDT) for BeiDou satellites
    if (strcmp(char(Eph(31)),'C'))
        time = time - 14;
    end
    
    dt = check_t(time - toc);
    corr = (af2 * dt + af1) * dt + af0;
end
