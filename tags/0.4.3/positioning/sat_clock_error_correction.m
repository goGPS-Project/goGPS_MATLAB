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
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
%----------------------------------------------------------------------------------------------

%if GLONASS
if (strcmp(char(Eph(31)),'R'))

    TauN     = Eph(2);
    GammaN   = Eph(3);
    ref_toe  = Eph(32);
    
    dt = check_t(time - ref_toe);
    corr = -TauN + GammaN*dt;
    
else %if GPS/Galileo/QZSS/BeiDou

    af2 = Eph(2);
    af0 = Eph(19);
    af1 = Eph(20);
    ref_toc = Eph(33);
    
    %consider BeiDou time (BDT) for BeiDou satellites
    if (strcmp(char(Eph(31)),'C'))
        time = time - 14;
    end
    
    dt = check_t(time - ref_toc);
    corr = (af2 * dt + af1) * dt + af0;
end
