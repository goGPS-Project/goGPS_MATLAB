function icol = find_eph(Eph, sat, time)

% SYNTAX:
%   icol = find_eph(Eph, sat, time);
%
% INPUT:
%   Eph = ephemerides matrix
%   sat = satellite index
%   time = GPS time
%
% OUTPUT:
%   icol = column index for the selected satellite
%
% DESCRIPTION:
%   Selection of the column corresponding to the specified satellite
%   (with respect to the specified GPS time) in the ephemerides matrix.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) Kai Borre
% Kai Borre and C.C. Goad 11-26-96
%
% Adapted by Mirko Reguzzoni, Eugenio Realini, 2009
%----------------------------------------------------------------------------------------------

isat = find(Eph(30,:) == sat);

n = size(isat,2);
if (n == 0)
    icol = [];
    return
end
icol = isat(1);

delta = 0;

%consider BeiDou time (BDT) for BeiDou satellites
if (strcmp(char(Eph(31)),'C'))
    delta = 14;
    time = time - delta;
end

time_eph = Eph(32,icol);
dtmin = time_eph - time;
for t = isat
    
    time_eph = Eph(32,t);
    dt = time_eph - time;
    if (abs(dt) < abs(dtmin))
        icol = t;
        dtmin = dt;
    end
end

%maximum interval from ephemeris reference time
fit_interval = Eph(29,icol);
if (fit_interval ~= 0)
    dtmax = fit_interval*3600/2;
else
    switch (char(Eph(31,icol)))
        case 'R' %GLONASS
            dtmax = 900;
        case 'J' %QZSS
            dtmax = 3600;
        otherwise
            dtmax = 7200;
    end
end

if (fix(abs(dtmin)) - delta > dtmax)
    icol = [];
    return
end

%check satellite health
%the second and third conditions are temporary (QZSS and Galileo health flag is kept on for tests)
if (Eph(27,icol) ~= 0 && ~strcmp(char(Eph(31,icol)),'J') && ~strcmp(char(Eph(31,icol)),'E'))
    icol = [];
    return
end
