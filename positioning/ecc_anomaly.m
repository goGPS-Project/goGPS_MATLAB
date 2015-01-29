function [Ek, n] = ecc_anomaly(time, Eph)

% SYNTAX:
%   [Ek, n] = ecc_anomaly(time, Eph);
%
% INPUT:
%   time = GPS time
%   Eph = ephemerides matrix
%
% OUTPUT:
%   Ek = eccentric anomaly
%   n = corrected mean motion [rad/sec]
%
% DESCRIPTION:
%   Computation of the eccentric anomaly.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Partially based on SATPOS.M (EASY suite) by Kai Borre
%----------------------------------------------------------------------------------------------

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

%get ephemerides
M0       = Eph(3);
roota    = Eph(4);
deltan   = Eph(5);
ecc      = Eph(6);
time_eph = Eph(32);

%cr = goGNSS.CIRCLE_RAD;
cr = 6.283185307179600;

A  = roota*roota;              %semi-major axis
tk = check_t(time - time_eph); %time from the ephemerides reference epoch
n0 = sqrt(GM/A^3);             %computed mean motion [rad/sec]
n  = n0 + deltan;              %corrected mean motion [rad/sec]
Mk = M0 + n*tk;                %mean anomaly
Mk = rem(Mk+cr,cr);
Ek = Mk;

max_iter = 12; %it was 10 when using only GPS (convergence was achieved at 4-6 iterations);
               % now it set to 12 because QZSS PRN 193 can take 11 iterations to converge

for i = 1 : max_iter
   Ek_old = Ek;
   Ek = Mk+ecc*sin(Ek);
   dEk = rem(Ek-Ek_old,cr);
   if abs(dEk) < 1.e-12
      break
   end
end

if (i == max_iter)
    fprintf('WARNING: Eccentric anomaly does not converge.\n')
end

Ek = rem(Ek+cr,cr);
