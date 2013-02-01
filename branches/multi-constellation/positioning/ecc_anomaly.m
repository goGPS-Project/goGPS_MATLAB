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
%                           goGPS v0.3.1 beta
%
% Copyright (C) 2009-2012 Mirko Reguzzoni, Eugenio Realini
%
% Partially based on SATPOS.M (EASY suite) by Kai Borre
%----------------------------------------------------------------------------------------------

global GM

%get ephemerides
M0       =   Eph(3);
roota    =   Eph(4);
deltan   =   Eph(5);
ecc      =   Eph(6);
toe      =   Eph(18);

A = roota*roota;            %semi-major axis
tk = check_t(time-toe);     %time from the ephemerides reference epoch
n0 = sqrt(GM/A^3);          %computed mean motion [rad/sec]
n = n0+deltan;              %corrected mean motion [rad/sec]
Mk = M0+n*tk;               %mean anomaly
Mk = rem(Mk+2*pi,2*pi);
Ek = Mk;

max_iter = 12; %it was 10 when using only GPS (convergence was achieved at 4-6 iterations);
               % now it set to 12 because QZSS PRN 193 can take 11 iterations to converge

for i = 1 : max_iter
   Ek_old = Ek;
   Ek = Mk+ecc*sin(Ek);
   dEk = rem(Ek-Ek_old,2*pi);
   if abs(dEk) < 1.e-12
      break
   end
end

if (i == max_iter)
    fprintf('WARNING: Eccentric anomaly does not converge.\n')
end

Ek = rem(Ek+2*pi,2*pi);
