function [Ek] = ecc_anomaly(t, Eph)

% SYNTAX:
%   [Ek] = ecc_anomaly(t, Eph);
%
% INPUT:
%   t = GPS time
%   Eph = ephemerides matrix
%
% OUTPUT:
%   Ek = eccentric anomaly
%
% DESCRIPTION:
%   Computation of the eccentric anomaly.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 beta
%
% Copyright (C) 2009-2010 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Graduate School for Creative Cities, Osaka City University, Japan
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
tk = check_t(t-toe);        %time from the ephemerides reference epoch
n0 = sqrt(GM/A^3);          %computed mean motion [rad/sec]
n = n0+deltan;              %corrected mean motion [rad/sec]
Mk = M0+n*tk;               %mean anomaly
Mk = rem(Mk+2*pi,2*pi);
Ek = Mk;

for i = 1:10
   Ek_old = Ek;
   Ek = Mk+ecc*sin(Ek);
   dEk = rem(Ek-Ek_old,2*pi);
   if abs(dEk) < 1.e-12
      break
   end
end

if (i == 10)
    fprintf('Eccentric anomaly does not converge!!\n')
end

Ek = rem(Ek+2*pi,2*pi);