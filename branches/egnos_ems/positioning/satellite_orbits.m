function [satp, satv] = satellite_orbits(t, Eph, sat, sbas)

% SYNTAX:
%   [satp, satv] = satellite_orbits(t, Eph, sat, sbas);
%
% INPUT:
%   t = clock-corrected GPS time
%   Eph  = ephemeris matrix
%   sat  = satellite PRN
%   sbas = SBAS corrections
%
% OUTPUT:
%   satp = satellite position (X,Y,Z)
%   satv = satellite velocity
%
% DESCRIPTION:
%   Computation of the satellite position (X,Y,Z) and velocity by means
%   of its ephemerides.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.1 beta
%
% Copyright (C) 2009-2012 Mirko Reguzzoni, Eugenio Realini
%
% Partially based on SATPOS.M (EASY suite) by Kai Borre
%----------------------------------------------------------------------------------------------


%-------------------------------------------------------------------------------
% VARIABLE INITIALIZATION
%-------------------------------------------------------------------------------

global Omegae_dot

%get ephemerides
roota     = Eph(4);
ecc       = Eph(6);
omega     = Eph(7);
cuc       = Eph(8);
cus       = Eph(9);
crc       = Eph(10);
crs       = Eph(11);
i0        = Eph(12);
IDOT      = Eph(13);
cic       = Eph(14);
cis       = Eph(15);
Omega0    = Eph(16);
Omega_dot = Eph(17);
toe       = Eph(18);

%-------------------------------------------------------------------------------
% ALGORITHM FOR THE COMPUTATION OF THE SATELLITE COORDINATES
%-------------------------------------------------------------------------------

[Ek, n] = ecc_anomaly(t, Eph);

A = roota*roota;            %semi-major axis
tk = check_t(t-toe);        %time from the ephemerides reference epoch
fk = atan2(sqrt(1-ecc^2)*sin(Ek), cos(Ek)-ecc);
phik = fk+omega;
phik = rem(phik,2*pi);
uk = phik               + cuc*cos(2*phik)+cus*sin(2*phik);
rk = A*(1-ecc*cos(Ek)) + crc*cos(2*phik)+crs*sin(2*phik);
ik = i0+IDOT*tk       + cic*cos(2*phik)+cis*sin(2*phik);
Omegak = Omega0+(Omega_dot-Omegae_dot)*tk-Omegae_dot*toe;
Omegak = rem(Omegak+2*pi,2*pi);
x1k = cos(uk)*rk;
y1k = sin(uk)*rk;
%satellite coordinates (X,Y,Z)
xk = x1k*cos(Omegak)-y1k*cos(ik)*sin(Omegak);
yk = x1k*sin(Omegak)+y1k*cos(ik)*cos(Omegak);
zk = y1k*sin(ik);

satp(1,1) = xk + sbas.dx(sat);
satp(2,1) = yk + sbas.dy(sat);
satp(3,1) = zk + sbas.dz(sat);

%-------------------------------------------------------------------------------
% ALGORITHM FOR THE COMPUTATION OF THE SATELLITE VELOCITY (as in Remondi,
% GPS Solutions (2004) 8:181-183 )
%-------------------------------------------------------------------------------
if (nargout > 1)
    Mk_dot = n;
    Ek_dot = Mk_dot/(1-ecc*cos(Ek));
    fk_dot = sin(Ek)*Ek_dot*(1+ecc*cos(fk)) / ((1-cos(Ek)*ecc)*sin(fk));
    phik_dot = fk_dot;
    uk_dot = phik_dot + 2*(cus*cos(2*phik)-cuc*sin(2*phik))*phik_dot;
    rk_dot = A*ecc*sin(Ek)*Ek_dot + 2*(crs*cos(2*phik)-crc*sin(2*phik))*phik_dot;
    ik_dot = IDOT + 2*(cis*cos(2*phik)-cic*sin(2*phik))*phik_dot;
    Omegak_dot = Omega_dot - Omegae_dot;
    x1k_dot = rk_dot*cos(uk) - y1k*uk_dot;
    y1k_dot = rk_dot*sin(uk) + x1k*uk_dot;
    xk_dot = x1k_dot*cos(Omegak) - y1k_dot*cos(ik)*sin(Omegak) + y1k*sin(ik)*sin(Omegak)*ik_dot - yk*Omegak_dot;
    yk_dot = x1k_dot*sin(Omegak) + y1k_dot*cos(ik)*cos(Omegak) - y1k*sin(ik)*ik_dot*cos(Omegak) + xk*Omegak_dot;
    zk_dot = y1k_dot*sin(ik) + y1k*cos(ik)*ik_dot;
    
    satv(1,1) = xk_dot;
    satv(2,1) = yk_dot;
    satv(3,1) = zk_dot;
end
