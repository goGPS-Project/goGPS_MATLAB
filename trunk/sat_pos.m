function [satp] = sat_pos(t, Eph)

% SYNTAX:
%   [satp] = sat_pos(t, Eph);
%
% INPUT:
%   t = GPS time
%   Eph = ephemerides matrix
%
% OUTPUT:
%   satp = satellite position (X,Y,Z)
%
% DESCRIPTION:
%   Computation of the satellite position (X,Y,Z) by means
%   of its ephemerides .

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.2 alpha
%
% Copyright (C) 2009-2010 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Graduate School for Creative Cities, Osaka City University, Japan
%
% Partially based on SATPOS.M (EASY suite) by Kai Borre
%----------------------------------------------------------------------------------------------


%-------------------------------------------------------------------------------
% VARIABLE INITIALIZATION
%-------------------------------------------------------------------------------

global Omegae_dot

%get ephemerides
roota    =   Eph(4);
ecc      =   Eph(6);
omega    =   Eph(7);
cuc      =   Eph(8);
cus      =   Eph(9);
crc      =  Eph(10);
crs      =  Eph(11);
i0       =  Eph(12);
idot     =  Eph(13);
cic      =  Eph(14);
cis      =  Eph(15);
Omega0   =  Eph(16);
Omegadot =  Eph(17);
toe      =  Eph(18);

%-------------------------------------------------------------------------------
% ALGORITHM FOR THE COMPUTATION OF THE SATELLITE COORDINATES
%-------------------------------------------------------------------------------

Ek = ecc_anomaly(t, Eph);

A = roota*roota;            %semi-major axis
tk = check_t(t-toe);        %time from the ephemerides reference epoch
fk = atan2(sqrt(1-ecc^2)*sin(Ek), cos(Ek)-ecc);
phi = fk+omega;
phi = rem(phi,2*pi);
u = phi               + cuc*cos(2*phi)+cus*sin(2*phi);
r = A*(1-ecc*cos(Ek)) + crc*cos(2*phi)+crs*sin(2*phi);
ik = i0+idot*tk       + cic*cos(2*phi)+cis*sin(2*phi);
Omega = Omega0+(Omegadot-Omegae_dot)*tk-Omegae_dot*toe;
Omega = rem(Omega+2*pi,2*pi);
x1 = cos(u)*r;
y1 = sin(u)*r;
%satellite coordinates (X,Y,Z)
satp(1,1) = x1*cos(Omega)-y1*cos(ik)*sin(Omega);
satp(2,1) = x1*sin(Omega)+y1*cos(ik)*cos(Omega);
satp(3,1) = y1*sin(ik);
