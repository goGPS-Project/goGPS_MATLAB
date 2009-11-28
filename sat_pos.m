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
%                           goGPS v0.1 alpha
%
% Copyright (C) Kai Borre 
% Kai Borre 04-09-96
%
% Adapted by Mirko Reguzzoni, Eugenio Realini, 2009
%
%----------------------------------------------------------------------------------------------


%-------------------------------------------------------------------------------
% VARIABLE INITIALIZATION
%-------------------------------------------------------------------------------

global GM Omegae_dot

%get ephemerides
svprn    =   Eph(1);		% not used
af2      =   Eph(2);		% not used
M0       =   Eph(3);
roota    =   Eph(4);
deltan   =   Eph(5);
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
af0      =  Eph(19);		% not used
af1      =  Eph(20);		% not used
tom      =  Eph(21);		% not used

%-------------------------------------------------------------------------------
% ALGORITHM FOR THE COMPUTATION OF THE SATELLITE COORDINATES
%-------------------------------------------------------------------------------

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
      break;
   end
end

if (i == 10)
    warning('Eccentric anomaly does not converge!!')
end

Ek = rem(Ek+2*pi,2*pi);
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
