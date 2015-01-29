function [satp, satv] = satellite_orbits(t, Eph, sat, sbas)

% SYNTAX:
%   [satp, satv] = satellite_orbits(t, Eph, sat, sbas);
%
% INPUT:
%   t = clock-corrected GPS time
%   Eph  = ephemeris matrix
%   sat  = satellite index
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
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%----------------------------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%----------------------------------------------------------------------------------------------

switch char(Eph(31))
    case 'G'
        Omegae_dot = goGNSS.OMEGAE_DOT_GPS;
    case 'R'
        Omegae_dot = goGNSS.OMEGAE_DOT_GLO;
    case 'E'
        Omegae_dot = goGNSS.OMEGAE_DOT_GAL;
    case 'C'
        Omegae_dot = goGNSS.OMEGAE_DOT_BDS;
    case 'J'
        Omegae_dot = goGNSS.OMEGAE_DOT_QZS;
    otherwise
        fprintf('Something went wrong in satellite_orbits.m\nUnrecongized Satellite system!\n');
        Omegae_dot = goGNSS.OMEGAE_DOT_GPS;
end

%consider BeiDou time (BDT) for BeiDou satellites
if (strcmp(char(Eph(31)),'C'))
    t = t - 14;
end

%GPS/Galileo/BeiDou/QZSS satellite coordinates computation
if (~strcmp(char(Eph(31)),'R'))
    
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
    time_eph  = Eph(32);
    
    %SBAS satellite coordinate corrections
    if (~isempty(sbas))
        dx_sbas = sbas.dx(sat);
        dy_sbas = sbas.dy(sat);
        dz_sbas = sbas.dz(sat);
    else
        dx_sbas = 0;
        dy_sbas = 0;
        dz_sbas = 0;
    end
    
    %-------------------------------------------------------------------------------
    % ALGORITHM FOR THE COMPUTATION OF THE SATELLITE COORDINATES (IS-GPS-200E)
    %-------------------------------------------------------------------------------
    
    %eccentric anomaly
    [Ek, n] = ecc_anomaly(t, Eph);
    
    %cr = goGNSS.CIRCLE_RAD;
    cr = 6.283185307179600;
    
    A = roota*roota;             %semi-major axis
    tk = check_t(t - time_eph);  %time from the ephemeris reference epoch
    
    fk = atan2(sqrt(1-ecc^2)*sin(Ek), cos(Ek) - ecc);    %true anomaly
    phik = fk + omega;                           %argument of latitude
    phik = rem(phik,cr);
    
    uk = phik                + cuc*cos(2*phik) + cus*sin(2*phik); %corrected argument of latitude
    rk = A*(1 - ecc*cos(Ek)) + crc*cos(2*phik) + crs*sin(2*phik); %corrected radial distance
    ik = i0 + IDOT*tk        + cic*cos(2*phik) + cis*sin(2*phik); %corrected inclination of the orbital plane
    
    %satellite positions in the orbital plane
    x1k = cos(uk)*rk;
    y1k = sin(uk)*rk;
    
    %if GPS/Galileo/QZSS or MEO/IGSO BeiDou satellite
    if (~strcmp(char(Eph(31)),'C') || (strcmp(char(Eph(31)),'C') && Eph(1) > 5))
        
        %corrected longitude of the ascending node
        Omegak = Omega0 + (Omega_dot - Omegae_dot)*tk - Omegae_dot*toe;
        Omegak = rem(Omegak + cr, cr);
        
        %satellite Earth-fixed coordinates (X,Y,Z)
        xk = x1k*cos(Omegak) - y1k*cos(ik)*sin(Omegak);
        yk = x1k*sin(Omegak) + y1k*cos(ik)*cos(Omegak);
        zk = y1k*sin(ik);
        
        %apply SBAS corrections (if available)
        satp(1,1) = xk + dx_sbas;
        satp(2,1) = yk + dy_sbas;
        satp(3,1) = zk + dz_sbas;
        
    else %if GEO BeiDou satellite (ranging code number <= 5)
        
        %corrected longitude of the ascending node
        Omegak = Omega0 + Omega_dot*tk - Omegae_dot*toe;
        Omegak = rem(Omegak + cr, cr);
        
        %satellite coordinates (X,Y,Z) in inertial system
        xgk = x1k*cos(Omegak) - y1k*cos(ik)*sin(Omegak);
        ygk = x1k*sin(Omegak) + y1k*cos(ik)*cos(Omegak);
        zgk = y1k*sin(ik);
        
        %store inertial coordinates in a vector
        Xgk = [xgk; ygk; zgk];
        
        %rotation matrices from inertial system to CGCS2000
        Rx = [1        0          0;
              0 +cosd(-5) +sind(-5);
              0 -sind(-5) +cosd(-5)];
        
        oedt = Omegae_dot*tk;
        
        Rz = [+cos(oedt) +sin(oedt) 0;
              -sin(oedt) +cos(oedt) 0;
              0           0         1];
        
        %apply the rotations
        Xk = Rz*Rx*Xgk;
        
        xk = Xk(1);
        yk = Xk(2);
        zk = Xk(3);
        
        %store CGCS2000 coordinates
        satp(1,1) = xk;
        satp(2,1) = yk;
        satp(3,1) = zk;
    end
    
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
    
else %GLONASS satellite coordinates computation (GLONASS-ICD 5.1)

    time_eph = Eph(32); %ephemeris reference time

    X   = Eph(5);  %satellite X coordinate at ephemeris reference time
    Y   = Eph(6);  %satellite Y coordinate at ephemeris reference time
    Z   = Eph(7);  %satellite Z coordinate at ephemeris reference time

    Xv  = Eph(8);  %satellite velocity along X at ephemeris reference time
    Yv  = Eph(9);  %satellite velocity along Y at ephemeris reference time
    Zv  = Eph(10); %satellite velocity along Z at ephemeris reference time

    Xa  = Eph(11); %acceleration due to lunar-solar gravitational perturbation along X at ephemeris reference time
    Ya  = Eph(12); %acceleration due to lunar-solar gravitational perturbation along Y at ephemeris reference time
    Za  = Eph(13); %acceleration due to lunar-solar gravitational perturbation along Z at ephemeris reference time
    %NOTE:  Xa,Ya,Za are considered constant within the integration interval (i.e. toe ?}15 minutes)
    
    %integration step
    int_step = 60; %[s]
    
    %time from the ephemeris reference epoch
    tk = check_t(t - time_eph);
    
    %number of iterations on "full" steps
    n = floor(abs(tk/int_step));

    %array containing integration steps (same sign as tk)
    ii = ones(n,1)*int_step*(tk/abs(tk));
    
    %check residual iteration step (i.e. remaining fraction of int_step)
    int_step_res = rem(tk,int_step);

    %adjust the total number of iterations and the array of iteration steps
    if (int_step_res ~= 0)
        n = n + 1;
        ii = [ii; int_step_res];
    end
    
    %numerical integration steps (i.e. re-calculation of satellite positions from toe to tk)
    pos = [X Y Z];
    vel = [Xv Yv Zv];
    acc = [Xa Ya Za];

    for s = 1 : n

        %Runge-Kutta numerical integration algorithm
        %
        %step 1
        pos1 = pos;
        vel1 = vel;
        [pos1_dot, vel1_dot] = satellite_motion_diff_eq(pos1, vel1, acc, goGNSS.ELL_A_GLO, goGNSS.GM_GLO, goGNSS.J2_GLO, goGNSS.OMEGAE_DOT_GLO);
        %
        %step 2
        pos2 = pos + pos1_dot*ii(s)/2;
        vel2 = vel + vel1_dot*ii(s)/2;
        [pos2_dot, vel2_dot] = satellite_motion_diff_eq(pos2, vel2, acc, goGNSS.ELL_A_GLO, goGNSS.GM_GLO, goGNSS.J2_GLO, goGNSS.OMEGAE_DOT_GLO);
        %
        %step 3
        pos3 = pos + pos2_dot*ii(s)/2;
        vel3 = vel + vel2_dot*ii(s)/2;
        [pos3_dot, vel3_dot] = satellite_motion_diff_eq(pos3, vel3, acc, goGNSS.ELL_A_GLO, goGNSS.GM_GLO, goGNSS.J2_GLO, goGNSS.OMEGAE_DOT_GLO);
        %
        %step 4
        pos4 = pos + pos3_dot*ii(s);
        vel4 = vel + vel3_dot*ii(s);
        [pos4_dot, vel4_dot] = satellite_motion_diff_eq(pos4, vel4, acc, goGNSS.ELL_A_GLO, goGNSS.GM_GLO, goGNSS.J2_GLO, goGNSS.OMEGAE_DOT_GLO);
        %
        %final position and velocity
        pos = pos + (pos1_dot + 2*pos2_dot + 2*pos3_dot + pos4_dot)*ii(s)/6;
        vel = vel + (vel1_dot + 2*vel2_dot + 2*vel3_dot + vel4_dot)*ii(s)/6;
    end

    %transformation from PZ-90.02 to WGS-84 (G1150)
    satp(1,1) = pos(1) - 0.36;
    satp(2,1) = pos(2) + 0.08;
    satp(3,1) = pos(3) + 0.18;
    
    %satellite velocity
    satv(1,1) = vel(1);
    satv(2,1) = vel(2);
    satv(3,1) = vel(3);
end
