function [doppler_app1, doppler_app2] = doppler_shift_approx(pos_R, vel_R, pos_S, vel_S, time, rec_clock_drift, sat, Eph, lambda)

% SYNTAX:
%   [doppler_app1, doppler_app2] = doppler_shift_approx(pos_R, vel_R, pos_S, vel_S, time, rec_clock_drift, sat, Eph, lambda);
%
% INPUT:
%   pos_R = approximate rover position
%   vel_R = approximate rover velocity
%   pos_S = satellite position at transmission time
%   vel_S = satellite velocity at transmission time
%   time = signal transmission GPS time
%   sat = satellite id number
%   rec_clock_drift = receiver clock drift
%   Eph = satellite ephemerides matrix
%   lambda = vector containing GNSS wavelengths for the selected satellite
%
% OUTPUT:
%   doppler_app1 = approximate doppler shift on L1
%   doppler_app2 = approximate doppler shift on L2
%
% DESCRIPTION:
%   Computation of an approximate value of the doppler shift observation.

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

v_light = goGNSS.V_LIGHT;

LOS  = pos_S - pos_R;             %receiver-satellite line-of-sight vector
LOSu = LOS / norm(LOS);           %receiver-satellite line-of-sight unit vector [= LOS / sqrt(LOS(1)^2 + LOS(2)^2 + LOS(3)^2)]
vrel = vel_S - vel_R;             %receiver-satellite relative velocity vector
radial_vel = dot(vrel,LOSu);      %receiver-satellite radial velocity [= vrel(1)*LOSu(1) + vrel(2)*LOSu(2) + vrel(3)*LOSu(3)]
k = find_eph(Eph, sat, time);
if (~isempty(k))
    af2 = Eph(2,k);
    af1 = Eph(20,k);
    toc = Eph(21,k);
    if (strcmp(char(Eph(31)),'C')); time = time - 14; end %consider BeiDou time (BDT) for BeiDou satellites
    dt = check_t(time - toc);
    sat_clock_drift = af1 + 2*af2*dt; %satellite clock drift
    doppler_app1 = -(radial_vel + v_light*(rec_clock_drift - sat_clock_drift)) / lambda(1);
    doppler_app2 = -(radial_vel + v_light*(rec_clock_drift - sat_clock_drift)) / lambda(2);
else
    doppler_app1 = 0;
    doppler_app2 = 0;
end
