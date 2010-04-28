function MQ_goGPS_SA_loop(time, Eph_R, pr1_R, pr2_R, snr_R, phase)

% SYNTAX:
%   MQ_goGPS_SA_loop(time, Eph_R, pr1_R, pr1_M, pr2_R, snr_R, phase);
%
% INPUT:
%   time = GPS time
%   Eph_R = satellite ephemerides
%   pr1_R = ROVER code observations (L1 carrier)
%   pr2_R = ROVER code observations (L2 carrier)
%   snr_R = ROVER-SATELLITE signal-to-noise ratio
%   phase = L1 carrier (phase=1), L2 carrier (phase=2)
%
% DESCRIPTION:
%   Computation of the rover ground position (X,Y,Z).
%   Standalone code positioning by least squares adjustment.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 beta
%
% Copyright (C) 2009-2010 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Graduate School for Creative Cities, Osaka City University, Japan
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

global sigmaq0
global cutoff o1 o2 o3

global Xhat_t_t Cee conf_sat conf_cs pivot pivot_old
global azR elR distR
global PDOP HDOP VDOP

%covariance matrix initialization
cov_pos_SA = [];

%----------------%
%--- BANCROFT ---%
%----------------%

%visible satellites (ROVER)
if (phase == 1)
   sat_pr = find(pr1_R ~= 0);
else
   sat_pr = find(pr2_R ~= 0);
end

if (size(sat_pr,1) >= 4)

   %ROVER positioning by means of Bancroft algorithm
   if (phase == 1)
      [pos_R, pos_SAT] = input_bancroft(pr1_R(sat_pr), sat_pr, time, Eph_R);
   else
      [pos_R, pos_SAT] = input_bancroft(pr2_R(sat_pr), sat_pr, time, Eph_R);
   end
   
   pos_R = pos_R(1:3);
   pos_SAT = pos_SAT(:,1:3);
   
   %-----------------------------------------------------------------------------------
   % CHECK SATELLITE ELEVATION, PIVOT AND CUT-OFF
   %-----------------------------------------------------------------------------------
   
   %initialization
   azR = zeros(32,1);
   elR = zeros(32,1);
   distR = zeros(32,1);
   
   %satellite azimuth, elevation, ROVER-SATELLITE distance
   [azR(sat_pr), elR(sat_pr), distR(sat_pr)] = topocent(pos_R, pos_SAT);
   
   %elevation cut-off
   sat_cutoff = find(elR > cutoff);
   sat_pr = intersect(sat_pr,sat_cutoff);
   
   %previous pivot
   pivot_old = 0;
   
   %actual pivot
   [null_max_elR, i] = max(elR(sat_pr)); %#ok<ASGLU>
   pivot = sat_pr(i);

   %--------------------------------------------------------------------------------------------
   % SATELLITE CONFIGURATION
   %--------------------------------------------------------------------------------------------
   
   %satellite configuration
   conf_sat = zeros(32,1);
   conf_sat(sat_pr,1) = +1;
   
   %no cycle-slips when working with code only
   conf_cs = zeros(32,1);
   
   if (phase == 1)
       [pos_SA, cov_pos_SA, PDOP, HDOP, VDOP] = code_SA(pos_R(1:3), pr1_R(sat_pr), snr_R(sat_pr), sat_pr, time, Eph_R);
   else
       [pos_SA, cov_pos_SA, PDOP, HDOP, VDOP] = code_SA(pos_R(1:3), pr2_R(sat_pr), snr_R(sat_pr), sat_pr, time, Eph_R);
   end
else
   fprintf('Less than 4 satellites\n');
   pos_SA = Xhat_t_t([1,o1+1,o2+1]);
end

if isempty(cov_pos_SA) %if it was not possible to compute the covariance matrix
    cov_pos_SA = sigmaq0 * eye(3);
end
sigmaq_pos_SA = diag(cov_pos_SA);

%-------------------------------------------------------------------------------

Xhat_t_t = zeros(o3,1);
Xhat_t_t(1) = pos_SA(1);
Xhat_t_t(o1+1) = pos_SA(2);
Xhat_t_t(o2+1) = pos_SA(3);
Cee(:,:) = zeros(o3);
Cee(1,1) = sigmaq_pos_SA(1);
Cee(o1+1,o1+1) = sigmaq_pos_SA(2);
Cee(o2+1,o2+1) = sigmaq_pos_SA(3);