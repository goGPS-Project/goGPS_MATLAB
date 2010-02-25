function MQ_goGPS_loop(time, Eph_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, snr_R, snr_M, phase)

% SYNTAX:
%   MQ_goGPS_loop(time, Eph_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, snr_R, snr_M, phase);
%
% INPUT:
%   time = GPS time
%   Eph_R = satellite ephemerides
%   pos_M = MASTER position
%   pr1_R = ROVER code observations (L1 carrier)
%   pr1_M = MASTER code observations (L1 carrier)
%   pr2_R = ROVER code observations (L2 carrier)
%   pr2_M = MASTER code observations (L2 carrier)
%   snr_R = ROVER-SATELLITE signal-to-noise ratio
%   snr_M = MASTER-SATELLITE signal-to-noise ratio
%   phase = L1 carrier (phase=1), L2 carrier (phase=2)
%
% DESCRIPTION:
%   Computation of the rover ground position (X,Y,Z).
%   Differential code positioning by least squares adjustment.

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

global cutoff
global Xhat_t_t Cee conf_sat pivot
global o1 o2 o3
global azR elR distR
global sigmaq0

cov_pos_MQ = [];

%----------------%
%--- BANCROFT ---%
%----------------%

%visible satellites (ROVER)
if (phase == 1)
   sat_R = find(pr1_R ~= 0);
else
   sat_R = find(pr2_R ~= 0);
end

if (size(sat_R,1) >= 4)

   %ROVER positioning by means of Bancroft algorithm
   if (phase == 1)
      [pos_BAN, pos_SAT] = input_bancroft(pr1_R(sat_R), sat_R, time, Eph_R);
   else
      [pos_BAN, pos_SAT] = input_bancroft(pr2_R(sat_R), sat_R, time, Eph_R);
   end

   %-----------------------%
   %---  LEAST SQUARES  ---%
   %-----------------------%

   %visible satellites (MASTER)
   if (phase == 1)
      sat_M = find(pr1_M ~= 0);
   else
      sat_M = find(pr2_M ~= 0);
   end

   %-----------------------------------------------------------------------------------
   % SATELLITE ELEVATION AND PIVOT
   %-----------------------------------------------------------------------------------

   %initialization
   azR = zeros(32,1);
   elR = zeros(32,1);
   distR = zeros(32,1);

   %satellites in common and pivot selection
   sat = [];
   pivot = 0;
   max_elR = 0;

   for k = 1 : length(sat_R)

      %azimuth, elevation and ROVER-SATELLITE distance computation
      [azR(sat_R(k)), elR(sat_R(k)), distR(sat_R(k))] = topocent(pos_BAN(1:3), pos_SAT(k,1:3));

      %satellites in common
      if (ismember(sat_R(k),sat_M))

         sat = [sat; sat_R(k)];

         %pivot satellite
         if (elR(sat_R(k)) > max_elR)
            pivot = sat_R(k);
            max_elR = elR(sat_R(k));
         end
      end
   end

   %just for skyplot display
   conf_sat = zeros(32,1);
   conf_sat(sat) = +1;
   conf_sat(elR<cutoff) = 0;

   if (size(sat,1) >= 4)

      if (phase == 1)
         [pos_MQ, cov_pos_MQ] = code_double_diff(pos_BAN(1:3), pr1_R(sat), snr_R(sat), pos_M, pr1_M(sat), snr_M(sat), time, sat, pivot, Eph_R);
      else
         [pos_MQ, cov_pos_MQ] = code_double_diff(pos_BAN(1:3), pr2_R(sat), snr_R(sat), pos_M, pr2_M(sat), snr_M(sat), time, sat, pivot, Eph_R);
      end
   else
      fprintf('Less than 4 satellites in common\n');
      pos_MQ = Xhat_t_t([1,o1+1,o2+1]);
   end
else
   fprintf('Less than 4 satellites\n');
   pos_MQ = Xhat_t_t([1,o1+1,o2+1]);
end

if isempty(cov_pos_MQ) %if it was not possible to compute the covariance matrix
    cov_pos_MQ = sigmaq0 * eye(3);
end
sigmaq_pos_MQ = diag(cov_pos_MQ);

%-------------------------------------------------------------------------------

Xhat_t_t = zeros(o3,1);
Xhat_t_t(1) = pos_MQ(1);
Xhat_t_t(o1+1) = pos_MQ(2);
Xhat_t_t(o2+1) = pos_MQ(3);
Cee(:,:) = zeros(o3);
Cee(1,1) = sigmaq_pos_MQ(1);
Cee(o1+1,o1+1) = sigmaq_pos_MQ(2);
Cee(o2+1,o2+1) = sigmaq_pos_MQ(3);