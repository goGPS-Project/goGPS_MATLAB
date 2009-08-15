function MQ_goGPS_SA_loop(time, Eph_R, pr1_R, pr2_R, phase)

% SYNTAX:
%   MQ_goGPS_loop(time, Eph_R, pr1_R, pr1_M, pr2_R, pr2_M, phase);
%
% INPUT:
%   time = GPS time
%   Eph_R = satellite ephemerides
%   pr1_R = ROVER code observations (L1 carrier)
%   pr2_R = ROVER code observations (L2 carrier)
%   phase = L1 carrier (phase=1), L2 carrier (phase=2)
%
% DESCRIPTION:
%   Computation of the rover ground position (X,Y,Z).
%   Standalone code positioning by least squares adjustment.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 pre-alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini*
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
%
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

global a f
global cutoff
global Xhat_t_t Cee conf_sat pivot
global o1 o2 o3
global azR elR distR

%parametro di cutoff per minimi quadrati
cutoff = 15;

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
   
   %-------------------------------------%
   %---  LEAST SQUARES (STAND-ALONE)  ---%
   %-------------------------------------%

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
       [azR(sat_R(k)), elR(sat_R(k)), distR(sat_R(k))] = topocent(pos_BAN(1:3), pos_SAT(k,1:3), a, f);

       %satellites in common
       sat = [sat; sat_R(k)];

       %pivot satellite
       if (elR(sat_R(k)) > max_elR)
           pivot = sat_R(k);
           max_elR = elR(sat_R(k));
       end

   end

   %just for skyplot display
   conf_sat = zeros(32,1);
   conf_sat(sat) = +1;
   conf_sat(elR<cutoff) = 0;
   
  if (size(sat,1) >= 4)

      if (phase == 1)
         [pos_SA, cov_pos_SA] = code_SA(pos_BAN(1:3), pr1_R, time, Eph_R);
      else
         [pos_SA, cov_pos_SA] = code_SA(pos_BAN(1:3), pr2_R, time, Eph_R);
      end

   else
      fprintf('Less than 4 common satellites\n');
	  pos_SA = Xhat_t_t([1,o1+1,o2+1]);
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