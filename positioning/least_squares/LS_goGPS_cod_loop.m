function LS_goGPS_cod_loop(time, Eph_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, snr_R, snr_M, iono, phase)

% SYNTAX:
%   LS_goGPS_cod_loop(time, Eph_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, snr_R, snr_M, iono, phase);
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
%   iono = ionosphere parameters
%   phase = L1 carrier (phase=1), L2 carrier (phase=2)
%
% DESCRIPTION:
%   Computation of the rover ground position (X,Y,Z).
%   Differential code positioning by least squares adjustment.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.2.0 beta
%
% Copyright (C) 2009-2011 Mirko Reguzzoni, Eugenio Realini
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
global azR elR distR azM elM distM
global PDOP HDOP VDOP

cov_pos_LS = [];

%--------------------------------------------------------------------------------------------
% SATELLITE SELECTION
%--------------------------------------------------------------------------------------------

if (length(phase) == 2)
    sat = find( (pr1_R ~= 0) & (pr1_M ~= 0) & ...
                (pr2_R ~= 0) & (pr2_M ~= 0) );
else
    if (phase == 1)
        sat = find( (pr1_R ~= 0) & (pr1_M ~= 0) );
    else
        sat = find( (pr2_R ~= 0) & (pr2_M ~= 0) );
    end
end

if (size(sat,1) >= 4)

   %ROVER positioning by means of Bancroft algorithm
   if (phase == 1)
      [pos_R, pos_S] = input_bancroft(pr1_R(sat), sat, time, Eph_R);
   else
      [pos_R, pos_S] = input_bancroft(pr2_R(sat), sat, time, Eph_R);
   end

   %-----------------------------------------------------------------------------------
   % CHECK SATELLITE ELEVATION, PIVOT AND CUT-OFF
   %-----------------------------------------------------------------------------------
   
   %initialization
   azR = zeros(32,1);
   elR = zeros(32,1);
   distR = zeros(32,1);
   azM = zeros(32,1);
   elM = zeros(32,1);
   distM = zeros(32,1);
   
   %satellite azimuth, elevation, ROVER-SATELLITE distance
   [azR(sat), elR(sat), distR(sat)] = topocent(pos_R, pos_S);
   
   %elevation cut-off
   sat_cutoff = find(elR > cutoff);
   sat = intersect(sat,sat_cutoff);
   
   %previous pivot
   pivot_old = 0;
   
   %actual pivot
   [null_max_elR, i] = max(elR(sat)); %#ok<ASGLU>
   pivot = sat(i);

   %--------------------------------------------------------------------------------------------
   % SATELLITE CONFIGURATION
   %--------------------------------------------------------------------------------------------
   
   %satellite configuration
   conf_sat = zeros(32,1);
   conf_sat(sat,1) = +1;
   
   %no cycle-slips when working with code only
   conf_cs = zeros(32,1);

   if (size(sat,1) >= 4)
       if (phase == 1)
           [pos_LS, cov_pos_LS, PDOP, HDOP, VDOP] = code_double_diff(pos_R(1:3), pr1_R(sat), snr_R(sat), pos_M, pr1_M(sat), snr_M(sat), time, sat, pivot, Eph_R, iono);
       else
           [pos_LS, cov_pos_LS, PDOP, HDOP, VDOP] = code_double_diff(pos_R(1:3), pr2_R(sat), snr_R(sat), pos_M, pr2_M(sat), snr_M(sat), time, sat, pivot, Eph_R, iono);
       end
   else
       if (~isempty(Xhat_t_t))
           pos_LS = Xhat_t_t([1,o1+1,o2+1]);
           pivot = 0;
       else
           pos_LS = NaN(3,1);
       end
   end
else
   if (~isempty(Xhat_t_t))
        pos_LS = Xhat_t_t([1,o1+1,o2+1]);
        pivot = 0;
    else
        pos_LS = NaN(3,1);
    end
end

if isempty(cov_pos_LS) %if it was not possible to compute the covariance matrix
    cov_pos_LS = sigmaq0 * eye(3);
end
sigmaq_pos_LS = diag(cov_pos_LS);

%-------------------------------------------------------------------------------

Xhat_t_t = zeros(o3,1);
Xhat_t_t(1) = pos_LS(1);
Xhat_t_t(o1+1) = pos_LS(2);
Xhat_t_t(o2+1) = pos_LS(3);
Cee(:,:) = zeros(o3);
Cee(1,1) = sigmaq_pos_LS(1);
Cee(o1+1,o1+1) = sigmaq_pos_LS(2);
Cee(o2+1,o2+1) = sigmaq_pos_LS(3);