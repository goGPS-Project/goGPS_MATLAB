function [matrix_oss, snr] = RINEX_get_obs(file_RINEX, num_sat, num_type_oss, col_L1)

% SYNTAX:
%   [matrix_oss, snr] = RINEX_get_obs(file_RINEX, num_sat, num_type_oss, col_L1);
%
% INPUT:
%   file_RINEX = observation RINEX file
%   num_sat = number of visible satellites
%   num_type_oss = number of types of observations
%   col_L1 = column index for the phase observation L1
%
% OUTPUT:
%   matrix_oss = observation matrix
%   snr = signal-to-noise ratio
%
% DESCRIPTION:
%   Acquisition of RINEX observation data (code, phase).
%   Acquisition of signal-to-noise ratio.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 alpha
%
% Copyright (C) 2009-2010 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Graduate School for Creative Cities, Osaka City University, Japan
%
% Partially based on GRABDATA.M (EASY suite) by Kai Borre
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

%variable initialization
matrix_oss = zeros(num_sat, num_type_oss);
snr = zeros(num_sat,1);

%less than 5 observation types --> just 1 line to be read
if num_type_oss <= 5
   for u = 1:num_sat
      lin = fgetl(file_RINEX);
      for k = 1:num_type_oss
          %data save
          if (length(lin) < 2+16*(k-1)) | (isempty(sscanf(lin(2+16*(k-1):16*k-2),'%f')))
              matrix_oss(u,k) = 0;
          else
              matrix_oss(u,k) = sscanf(lin(2+16*(k-1):16*k-2),'%f');
              if (k == col_L1) & (numel(lin)>=16*k)
                  snr(u,1) = sscanf(lin(16*k),'%f');
              end
          end
      end
   end

else

   %more than 5 observation types --> 2 lines to be read
   matrix_oss = matrix_oss(:,[1 2 3 4 5]);
   num_type_oss = 5;
   for u = 1:num_sat
      lin = fgetl(file_RINEX);
      lin_doppler = fgetl(file_RINEX); %#ok<NASGU>
      for k = 1:num_type_oss
         if (length(lin) < 1+16*(k-1)) | (isempty(sscanf(lin(1+16*(k-1):16*k-2),'%f')))
             matrix_oss(u,k) = 0;
         else
             %data save
             matrix_oss(u,k) = sscanf(lin(1+16*(k-1):16*k-2),'%f');
             if (k == col_L1) & (numel(lin)>=16*k)
                 snr(u,1) = sscanf(lin(16*k),'%f');
             end
         end
      end
   end

end
