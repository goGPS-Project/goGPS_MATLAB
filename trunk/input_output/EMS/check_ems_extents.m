function [ems_data_available] = check_ems_extents(time_R, pr, snr, Eph, iono, sbas)

% SYNTAX:
%   [ems_data_available] = check_ems_extents(time_R, pr, snr, Eph, iono, sbas);
%
% INPUT:
%   time_R = reference vector of GPS time of week
%   pr     = pseudorange
%   snr    = signal-to-noise ratio
%   Eph    = ephemerides
%   iono   = ionospheric parameters (Klobuchar)
%   sbas   = SBAS corrections
%
% OUTPUT:
%   ems_data_available = boolean flag for data availability check
%
% DESCRIPTION:
%   Function that check that the approximate position of the receiver
%   (first available positioning epoch) is within the EMS grids.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.1 beta
%
% Copyright (C) 2009-2012 Mirko Reguzzoni, Eugenio Realini
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

ems_data_available = 0;

fprintf('Checking that the receiver approximate position falls within the available EMS grids... ');

pos_R = zeros(3,1);

if (~isempty(find(Eph(1,:,:) ~= 0, 1)))
    
    cutoff = 15;
    snr_threshold = 0;
    
    i = 1;
    
    while (sum(abs((pos_R))) == 0 & i <= length(time_R))
        
        satObs = find(pr(:,i) ~= 0);
        
        Eph_t  = rt_find_eph (Eph, time_R(i));
        
        satEph = find(Eph_t(1,:) ~= 0);
        satAvail = intersect(satObs,satEph)';
        
        if (length(satAvail) >=4)
            pos_R = init_positioning(time_R(i), pr(satAvail,i), snr(satAvail,i), Eph_t(:,:), [], [], [], iono(:,i), [], [], [], [], satAvail, cutoff, snr_threshold, 0, 0);
        end
        
        i = i + 1;
        
    end
end

if (sum(abs((pos_R))) ~= 0)
    
    [lat_R, lon_R] = cart2geod(pos_R(1), pos_R(2), pos_R(3));
    
    igp4 = sel_igp(lat_R, lon_R, sbas.igp, sbas.lat_igp, sbas.lon_igp);
    
    if(isempty(igp4))
        fprintf('FALSE\n');
    else
        fprintf('TRUE\n');
        fprintf('EMS files successfully read. Applying SBAS corrections.\n');
        ems_data_available = 1;
    end
else
    fprintf('\n');
    fprintf('Positioning not possible. Processing stopped.\n');
    return
end
