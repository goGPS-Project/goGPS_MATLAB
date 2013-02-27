function [XS, dtS, XS_tx, VS_tx, time_tx, no_eph] = satellite_positions(time_rx, pseudorange, sat, Eph, SP3_time, SP3_coor, SP3_clck, sbas, err_tropo, err_iono, dtR)

% SYNTAX:
%   [XS, dtS, XS_tx, VS_tx, time_tx, no_eph] = satellite_positions(time_rx, pseudorange, sat, Eph, SP3_time, SP3_coor, SP3_clck, sbas, err_tropo, err_iono, dtR);
%
% INPUT:
%   time_rx     = reception time
%   pseudorange = observed code pseudoranges
%   sat         = available satellite PRNs
%   Eph         = ephemeris
%   SP3_time    = precise ephemeris time
%   SP3_coor    = precise ephemeris coordinates
%   SP3_clck    = precise ephemeris clocks
%   sbas        = SBAS corrections
%   err_tropo   = tropospheric delays
%   err_iono    = ionospheric delays
%   dtR         = receiver clock offset
%
% OUTPUT:
%   XS      = satellite position at transmission time in ECEF(time_rx) (X,Y,Z)
%   dtS     = satellite clock error (vector)
%   XS_tx   = satellite position at transmission time in ECEF(time_tx) (X,Y,Z)
%   VS_tx   = satellite velocity at transmission time in ECEF(time_tx) (X,Y,Z)
%   time_tx = transmission time (vector)
%   no_eph  = satellites with no ephemeris available (vector) (0: available, 1: not available)
%
% DESCRIPTION:

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

nsat = length(sat);

time_tx = zeros(nsat,1);
dtS     = zeros(nsat,1);
XS      = zeros(nsat,3);
XS_tx   = zeros(nsat,3);
VS_tx   = zeros(nsat,3);

%satellites with no ephemeris available
no_eph  = zeros(nsat,1);

for i = 1 : nsat
    
    k = find_eph(Eph, sat(i), time_rx);
    
    if (isempty(k) & isempty(SP3_time))
        no_eph(i) = 1;
        continue
    end
    
    %compute signal transmission time
    [time_tx(i,1), dtS(i,1)] = transmission_time(time_rx, pseudorange(i), sat(i), Eph(:,k), SP3_time, SP3_clck, sbas, err_tropo(i), err_iono(i), dtR);
    
    if (dtS(i,1) == 0)
        no_eph(i) = 1;
        continue
    end
    
    if (isempty(SP3_time))
        
        %compute satellite position and velocity at transmission time
        [XS_tx(i,:), VS_tx(i,:)] = satellite_orbits(time_tx(i,1), Eph(:,k), sat(i), sbas);
    else
        %interpolate SP3 coordinates at transmission time
        [XS_tx(i,:), VS_tx(i,:)] = interpolate_SP3_coord(time_tx(i,1), SP3_time, SP3_coor(:,sat(i),:));

        %relativistic correction term
        dtrel = relativistic_error_correction(time_tx(i,1), Eph, XS_tx(i,:), VS_tx(i,:));
        time_tx(i,1) = time_tx(i,1) - dtrel;
        
        %second iteration for taking into account the relativistic effect
        [XS_tx(i,:), VS_tx(i,:)] = interpolate_SP3_coord(time_tx(i,1), SP3_time, SP3_coor(:,sat(i),:));
    end
    
    %computation of ECEF satellite position at time_rx
    traveltime = time_rx - time_tx(i,1);
    XS(i,:) = earth_rotation_correction(traveltime, XS_tx(i,:));
end
