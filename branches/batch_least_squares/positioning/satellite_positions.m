function [XS, dtS, XS_tx, VS_tx, time_tx, no_eph, sys_idx] = satellite_positions(time_rx, pseudorange, sat, Eph, SP3, sbas, err_tropo, err_iono, dtR)

% SYNTAX:
%   [XS, dtS, XS_tx, VS_tx, time_tx, no_eph, sys_idx] = satellite_positions(time_rx, pseudorange, sat, Eph, SP3, sbas, err_tropo, err_iono, dtR);
%
% INPUT:
%   time_rx     = reception time
%   pseudorange = observed code pseudoranges
%   sat         = available satellite indexes
%   Eph         = ephemeris
%   SP3         = structure containing precise ephemeris data
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
%   sys_idx = array with different values for different systems
%
% DESCRIPTION:

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

nsat = length(sat);

time_tx = zeros(nsat,1);
dtS     = zeros(nsat,1);
XS      = zeros(nsat,3);
XS_tx   = zeros(nsat,3);
VS_tx   = zeros(nsat,3);

%satellites with no ephemeris available
no_eph  = zeros(nsat,1);

%system array
sys_idx = zeros(nsat,1);

for i = 1 : nsat
    
    k = find_eph(Eph, sat(i), time_rx);
    
    if (isempty(k) & isempty(SP3))
        no_eph(i) = 1;
        continue
    end

    %compute signal transmission time
    [time_tx(i,1), dtS(i,1)] = transmission_time(time_rx, pseudorange(i), sat(i), Eph(:,k), SP3, sbas, err_tropo(i), err_iono(i), dtR);

    if (isempty(time_tx(i,1)))
        no_eph(i) = 1;
        continue
    end
    
    if (isempty(SP3))
        
        %compute satellite position and velocity at transmission time
        [XS_tx(i,:), VS_tx(i,:)] = satellite_orbits(time_tx(i,1), Eph(:,k), sat(i), sbas);
        
        %detect satellite constellation
        sys = Eph(31,k);
    else
        %interpolate SP3 coordinates at transmission time
        [XS_tx(i,:), VS_tx(i,:)] = interpolate_SP3_coord(time_tx(i,1), SP3, sat(i));

        %relativistic correction term
        dtrel = relativistic_error_correction(time_tx(i,1), Eph, XS_tx(i,:), VS_tx(i,:));
        time_tx(i,1) = time_tx(i,1) - dtrel;
        dtS(i,1) = dtS(i,1) + dtrel;
        
        %second iteration for taking into account the relativistic effect
        [XS_tx(i,:), VS_tx(i,:)] = interpolate_SP3_coord(time_tx(i,1), SP3, sat(i));
        
        %detect satellite constellation
        sys = SP3.sys(sat(i));
    end
    
    %computation of ECEF satellite position at time_rx
    traveltime = time_rx - time_tx(i,1);
    switch char(sys)
        case 'G'
            Omegae_dot = goGNSS.OMEGAE_DOT_GPS; sys_idx(i,1) = 1;
        case 'R'
            Omegae_dot = goGNSS.OMEGAE_DOT_GLO; sys_idx(i,1) = 2;
        case 'E'
            Omegae_dot = goGNSS.OMEGAE_DOT_GAL; sys_idx(i,1) = 3;
        case 'C'
            Omegae_dot = goGNSS.OMEGAE_DOT_BDS; sys_idx(i,1) = 4;
        case 'J'
            Omegae_dot = goGNSS.OMEGAE_DOT_QZS; sys_idx(i,1) = 5;
        otherwise
            fprintf('Something went wrong in satellite_positions.m\nUnrecognized Satellite system!\n');
            Omegae_dot = goGNSS.OMEGAE_DOT_GPS;
    end
    XS(i,:) = earth_rotation_correction(traveltime, XS_tx(i,:), Omegae_dot);
end
