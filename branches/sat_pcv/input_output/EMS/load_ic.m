function [igp, m_ivd, lat_igp, lon_igp, GPS_time26] = load_ic(iodi_mask, band_mask, igp_mask, MT, msg, GPS_time)

% SYNTAX:
%   [igp, m_ivd, lat_igp, lon_igp, GPS_time26] = load_ic(iodi_mask, band_mask, igp_mask, MT, msg, GPS_time);
%
% INPUT:
%   iodi_mask = IGP mask IODI [vector]
%   band_mask = band numbers [vector]
%   igp_mask  = IGP masks of the transmitted bands [matrix]
%   MT  = message types [vector]
%   msg = EGNOS message hexadecimal strings [matrix]
%   GPS_time = GPS time of the EGNOS messages [matrix]
%
% OUTPUT:
%   igp   = IDs of the IGP nodes of the global grid [vector]
%   m_ivd = ionosphere vertical delays at IGP nodes [matrix]
%   lat_igp  = IGP node latitude values [vector]
%   lon_igp  = IGP node longitude values [vector]
%   GPS_time = GPS time for the m_ivd matrix [matrix]
%
% DESCRIPTION:
%   Ionospheric vertical delay correction.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Code contributed by Giuliano Sironi, 2011
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

%load the global grid coordinates
[grid_igp] = global_IGP_grid;

%band * 1000 + IGP_node_number
grid_id = grid_igp(:,1) * 1000 + grid_igp(:,2);

%keep only the MTs that contain iono vertical delays
% WARNING: messages are not provided at regular intervals
r_MT = find(MT == 26);

%number of MT 26
n26 = length(r_MT);

%initialization
band      = NaN(n26,1);
block     = NaN(n26,1);
ivd_old   = NaN(n26,15);
givei_old = NaN(n26,15);
%IODI     = NaN(n26,1);
nodes_igp = NaN(n26,15);
nodes_id  = NaN(n26,15);
ivd       = NaN(n26,15);
givei     = NaN(n26,15);
n_nodes    = NaN(n26,1);

%GPS time of the messages with iono vertical delays
GPS_time26 = GPS_time(r_MT,:);

for i = 1 : n26
    
    [banda26, block26, ivd26, givei26, igp_temp] = ems2idc(msg(r_MT(i),:), iodi_mask, band_mask, igp_mask);
    
    band(i)        = banda26;
    block(i)       = block26;
    ivd_old(i,:)   = ivd26;     %IGP vertical delay [1x15]
    givei_old(i,:) = givei26;   %[1x15]
    igp_mask26     = igp_temp';
    %IODI(i)       = IODI26;

    %band and block values <--> IGP nodes of the IGP mask
    %r_b = find(bands == band(i));
    %igp_mask26 = igp_mask(:,r_b);
    
    %block is from 0 to 13, thus nodes_mask has max value = 210, but the
    %mask max value is 201. However it should not be a problem, since only
    %the grid nodes over Europe are provided, thus the mask is never "full"
    %and not all block for the band are available
    nodes_mask = [1:15] + 15 * block(i);
    
    %NOTE: some nodes of a block may not be valid, and they must be removed
    %(e.g. I have 9 valid IGP points (IGP>0), but the block is always of 15
    %nodes)
    r_im = find(igp_mask26(nodes_mask));
    nodes_igp(i,1:length(r_im)) = igp_mask26(nodes_mask(r_im));
    
    n_nodes(i) = length(r_im);
    
    % >>> new ID method for IGP nodes <<<
    % from band (1 digit) and node number (3 digits) to only one number of
    % 4 digits; the first digit indicates the band, the other 3 the node.
    % In this way the search is done on one dataset only.
    nodes_id(i,1:length(r_im)) = band(i) * 1000 + nodes_igp(i,1:length(r_im));
    
    %IGP matrix of the iono vertical delays
    ivd(i,1:length(r_im)) = ivd_old(i,1:length(r_im));
    
    %givei matrix
    givei(i,1:length(r_im)) = givei_old(i,1:length(r_im));
end    


%find all the valid IGP points
igp = unique(nodes_id);

%remove NaNs
igp = igp(~isnan(igp))';

%ivd matrix allocation
m_ivd = NaN(n26, length(igp));

%initialize the ivd matrix
for k = 1 : n_nodes(1)
    i_igp = find(igp == nodes_id(1, k));
    if givei(1,k) < 15 & ivd(1,k) < 63.875 %"dont'use" value
        m_ivd(1, i_igp) = ivd(1,k); %#ok<FNDSB>
    end
end

%number of elements ~= NaN on each line
%n_nigp(1) = sum(~isnan(m_ivd(1,:)));

for i = 2 : n26
    m_ivd(i,:) = m_ivd(i-1,:);
    for j = 1 : n_nodes(i)
        i_igp = find(igp == nodes_id(i, j));
        if givei(i,j) < 15 & ivd(i,j) < 63.875 %"dont'use" value
            m_ivd(i, i_igp) = ivd(i,j);
        else 
            m_ivd(i, i_igp) = NaN;
        end
    end
    
    %number of elements ~= NaN on each line
    %n_nigp(i) = sum(~isnan(m_ivd(i,:)));
end

%NOTE: sometimes the rows of matrix m_ivd are duplicated because the new message
%has no effect on the values of the iono vertical delays on the valid IGPs

%latitude and longitude vectors for the IGP nodes
lat_igp = NaN(1, length(igp));
lon_igp = NaN(1, length(igp));

for i = 1 : length(igp)
    r_lat = find(grid_id == igp(i));
    lat_igp(i) = grid_igp(r_lat, 3);
    lon_igp(i) = grid_igp(r_lat, 4);
end
