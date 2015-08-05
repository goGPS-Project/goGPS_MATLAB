function [dx_E, dy_E, dz_E, doffset_E, iode_E, GPS_time25] = load_ltc(iodp_mask, prn_mask, MT, msg, GPS_time)

% SYNTAX:
%   [dx_E, dy_E, dz_E, doffset_E, iode_E, GPS_time25] = load_ltc(iodp_mask, prn_mask, MT, msg, GPS_time);
%
% INPUT:
%   iodp_mask = IODP masks [vector]
%   prn_mask  = PRN masks [matrix]
%   MT  = message types [vector]
%   msg = EGNOS message strings [matrix]
%   GPS_time = GPS time of the messages
%
% OUTPUT:
%   dx_E = X ECEF satellite position corrections [matrix]
%   dy_E = Y ECEF satellite position corrections [matrix]
%   dz_E = Z ECEF satellite position corrections [matrix]
%   doffset_E = satellite clock corrections [matrix]
%   iode_E = IODE [matrix]
%   GPS_time25 = GPS time of the correction matrices
%
% DESCRIPTION:
%   Load the long term corrections (LTC).

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
        
%keep only the MTs that contain long term corrections
% WARNING: the messages are not provided at regular time intervals
r_MT = find(MT == 25 | MT == 24);

%number of message types
n25 = length(r_MT);

%number of GPS satellites
nGPSsat = 32;

%initialization
flag = NaN(n25,4);
prn  = NaN(n25,4);
IODe = NaN(n25,4);
d_x  = NaN(n25,4);
d_y  = NaN(n25,4);
d_z  = NaN(n25,4);
d_offset = NaN(n25,4);

nsv = NaN(n25,1);

%GPS time of the messages
GPS_time25 = GPS_time(r_MT,:);

for i = 1 : n25
    
    [mt, flag25, prn25, IODe25, d_x25, d_y25, d_z25, d_offset25] = ems2ltc(msg(r_MT(i),:), iodp_mask, prn_mask); %#ok<ASGLU>
    
    %MT25(i)       = mt;
    flag(i,:)     = flag25;
    prn(i,:)      = prn25;
    IODe(i,:)     = IODe25;
    d_x(i,:)      = d_x25;
    d_y(i,:)      = d_y25;
    d_z(i,:)      = d_z25;
    d_offset(i,:) = d_offset25;
    
    nsv(i) = length(find(~isnan(prn(i,:)))); %max 4 SV
end

%gather the data for all satellites
dx_E = zeros(n25, nGPSsat);
dy_E = zeros(n25, nGPSsat);
dz_E = zeros(n25, nGPSsat);
doffset_E = zeros(n25, nGPSsat);
iode_E = zeros(n25, nGPSsat);

i_sv = find(~isnan(prn(1,:)));
for k = 1 : nsv(1)
                
     %i_dx = find(SV == prn(1, i_sv(k)));
        
     %store data ordered by PRN (column-wise)
     dx_E(1, prn(1, i_sv(k))) = d_x(1,i_sv(k));
     dy_E(1, prn(1, i_sv(k))) = d_y(1,i_sv(k));
     dz_E(1, prn(1, i_sv(k))) = d_z(1,i_sv(k));
     doffset_E(1, prn(1, i_sv(k))) = d_offset(1,i_sv(k));
     iode_E(1, prn(1, i_sv(k))) = IODe(1,i_sv(k));
           
end

%store the corrections for the subsequent epochs
for j = 2 : n25
    
    dx_E(j,:) = dx_E(j-1,:);
    dy_E(j,:) = dy_E(j-1,:);
    dz_E(j,:) = dz_E(j-1,:);
    doffset_E(j,:) = doffset_E(j-1,:);
    iode_E(j,:) = iode_E(j-1,:);
    
    i_sv = find(~isnan(prn(j,:)));
    
    for k = 1 : nsv(j)
                
        %i_dx = find(SV == prn(j, i_sv(k)));
        
        %store data ordered by PRN (column-wise)
        dx_E(j, prn(j, i_sv(k))) = d_x(j,i_sv(k));
        dy_E(j, prn(j, i_sv(k))) = d_y(j,i_sv(k));
        dz_E(j, prn(j, i_sv(k))) = d_z(j,i_sv(k));
        doffset_E(j, prn(j, i_sv(k))) = d_offset(j,i_sv(k)); 
        iode_E(j, prn(j, i_sv(k))) = IODe(j,i_sv(k));
    end
end
