function [flag, prn, IODe, delta_x, delta_y, delta_z, delta_offset] = v_code0(half_msg, iodp_mask, prn_mask)

% SYNTAX:
%   [flag, prn, IODe, delta_x, delta_y, delta_z, delta_offset] = v_code0(half_msg, iodp_mask, prn_mask);
%
% INPUT:
%   half_msg = binary string from ems2ltc.m
%   iodp_mask = PRN mask IODPs [vector]
%   prn_mask  = PRN masks [matrix]
%
% OUTPUT:
%   flag = flag data/no data (see note 2 p.25 App.A RTCA DO 229C) [vector]
%   prn  = PRN of the satellites corresponding to the delta corrections [vector]
%   IODe = ephemeris IOD [vector]
%   delta_x = delta X ECEF of the SVs in prn [vector]
%   delta_y = delta Y ECEF of the SVs in prn [vector]
%   delta_z = delta Z ECEF of the SVs in prn [vector]
%   delta_offset = delta offset of the clock of the SVs in prn [vector]
%
% DESCRIPTION:
%   

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

flag    = NaN(1,2);
prn     = NaN(1,2); 
IODe    = NaN(1,2); 
delta_x = NaN(1,2); 
delta_y = NaN(1,2); 
delta_z = NaN(1,2); 
delta_offset = NaN(1,2);

%resolution of the ECEF coordinates
res_ecef = 2^(-3); %LSB value = 0.125 (m)

%resolution of the clock offset
res_offset = 2^(-31); %LSB value (s)

%in the message each delta ECEF is 9 bits
param_8 = 2^(8) - 1; %2^(n_bit-1) - 1

% for twos complement
param_9 = 2^(9) - 1; %2^(n_bit-1) - 1

% further split the half messages in two parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1st part

%PRN mask number
prn_n = fbin2dec(half_msg(1:6));

if (prn_n == 0)
    
    flag(1)    = 1;   %null message
    prn(1)     = NaN;
    IODe(1)    = NaN;
    delta_x(1) = NaN;
    delta_y(1) = NaN;
    delta_z(1) = NaN;
    delta_offset(1) = NaN;
    
else
    
    flag(1) = 0; %valid message
    
	iodp = fbin2dec(half_msg(103:104));

    %find the right mask
    i_iodp = find(iodp_mask == iodp);
    i_prn_mask = prn_mask(i_iodp,:); %#ok<FNDSB>

    %PRN SV
    prn(1) = i_prn_mask(prn_n);

    IODe(1) = fbin2dec(half_msg(7:14));
    
    %delta position
    d_x = fbin2dec(half_msg(15:23));
    d_y = fbin2dec(half_msg(24:32));
    d_z = fbin2dec(half_msg(33:41));
    
    vect_d = [d_x d_y d_z];
    
    %twos complement
    p8 = vect_d > param_8;
    vect_d = vect_d - 2^(9) * p8;
    
    %multiply by the LSB value
    vect_delta = vect_d * res_ecef; % m 
    delta_x(1) = vect_delta(1);
    delta_y(1) = vect_delta(2);
    delta_z(1) = vect_delta(3);

    %clock offset
    d_offset = fbin2dec(half_msg(42:51));
    
    %twos complement
    p9 = d_offset > param_9;
    d_offset = d_offset - 2^(10) * p9;
    
    %multiply by the LSB value
    delta_offset(1) = d_offset * res_offset; % s
    
end

%2nd part

%PRN mask number
prn_n = fbin2dec(half_msg(52:57));

if (prn_n == 0)
    
    flag(2)    = 1;   %null message
    prn(2)     = NaN;
    IODe(2)    = NaN;
    delta_x(2) = NaN;
    delta_y(2) = NaN;
    delta_z(2) = NaN;
    delta_offset(2) = NaN;
    
else
    
    flag(2) = 0; %valid message
    
% 	iodp = fbin2dec(half_msg(103:104));
% 
%     %trovo la mask corretta
%     i_iodp = find(iodp_mask == iodp);
%     i_prn_mask = prn_mask(i_iodp,:);

    %PRN SV
    prn(2) = i_prn_mask(prn_n);

    IODe(2) = fbin2dec(half_msg(58:65));
    
    %delta position
    d_x = fbin2dec(half_msg(66:74));
    d_y = fbin2dec(half_msg(75:83));
    d_z = fbin2dec(half_msg(84:92));
    
    vect_d = [d_x d_y d_z];
    
    %twos complement
    p8 = vect_d > param_8;
    vect_d = vect_d - 2^(9) * p8;
    
    %multiply by LSB value
    vect_delta = vect_d * res_ecef; % m 
    delta_x(2) = vect_delta(1);
    delta_y(2) = vect_delta(2);
    delta_z(2) = vect_delta(3);

    %clock offset
    d_offset = fbin2dec(half_msg(93:102));
    
    %twos complement
    p9 = d_offset > param_9;
    d_offset = d_offset - 2^(10) * p9;
    
    %multiply by LSB value
    delta_offset(2) = d_offset * res_offset; % s
 
end
