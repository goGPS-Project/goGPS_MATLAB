function [MT, IODF, prc, udrei, sv, n_sv] = ems2prc(msg, iodp_mask, prn_mask)

% SYNTAX:
%   [MT, IODF, prc, udrei, sv, n_sv] = ems2prc(msg, iodp_mask, prn_mask);
%
% INPUT:
%   msg = hexadecimal string for 1 EGNOS message (from .ems files)
%         eg: msg = ('C60C7FD0000000003FCC003FC8000003FB4007FE0029BBBBB9BB9599F32C2A40');
%   iodp_mask = PRN mask IODPs [vector]
%   prn_mask  = PRN masks [matrix]
%
% OUTPUT:
%   MT = fast corrections (FC) message types
%   IODF  = quality flag (0,1,2: OK, 3: problem with one or more SVs)
%   prc   = pseudorange corrections [vector]
%   udrei = UDREI (User Differential Range Error Indicator) [vector]
%   sv    = PRN of the satellites corresponding to the PRCs
%   n_sv  = number of previous SVs
%
% DESCRIPTION:
%   PRC from EGNOS messages MT 0, 2, 3, 4, 5, 24.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Code contributed by Giuliano Sironi, 2011
% Portions of code contributed by Antonio Herrera Olmo, 2012
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

%msg = ('C60C7FD0000000003FCC003FC8000003FB4007FE0029BBBBB9BB9599F32C2A40');
%msg = ('530D7FD0000000003FCC003FC8000003FB4007FE4029BBBBB9BB9599F4D92900');

%resolution of EGNOS differential corrections
res = 2^(-3); % LSB value = 0.125 (m)
% NOTE: differential corrections range: da -256.000 a +255.875 m

%each PRC is given by 12 bits
param_11 = 2^(11) - 1; %2^(n_bit-1) - 1
% integer to signed integer in twos complement

%EGNOS message length: 250 bits
% usually in .ems files there are 64 hexadecimal values (256 bits), thus
% the last 6 bits are not used

h2d = hex2dec(msg');         %decimal to hexadecimal conversion
d2b = dec2bin(h2d, 4);       %4-bit binary number
s   = reshape(d2b', 1, 256); %binary string representing all the 256 bits of the message

%start the conversion from binary to decimal
MT = fbin2dec(s(9:14));

%distinguish MT 0, 2, 3 and 4 from MT 24
if (MT == 2 | MT == 3 | MT == 4 | MT == 0)

    IODF = fbin2dec(s(15:16));
    %NOTE:
    % IODF = 0,1,2: OK
    %      = 3:     problem with one or more SVs

    IODP = fbin2dec(s(17:18));
    
    %find the right PRN mask
    i_iodp = find(iodp_mask == IODP);
    i_prn_mask = prn_mask(i_iodp,:); %#ok<FNDSB>
   
    %compute the PRC
    %start bit and end bit of the 13 PRCs in MT 2, 3, 4
    start_bit = [19 : 12 : 163];
    end_bit   = [30 : 12 : 174];
    
    %allocate the 13 binary numbers that correspond to the PRCs
    for j = 1 : 13
        num_bin(j,:) = s(start_bit(j) : end_bit(j));
    end

    %convert to integers
    num_int = fbin2dec(num_bin);
    
    %convert in twos complement
    p11 = num_int > param_11;
    num_int = num_int - 2^(12) * p11;
    
    %multiply by the LSB value
    PRC = num_int * res; % m 

    %compute the UDREI
    %start bit and end bit of the 13 UDREIs in MT 2, 3, 4
    start_bit_u = [175 : 4 : 223];
    end_bit_u   = [178 : 4 : 226];
    
    %allocate the 13 binary numbers that correspond to the UDREIs
    for j = 1 : 13
        num_bin_u(j,:) = s(start_bit_u(j) : end_bit_u(j));
    end

    %convert to integers
    UDREI = fbin2dec(num_bin_u);
    % NOTE:
    % UDREI = 0 - 13: OK
    %       = 14:     SV non monitored
    %       = 15:     SV must not be used
    
    if MT == 2 | MT == 0
        SV = [1 : 13];
    elseif MT == 3
        SV = [14 : 26];
        elseif MT == 4
        SV = [27 : 39];
        if max(prn_mask) < 39
            SV = [27 : max(prn_mask)];  
        end
    end

elseif (MT == 24)

    %comute the PRCs
    %start bit and end bit of the 6 PRCs in MT 24
    start_bit = [15 : 12 : 75];
    end_bit   = [26 : 12 : 86];
    
    %allocate the 13 binary numbers that correspond to the PRCs
    for j = 1 : 6
        num_bin(j,:) = s(start_bit(j) : end_bit(j));
    end

    %convert to integers
    num_int = fbin2dec(num_bin);
    
    %twos complement
    p11 = num_int > param_11;
    num_int = num_int - 2^(12) * p11;
    
    %multiply by the LSB value
    PRC = num_int * res; % m

    %compute the UDREI
    %start bit and end bit of the 6 UDREIs in MT 24
    start_bit_u = [87 : 4 : 107];
    end_bit_u   = [90 : 4 : 110];
    
    %allocate the 6 binary numbers that correspond to the UDREIs
    for j = 1 : 6
        num_bin_u(j,:) = s(start_bit_u(j) : end_bit_u(j));
    end

    %convert to integers
    UDREI = fbin2dec(num_bin_u);
    % NOTE:
    % UDREI = 0 - 13: OK
    %       = 14:     SV non monitored
    %       = 15:     SV must not be used
    
    IODP = fbin2dec(s(111:112));
    
    %find the right PRN mask
    i_iodp = find(iodp_mask == IODP);
    i_prn_mask = prn_mask(i_iodp,:); %#ok<FNDSB>
    
    block_ID = fbin2dec(s(113:114));
      
    SV_block = [1 : 6; ...
                14 : 19; ...
                27 : 32];
   
    SV = SV_block(block_ID + 1, :);         
    
    IODF = fbin2dec(s(115:116));
    %NOTE:
    % IODF = 0,1,2: OK
    %      = 3:     problem with one or more SVs
end

%select only the valid daata (UDREI < 14)
%data_ok = find(UDREI < 14);

prc   = NaN(1, 13);
udrei = NaN(1, 13);
sv    = NaN(1, 13);

prc(1 : length(PRC))     = PRC;
udrei(1 : length(UDREI)) = UDREI;

sv(1 : length(SV)) = i_prn_mask(SV);
n_sv               = length(SV);
