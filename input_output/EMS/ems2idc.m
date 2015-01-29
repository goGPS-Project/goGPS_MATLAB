function [band, block, ivd, givei, igp] = ems2idc(msg, iodi_mask, band_mask, igp_mask)

% SYNTAX:
%   [band, block, ivd, givei, igp] = ems2idc(msg, iodi_mask, band_mask, igp_mask);
%
% INPUT:
%   msg = hexadecimal string for 1 EGNOS message (from .ems files)
%         eg: msg = ('536A53FDFFEF05C02DC16E0A7053829C14FFEFFF7FFBFFDFFEFFF7803BEF9280');
%   iodi_mask  = mask IODI
%   band_mask = mask band
%   igp_mask   = mask IGP
%
% OUTPUT:
%   band  = output band
%   block = output block
%   ivd   = ionospheric vertical delays
%   givei = output givei
%   igp   = output IGP
%
% DESCRIPTION:
%   Find the ionospheric vertical delays.

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

%msg = ('536A53FDFFEF05C02DC16E0A7053829C14FFEFFF7FFBFFDFFEFFF7803BEF9280');

%resolution of the EGNOS IGP vertical delays
risoluzione = 2^(-3); %LSB value = 0.125 (m)
% NOTE: range from 0 to +63.875 m

h2d = hex2dec(msg');         %decimal to hexadecimal conversion
d2b = dec2bin(h2d, 4);       %4-bit binary number
s   = reshape(d2b', 1, 256); %binary string representing all the 256 bits of the message

%start the conversion from binary to decimal
%MT = fbin2dec(s(9:14));

band = fbin2dec(s(15:18));

block = fbin2dec(s(19:22));

%Ionospheric vertical delays
%start bit and end bit of the 15 ivd
start_bit = [23 : 13 : 205];
end_bit   = [31 : 13 : 213];
    
%allocate the 15 binary numbers that correspond to the iono vertical delays
for j = 1 : 15
    num_bin(j,:) = s(start_bit(j) : end_bit(j));
end

%convert to decimal
num_int = fbin2dec(num_bin);

%multiply by the LSB value
ivd = num_int' * risoluzione; % m

%GIVEI
%start bit and end bit of the 15 GIVEIs
start_bit_g = [32 : 13 : 214];
end_bit_g   = [35 : 13 : 217];
    
%allocate the 15 binary numbers that correspond to the GIVEIs
for j = 1 : 15
    num_bin_g(j,:) = s(start_bit_g(j) : end_bit_g(j));
end

%convert to decimal
givei = fbin2dec(num_bin_g)';

IODI = fbin2dec(s(218:219));

%find the right IGP mask
i_band = find(band_mask == band);
    
i_iodi = find(iodi_mask(i_band) == IODI);
igp = igp_mask(i_band(i_iodi),:); %#ok<FNDSB>
