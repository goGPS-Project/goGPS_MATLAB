function [n_bands, band, iodi, igp_mask] = ems2igpmask(msg)

% SYNTAX:
%   [n_bands, band, iodi, igp_mask] = ems2igpmask(msg);
%
% INPUT:
%   msg = hexadecimal string for 1 EGNOS message (from .ems files)
%         eg: msg = ('534964000000FFFFFFC00000001FFFC00001FFFC00003FFFE001FC003C15ABC0');
%
% OUTPUT:
%   n_bands = number of transmitted bands
%   band = band number
%   iodi = band IODI
%   igp_mask = IGP mask
%
% DESCRIPTION:
%   The IGP mask (MT 18) connects the ionospheric delay for each node
%   of the MT 26 to the right IGP (Iono Grid Point).

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

%msg = ('534964000000FFFFFFC00000001FFFC00001FFFC00003FFFE001FC003C15ABC0');

%EGNOS message length: 250 bits
% usually in .ems files there are 64 hexadecimal values (256 bits), thus
% the last 6 bits are not used

h2d = hex2dec(msg');         %decimal to hexadecimal conversion
d2b = dec2bin(h2d, 4);       %4-bit binary number
s   = reshape(d2b', 1, 256); %binary string representing all the 256 bits of the message

%start the conversion from binary to decimal
% MT = fbin2dec(s(9:14)); %must be = 18

%number of bands
n_bands = fbin2dec(s(15:18)); %range: 0 - 11

%band number
band = fbin2dec(s(19:22)); %range: 0 - 10

%issue of data ionospheric
iodi = fbin2dec(s(23:24)); %range: 0 - 3

%bit mask for the IGP of the corresponding band
mask = fbin2dec(s(25 : 25 + 200)')'; %vector with the mask bits

%IGP vector
igp = [1 : 201];       

%%keep only the IGPs with valid mask (=1)
igp_mask = find (mask .* igp);  

igp_mask = [igp_mask, zeros(1, 201-length(igp_mask))];
