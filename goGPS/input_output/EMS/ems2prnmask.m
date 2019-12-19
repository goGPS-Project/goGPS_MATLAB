function [iodp_mask, prn_mask] = ems2prnmask(msg)

% SYNTAX:
%   [iodp_mask, prn_mask] = ems2prnmask(msg);
%
% INPUT:
%   msg = hexadecimal string for 1 EGNOS message (from .ems files)
%         eg: msg = ('9A05FFFFFFF8000000000000000000000440000000000000000000007853ED40');
%
% OUTPUT:
%   iodp_mask = MT1 IODP
%   prn_mask  = MT1 PRN mask [vector]
%
% DESCRIPTION:
%   Decode the PRN mask from EGNOS MT 1.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Giuliano Sironi 2011
%  Contributors:     Giuliano Sironi 2011, ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

% the PRN mask and the IODP mask connect:
%   -the PRC for each SV ID of the MT 0, 2, 3, 4 and 24
%   -the long term corrections of MT 24 and 25
% to the right GPS satellite.
% The association may change when a new SV is launched / an old SV is deactivated

%msg = ('9A05FFFFFFF8000000000000000000000440000000000000000000007853ED40');

%EGNOS MT1 message structure (prn mask):
% 8 bit: preamble
% 6 bit: message type (MT)
% 32 bit: bit mask for GPS satellites, ordered from PRN 1 to 32
%   se bit = 0 -> SV to be disregarded
%   se bit = 1 -> SV to be used
%       es: 01111111111111111111111111111110
%       PRN 1  NOT available (bit 1 = 0)
%       PRN 2  available     (bit 2 = 1)
%       ...
%       PRN 31 available     (bit 31 = 1)
%       PRN 32 NOT available (bit 32 = 0)
% 180 bit: bit mask for other satellites (eg: GLONASS)
% 24 bit: parity
% TOTAL: 250 bit

%EGNOS message length: 250 bits
% usually in .ems files there are 64 hexadecimal values (256 bits), thus
% the last 6 bits are not used

h2d = hex2dec(msg');         %decimal to hexadecimal conversion
d2b = dec2bin(h2d, 4);       %4-bit binary number
s   = reshape(d2b', 1, 256); %binary string representing all the 256 bits of the message

%start the conversion from binary to decimal
% MT = fbin2dec(s(9:14)); %must be = 1

%bits that compose the GPS satellites mask (PRN from 1 to 32)
mask = fbin2dec(s(15 : 15 + 31)')'; %vector with the mask bits

%vector with the PRNs of the GPS SVs
%sv_gps = [1 : 32];

%keep only the GPS SVs with valid mask (=1)
%prn_mask = find (mask .* sv_gps);
prn_mask = find (mask);

%vector the PRNs of active GPS satellites
prn_mask = [prn_mask, zeros(1, 32-length(prn_mask))];

%IODP
iodp_mask = fbin2dec(s(225:226));
