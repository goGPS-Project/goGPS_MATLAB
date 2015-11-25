function [MT, flag, prn, IODe, delta_x, delta_y, delta_z, delta_offset] = ems2ltc(msg, iodp_mask, prn_mask)

% SYNTAX:
%   [MT, flag, prn, IODe, delta_x, delta_y, delta_z, delta_offset] = ems2ltc(msg, iodp_mask, prn_mask);
%
% INPUT:
%   msg = hexadecimal string for 1 EGNOS message (from .ems files)
%         eg: msg = ('9A6430A00C0DFE01101EFE7DFE5FF2169415FF817FB4C8419F30F0168FC02D40');
%   iodp_mask = PRN mask IODPs [vector]
%   prn_mask  = PRN masks [matrix]
%
% OUTPUT:
%   MT = long term correction message types
%   flag = flag data/no data (see note 2 p.25 App.A RTCA DO 229C) [vector]
%   prn  = PRN of the satellites corresponding to the delta corrections [vector]
%   IODe = ephemeris IOD [vector]
%   delta_x = delta X ECEF of the SVs in prn [vector]
%   delta_y = delta Y ECEF of the SVs in prn [vector]
%   delta_z = delta Z ECEF of the SVs in prn [vector]
%   delta_offset = delta offset of the clock of the SVs in prn [vector]
%
% DESCRIPTION:
%   Extract the long term corrections from MT 24 and 25.

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

%MT24: msg = ('C660000003FE8004000003BA5BB9903A640FFF0400D04704FF504FF4BB729AC0');
%MT25: msg = ('9A6430A00C0DFE01101EFE7DFE5FF2169415FF817FB4C8419F30F0168FC02D40');

flag    = NaN(1,4);
prn     = NaN(1,4); 
IODe    = NaN(1,4); 
delta_x = NaN(1,4); 
delta_y = NaN(1,4); 
delta_z = NaN(1,4); 
delta_offset = NaN(1,4);

%EGNOS message length: 250 bits
% usually in .ems files there are 64 hexadecimal values (256 bits), thus
% the last 6 bits are not used

h2d = hex2dec(msg');         %decimal to hexadecimal conversion
d2b = dec2bin(h2d, 4);       %4-bit binary number
s   = reshape(d2b', 1, 256); %binary string representing all the 256 bits of the message

%start the conversion from binary to decimal
MT = fbin2dec(s(9:14));

%distinguish MT 25 and 24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (MT == 25)

   % 1st half

   vel_code1 = fbin2dec(s(15));
   
   %half message, 105 bits long. NOT 106 since we are not considering the velocity bit
   half_mess1 = s(16:120);
   
   if (vel_code1 == 0)
       
       [flag1, prn1, IODe1, delta_x1, delta_y1, delta_z1, delta_offset1] = v_code0(half_mess1, iodp_mask, prn_mask);
       
       flag(1,1:2)    = flag1;
       prn(1,1:2)     = prn1;
       IODe(1,1:2)    = IODe1;
       delta_x(1,1:2) = delta_x1;
       delta_y(1,1:2) = delta_y1;
       delta_z(1,1:2) = delta_z1;
       delta_offset(1,1:2) = delta_offset1;
       
%    elseif (vel_code1 == 1)
%        
%        [f] = v_code1(half_mess1, iodp_mask, prn_mask);
%        
%        flag(1,1)    = flag1;
%        prn(1,1)     = prn1;
%        IODe(1,1)    = IODe1;
%        delta_x(1,1) = delta_x1;
%        delta_y(1,1) = delta_y1;
%        delta_z(1,1) = delta_z1;
%        delta_offset(1,1) = delta_offset1;
       
   end

   % 2nd half

   vel_code2 = fbin2dec(s(121));
   
   %half message, 105 bits long. NOT 106 since we are not considering the velocity bit
   half_mess2 = s(122:226);
   
   if vel_code2 == 0
       
       [flag2, prn2, IODe2, delta_x2, delta_y2, delta_z2, delta_offset2] = v_code0(half_mess2, iodp_mask, prn_mask);
       
       flag(1,3:4)    = flag2;
       prn(1,3:4)     = prn2;
       IODe(1,3:4)    = IODe2;
       delta_x(1,3:4) = delta_x2;
       delta_y(1,3:4) = delta_y2;
       delta_z(1,3:4) = delta_z2;
       delta_offset(1,3:4) = delta_offset2;
       
%    elseif vel_code2 == 1
%        
%        [f] = v_code1(half_mess2, iodp_mask, prn_mask);
%        
%        flag(1,3)    = flag2;
%        prn(1,3)     = prn2;
%        IODe(1,3)    = IODe2;
%        delta_x(1,3) = delta_x2;
%        delta_y(1,3) = delta_y2;
%        delta_z(1,3) = delta_z2;
%        delta_offset(1,3) = delta_offset2;
       
   end

elseif (MT == 24)

   % ONLY 2nd half

   vel_code1 = fbin2dec(s(121));
   
   %half message, 105 bits long. NOT 106 since we are not considering the velocity bit
   half_mess1 = s(122:226);
   
   if (vel_code1 == 0)
       
       [flag1, prn1, IODe1, delta_x1, delta_y1, delta_z1, delta_offset1] = v_code0(half_mess1, iodp_mask, prn_mask);
       
       flag(1,1:2)    = flag1;
       prn(1,1:2)     = prn1;
       IODe(1,1:2)    = IODe1;
       delta_x(1,1:2) = delta_x1;
       delta_y(1,1:2) = delta_y1;
       delta_z(1,1:2) = delta_z1;
       delta_offset(1,1:2) = delta_offset1;
       
%    elseif (vel_code1 == 1)
%        
%        [f] = v_code1(half_mess1, iodp_mask, prn_mask);
%        
%        flag(1,1)    = flag1;
%        prn(1,1)     = prn1;
%        IODe(1,1)    = IODe1;
%        delta_x(1,1) = delta_x1;
%        delta_y(1,1) = delta_y1;
%        delta_z(1,1) = delta_z1;
%        delta_offset(1,1) = delta_offset1;
       
   end
end
