function [dec_out, pos_out, part] = FTX_TypeConv(type_in, bin_msg, pos_in)

% SYNTAX:
%   [dec_out, pos_out, part] = FTX_TypeConv(type_in, bin_msg, pos_in);
%
% INPUT:
%   type_in = data type in input [string]
%   bin_msg = binary message [string]
%   pos_in  = starting position in the binary message
%
% OUTPUT:
%   dec_out = data output in decimal format
%   pos_out = ending position in the binary message
%   part    = bytes composing the data output in decimal format
%
% DESCRIPTION:
%   Fastrax ITALK data types conversion tool.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Code contributed by Ivan Reguzzoni
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

if nargin < 2 || isempty(type_in) || isempty(bin_msg)
    error('stats:corr:TooFewInputs', ...
          'Requires other data.');
end

if ( ischar(type_in) == 0 )
    error('First element is not a string!')
end

if ( ischar(bin_msg) == 0 )
    error('Binary message must be a string type')
end

if (nargin == 2) && ischar(type_in)
    pos_in = 1;
end

switch type_in
case 'BOOL'
        % No definition
        error('Type not valid.');
case 'INT16'    % 16 bits   (range: ?2768 through 32767 decimal)
        temp1 = bin_msg(pos_in:pos_in+7);  pos_in = pos_in + 8;
        temp2 = bin_msg(pos_in:pos_in+7);  pos_in = pos_in + 8;
        dec_out = [temp2,temp1];
        dec_out = twos_complement_inside(dec_out);
        part(1) = fbin2dec(temp1); part(2) = fbin2dec(temp2);
case 'INT32'    % 32 bits   (range: ?147483648 through 2147483647 decimal)
        temp1 = bin_msg(pos_in:pos_in+7);  pos_in = pos_in + 8;
        temp2 = bin_msg(pos_in:pos_in+7);  pos_in = pos_in + 8;
        temp3 = bin_msg(pos_in:pos_in+7);  pos_in = pos_in + 8;
        temp4 = bin_msg(pos_in:pos_in+7);  pos_in = pos_in + 8;
        dec_out = [temp4,temp3,temp2,temp1];
        dec_out = twos_complement_inside(dec_out);
        part(1) = fbin2dec(temp1); part(2) = fbin2dec(temp2);
        part(3) = fbin2dec(temp3); part(4) = fbin2dec(temp4);
case 'WORD'     % 16 bits   (range: 0 through 65535 decimal)
        part(1) = fbin2dec(bin_msg(pos_in:pos_in+7));  pos_in = pos_in + 8;
        part(2) = fbin2dec(bin_msg(pos_in:pos_in+7));  pos_in = pos_in + 8;
        dec_out = part(1) + (part(2) * 2^8);  % little endian
case 'DWORD'    % 32 bits   (range: 0 through 4294967295 decimal)
        part(1) = fbin2dec(bin_msg(pos_in:pos_in+7));  pos_in = pos_in + 8;
        part(2) = fbin2dec(bin_msg(pos_in:pos_in+7));  pos_in = pos_in + 8;
        part(3) = fbin2dec(bin_msg(pos_in:pos_in+7));  pos_in = pos_in + 8;
        part(4) = fbin2dec(bin_msg(pos_in:pos_in+7));  pos_in = pos_in + 8;
        dec_out = part(1) + (part(2) * 2^8) + (part(3) * 2^16) + (part(4) * 2^24);  % little endian
case 'FLOAT'    % 32 bits
        % Original C Code 
        % INT16 mantissa = pW[0];
        % INT16 exponent = pW[1];
        % return mantissa ? (float)ldexp(mantissa, exponent - 15) : 0.0f;
        % ---------------
        Mantissa1 = bin_msg(pos_in:pos_in+7);  pos_in = pos_in + 8;
        Mantissa2 = bin_msg(pos_in:pos_in+7);  pos_in = pos_in + 8;
        Mantissa  = [Mantissa2,Mantissa1];
        Exponent1 = bin_msg(pos_in:pos_in+7);  pos_in = pos_in + 8;
        Exponent2 = bin_msg(pos_in:pos_in+7);  pos_in = pos_in + 8;
        Exponent  = [Exponent2,Exponent1];
        
        Mantissa = twos_complement_inside(Mantissa);
        Exponent = twos_complement_inside(Exponent);

        dec_out = Mantissa * 2^ (Exponent-15);
        part(1) = Mantissa;
        part(2) = Exponent;
        clear Mantissa1 Mantissa2 Exponent1 Exponent2 Mantissa Exponent        
case 'DOUBLE'   % 48 bits
        % Original C Code 
        % INT32 mantissa = MAKELONG(pW[0], pW[1]);
        % INT16 exponent = pW[2];
        % return mantissa ? ldexp(mantissa, exponent - 31) : 0;
        % ---------------
        Mantissa1 = bin_msg(pos_in:pos_in+7);  pos_in = pos_in + 8;
        Mantissa2 = bin_msg(pos_in:pos_in+7);  pos_in = pos_in + 8;
        Mantissa3 = bin_msg(pos_in:pos_in+7);  pos_in = pos_in + 8;
        Mantissa4 = bin_msg(pos_in:pos_in+7);  pos_in = pos_in + 8;
        Mantissa  = [Mantissa4,Mantissa3,Mantissa2,Mantissa1];
        Exponent1 = bin_msg(pos_in:pos_in+7);  pos_in = pos_in + 8;
        Exponent2 = bin_msg(pos_in:pos_in+7);  pos_in = pos_in + 8;
        Exponent  = [Exponent2,Exponent1];
        
        Mantissa = twos_complement_inside(Mantissa);
        Exponent = twos_complement_inside(Exponent);

        dec_out = Mantissa * 2^ (Exponent-31);
        part(1) = Mantissa;
        part(2) = Exponent;
        clear Mantissa1 Mantissa2 Mantissa3 Mantissa4 Exponent1 Exponent2 Mantissa Exponent
otherwise
    error('Type not valid.');
end

pos_out = pos_in;

end
%--------------------------------------------------------------------------

function out = twos_complement_inside(in)
    if in(1) == '1'
        in(:) = num2str(not(str2num(in(:))));         % logical not
        out = (-1) * (fbin2dec(in) + 1);
    else
        out = fbin2dec(in);
    end    
end