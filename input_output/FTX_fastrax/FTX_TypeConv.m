function [dec_out, pos_out, part] = FTX_TypeConv(type_in, bin_msg, pos_in)

%
% Code contributed by Ivan Reguzzoni

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
case 'INT16'    % 16 bits   (range: –32768 through 32767 decimal)
        temp1 = bin_msg(pos_in:pos_in+7);  pos_in = pos_in + 8;
        temp2 = bin_msg(pos_in:pos_in+7);  pos_in = pos_in + 8;
        dec_out = [temp2,temp1];
        dec_out = twos_complement_inside(dec_out);
        part(1) = bin2dec(temp1); part(2) = bin2dec(temp2);
case 'INT32'    % 32 bits   (range: –2147483648 through 2147483647 decimal)
        temp1 = bin_msg(pos_in:pos_in+7);  pos_in = pos_in + 8;
        temp2 = bin_msg(pos_in:pos_in+7);  pos_in = pos_in + 8;
        temp3 = bin_msg(pos_in:pos_in+7);  pos_in = pos_in + 8;
        temp4 = bin_msg(pos_in:pos_in+7);  pos_in = pos_in + 8;
        dec_out = [temp4,temp3,temp2,temp1];
        dec_out = twos_complement_inside(dec_out);
        part(1) = bin2dec(temp1); part(2) = bin2dec(temp2);
        part(3) = bin2dec(temp3); part(4) = bin2dec(temp4);
case 'WORD'     % 16 bits   (range: 0 through 65535 decimal)
        part(1) = bin2dec(bin_msg(pos_in:pos_in+7));  pos_in = pos_in + 8;
        part(2) = bin2dec(bin_msg(pos_in:pos_in+7));  pos_in = pos_in + 8;
        dec_out = part(1) + (part(2) * 2^8);  % little endian
case 'DWORD'    % 32 bits   (range: 0 through 4294967295 decimal)
        part(1) = bin2dec(bin_msg(pos_in:pos_in+7));  pos_in = pos_in + 8;
        part(2) = bin2dec(bin_msg(pos_in:pos_in+7));  pos_in = pos_in + 8;
        part(3) = bin2dec(bin_msg(pos_in:pos_in+7));  pos_in = pos_in + 8;
        part(4) = bin2dec(bin_msg(pos_in:pos_in+7));  pos_in = pos_in + 8;
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
        out = (-1) * (bin2dec(in) + 1);
    else
        out = bin2dec(in);
    end    
end