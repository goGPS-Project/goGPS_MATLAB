function checksum = checksumFTX (bit_msg, len)

% initialization
checksum = 0;
pos = 1;

for i = 0 : (len-1)
    tmp = (checksum + 1) * (bin2dec([bit_msg(pos+8:pos+15) bit_msg(pos:pos+7)])+i);
    sh_tmp = bitshift(tmp,-16);
    tmp = bitxor(tmp, sh_tmp);
    checksum = bitxor(checksum, tmp);
    % optimization
    if (checksum > 65535)
        checksum = checksum - (floor(checksum/65536))*65536;
%         checksum = dec2bin(checksum,16);
%         checksum = bin2dec(checksum(end-15:end));
    end
    pos = pos+16;
end

checksum = dec2bin(checksum,16);
checksum = [checksum(9:16) checksum(1:8)];

% ---- Test ----
% dec2hex(bin2dec(checksum))
