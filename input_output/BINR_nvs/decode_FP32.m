function [out, flag] = decode_FP32(in)

flag = 1;

sign = str2num(in(1));
esp  = fbin2dec(in(2:9));
mant = fbin2dec(in(10:32)) / 2^23;

if (esp > 0 && esp < 255)
    out = (-1)^sign * (2^(esp - 127)) * (1 + mant);
elseif (esp == 0 && mant ~= 0)
    out = (-1)^sign * (2^(-126)) * mant;
elseif (esp == 0 && mant == 0)
    out = (-1)^sign * 0;
elseif (esp == 255 && mant == 0)
    out = (-1)^sign * Inf;
    flag = 0;
    keyboard
elseif (esp == 255 && mant ~= 0)
    out = NaN;
    flag = 0;
    keyboard
end