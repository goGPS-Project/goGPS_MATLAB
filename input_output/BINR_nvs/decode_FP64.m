function [out, flag] = decode_FP64(in)

flag = 1;

sign = str2num(in(1));
esp  = fbin2dec(in(2:12));
mant = fbin2dec(in(13:64)) / 2^52;

if (esp > 0 && esp < 2047)
    out = (-1)^sign * (2^(esp - 1023)) * (1 + mant);
elseif (esp == 0 && mant ~= 0)
    out = (-1)^sign * (2^(-1022)) * mant;
elseif (esp == 0 && mant == 0)
    out = (-1)^sign * 0;
elseif (esp == 2047 && mant == 0)
    out = (-1)^sign * Inf;
    flag = 0;
    keyboard
elseif (esp == 2047 && mant ~= 0)
    out = NaN;
    flag = 0;
    keyboard
end