function [ frac_part ] = frac( num )
%return fractional part of a number
%
% SYNTAX:
% [ frac_part ] = frac( num )
frac_part = num - floor(num);
end

