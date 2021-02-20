function [ frac_part ] = fracFNI( num )
%return fractional distance to nearest integer
%
% SYNTAX:
% [ frac_part ] = fracFNI( num )
frac_part = num - round(num);
end