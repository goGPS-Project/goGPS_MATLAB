function [ std_out ] = nan_std( data )
% computue std ignoring nan
%
% SYNTAX
%  [ std ] = nan_std( data )
std_out = std(noNaN(data));
end

