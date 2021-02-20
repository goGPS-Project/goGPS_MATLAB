function [ avg_out ] = nan_mean( data )
% computue mean ignoring nan
%
% SYNTAX
%  [ std ] = nan_std( data )
avg_out = mean(noNaN(data));
end
