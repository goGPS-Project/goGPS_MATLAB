function [ x_par, Cxx  ] = LSsolver(x0,y0,xm,q);
   
% normal matrix
N = [sum(q), sum(q .* (x0-xm)); sum(q .* (x0-xm )), sum((q .* (x0-xm).^2))];

% known term
Y = [sum(q .* y0); sum(q .* y0 .* (x0-xm))];

% estimated parameters
x_par = inv(N) * Y;

% Cxx = sigma02 * inv(N); % covariance matrix