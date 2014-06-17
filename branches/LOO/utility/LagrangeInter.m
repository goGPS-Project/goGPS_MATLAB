function yi = LagrangeInter(x,y,xi)

% SYNTAX:
%   yi = LagrangeInter(x,y,xi);
%
% INPUT:
%   x,y = row-vectors of (n+1) data values (x,y)
%   xi  = a row-vector of x-values, where interpolation is to be found (could be a single value)
%
% OUTPUT:
%   yi = a row-vector of interpolated y-values
%
% DESCRIPTION:
%   Lagrange interpolation algorithm.

%----------------------------------------------------------------------------------------------
%
% Author: Dmitry Pelinovsky
% Available online at: http://dmpeli.mcmaster.ca/Matlab/Math4Q3/Lecture2-1/LagrangeInter.m
%
%----------------------------------------------------------------------------------------------

n = length(x) - 1; % the degree of interpolation polynomial
ni = length(xi); % the number of x-values, where interpolation is to be found

L = ones(n+1,ni); % the matrix for Lagrange interpolating polynomials L_(n,k)(x)
                    % has (n+1) rows for each polynomial at k = 0,1,...,n
                    % has ni column for each x-value of xi
                    
% Note: the algorithm uses the MATLAB capacities for matrix handling!
% The two nested loops below are designed to compute in parallel 
% the values of Lagrange interpolating polynomials at each x-value of xi !

for k = 0 : n % start the outer loop through the data values for x
    
    for kk = 0 : (k-1) % start the inner loop through the data values for x (if k = 0 this loop is not executed)
        L(kk+1,:) = L(kk+1,:).*(xi - x(k+1))/(x(kk+1)-x(k+1)); % see the Lagrange interpolating polynomials
    end % end of the inner loop
    
    for kk = k+1 : n % start the inner loop through the data values (if k = n this loop is not executed)
        L(kk+1,:) = L(kk+1,:).*(xi - x(k+1))/(x(kk+1)-x(k+1)); 
    end % end of the inner loop
    
end % the end of the outer loop

% Now multiply the values for Lagrange interpolating polynomials by the data values for y
yi = y * L; % use matrix multiplication of row-vector y and the matrix L, the output is the row-vector
