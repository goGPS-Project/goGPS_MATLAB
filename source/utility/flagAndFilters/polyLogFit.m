% Fits data to a polynomial of a specified degree plus a logarithmic term.
%
% [fitFunc, coeffs, x_out] = polyLogFit(x, y, deg, x_out) fits the input 
% data (x,y) to a polynomial of degree 'deg' plus an additional logarithmic 
% term. The output fitFunc is the fitted function value computed at the points 
% x_out. If x_out is not provided, it defaults to a logarithmic spacing 
% from the minimum to the maximum of x with 100 points. The function returns 
% the coefficients of the fitted polynomial plus the logarithmic term in 
% 'coeffs'.
%
% Inputs:
% - x: Vector of independent data points.
% - y: Vector of dependent data points.
% - deg: Degree of polynomial to fit.
% - x_out: (optional) Vector of points at which to compute the fitted values.
%
% Outputs:
% - fitFunc: Vector of fitted function values computed at the points in x_out.
% - coeffs: Vector of fitted coefficients for the polynomial and the logarithmic term.
% - x_out: Vector of points at which the fitted values were computed.
%
% Example:
% [y_fit, p, x_vals] = polyLogFit([1 2 3 4], [1.1 2.9 6.2 10.1], 2);
%
% Notes:
% 1. The function fits data to a polynomial of the form: 
%    a_n*x^n + a_{n-1}*x^{n-1} + ... + a_1*x + a_0 + b*log(x), where 'n' 
%    is the degree specified and 'b' is the coefficient of the logarithmic term.
% 2. The coefficients are returned in descending order of power of 'x', 
%    followed by the logarithmic term coefficient.
%

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0
%
%--------------------------------------------------------------------------
%  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti, ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

function [fitFunc, coeffs, x_out] = polyLogFit(x, y, deg, x_out)
    % Function to fit data to a polynomial of degree 'deg' plus a logarithmic term
    
    if nargin < 4
        x_out = logspace(log10(min(x)), log10(max(x)), 100);
    end

    % Create a matrix with columns [x^n, x^{n-1}, ..., x, 1, log(x)]
    X = zeros(length(x), deg+2);
    for i = 1:deg
        X(:,i) = x.^(deg+1-i);
    end
    X(:,end-1) = ones(size(x));
    X(:,end) = log(x);

    % Solve for the coefficients using a least-squares approach
    coeffs = X \ y;

    % Generate the fit for the provided x_out values
    X_out = zeros(length(x_out), deg+2);
    for i = 1:deg
        X_out(:,i) = x_out.^(deg+1-i);
    end
    X_out(:,end-1) = ones(size(x_out));
    X_out(:,end) = log(x_out);
    fitFunc = X_out * coeffs;
end
