% logFit Fit the data (x, y) using a logarithmic function + constant
% and evaluates it at x_out.
%
% Parameters:
%   x: Input data for x
%   y: Input data for y
%   x_out: Points at which the fitted function should be evaluated. If
%          not provided, it defaults to a linspace from 1 to max(x).
%
% Returns:
%   line_fit: Evaluated values of the fitted function at x_out

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

function line_fit = logFit(x, y, x_out)
    % Initial fit for threshold estimation
    params_initial = fitLogModel(x, y);
    residuals = y - (params_initial(1) * log(params_initial(2) * x + 1) + params_initial(3));
    residuals = residuals(~isnan(residuals) & ~isinf(residuals)); % Remove NaNs and Infs
    threshold = 1.4826 * mad(residuals); % Median Absolute Deviation

    % RANSAC fit
    [bestParams, ~] = ransacLogFit(x, y, 1000, threshold);
    
    % Predict using the best parameters from RANSAC
    line_fit = bestParams(1) * log(bestParams(2) * x_out + 1) + bestParams(3);
end

function [bestParams, bestInliers] = ransacLogFit(x, y, iterNum, threshold)
    maxInliers = 0;
    bestParams = [];
    bestInliers = [];

    for i = 1:iterNum
        % Randomly select subset of data
        indices = randperm(length(x), 3);  % 3 is the minimum number of points to fit this model
        subsetX = x(indices);
        subsetY = y(indices);
        
        % Fit the model to this subset
        params = fitLogModel(subsetX, subsetY);
        
        % Calculate the errors for the entire dataset
        predictedY = params(1) * log(params(2) * x + 1) + params(3);
        errors = abs(predictedY - y);
        
        % Find inliers
        inliers = find(errors < threshold);
        
        % Check if this model is better than previous best
        if length(inliers) > maxInliers
            maxInliers = length(inliers);
            bestParams = params;
            bestInliers = inliers;
        end
    end
end

function params = fitLogModel(x, y)
    % Define the model: y = a * log(b * x + 1) + c
    logModel = @(params, x) params(1) * log(params(2) * x + 1) + params(3);

    % Initial guess based on the data (you can refine this as needed)
    A_init = (max(y) - min(y)) / log(max(x));
    B_init = 1 / max(x);
    C_init = min(y);
    initial_params = [A_init, B_init, C_init];
    
    % Bounds (to ensure positive values for a and b and no bounds for c)
    lb = [0, 0, -Inf];
    ub = [Inf, Inf, Inf];

    % Use lsqcurvefit to estimate the parameters
    params = lsqcurvefit(logModel, initial_params, x, y, lb, ub);
end

function params = lsqcurvefit(fun, x0, xdata, ydata, lb, ub)
    % Parameters for the algorithm
    maxIter = 1000;
    tol = 1e-6;
    lambda = 0.001;
    nu = 2.0; % Factor for increasing/decreasing lambda

    params = x0;
    n = length(params);
    I = eye(n);

    r = ydata - fun(params, xdata);
    prev_cost = r' * r;

    for iter = 1:maxIter
        % Compute Jacobian
        J = jacobian(fun, params, xdata, ydata);
        
        while true
            % Adjust and solve the normal equations
            delta = (J' * J + lambda * I) \ (J' * r);
            
            % Test update
            new_params = params + delta';
            new_params = max(new_params, lb);
            new_params = min(new_params, ub);
            
            % Compute new residuals
            new_r = ydata - fun(new_params, xdata);
            new_cost = new_r' * new_r;
            
            if new_cost < prev_cost
                % If cost reduction, accept the step and reduce lambda
                params = new_params;
                prev_cost = new_cost;
                lambda = lambda / nu;
                break;
            else
                % If no cost reduction, increase lambda and try again
                lambda = lambda * nu;
                if lambda > 1e10
                    return;
                end
            end
        end
        
        % Gradient-based stopping criterion
        grad = J' * r;
        if max(abs(grad)) < tol
            break;
        end
    end
end


function J = jacobian(fun, params, xdata, ydata)
    % Compute Jacobian using central differences
    eps = 1e-4; % A small value
    f0 = fun(params, xdata);
    J = zeros(length(f0), length(params));
    
    for i = 1:length(params)
        params1 = params;
        params2 = params;
        params1(i) = params1(i) - eps;
        params2(i) = params2(i) + eps;
        
        f1 = fun(params1, xdata);
        f2 = fun(params2, xdata);
        
        J(:, i) = (f2 - f1) / (2 * eps);
    end
end

function m = mad(x)
    % mad calculates the median absolute deviation of an array, ignoring NaNs.
    % It's implemented to avoid dependence on the Statistics and Machine Learning Toolbox.
    %x = x(~isnan(x));
    %m = median(abs(x - median(x)));
    % faster but risky:
    omitnan = false;
    m = matlab.internal.math.columnmedian(abs(x(:) - matlab.internal.math.columnmedian(x,omitnan)),omitnan);
end

