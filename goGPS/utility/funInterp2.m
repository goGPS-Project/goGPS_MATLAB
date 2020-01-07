function data = funInterp2(x_pred, y_pred, x_obs, y_obs, data_obs, fun)
%
% SINTAX:
%   data = funInterp2(x_pred, y_pred, x_obs, y_obs, data_obs, fun)
%
% INPUT:
%   x_pred , y_pred     coordinates of the prediction point (arrays [n x 1])
%   x_obs, y_obs        coordinates of the observation (array [n x 1], sparse points)
%   data_obs            data observation array [n x n_epochs]
%   fun                 virtual function f(dist) as function of the distance
%
% OUTPUT:
%   data                interpolated data in the point of interest [n x n_epochs]
%
% DEFAULT VALUES:
%   generic interpolator using as correlation function fun
%   <default: fun = @(dist) exp(-dist)>
%
% EXAMPLE:
%   fun = @(dist) 0.2 * exp(-(dist/1e4)) + exp(-(dist/6e3).^2);
%   temp = funInterp2(ep(:), np(:), e_obs(:), n_obs(:), td_obs(:,:), fun);
%

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti ...
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

    narginchk(5, 6);

    if nargin < 6
        % Correlation function
        fun = @(dist) exp(-dist);
    end

    % Init out data
    data = nan(size(data_obs, 2), numel(x_pred));

    [x_mesh, y_mesh] = meshgrid(x_obs, y_obs);
    d_obs = sqrt(abs(x_mesh - x_mesh').^2 + abs(y_mesh - y_mesh').^2);
    q_fun_obs = fun(d_obs);
    %q_fun_obs = exp(-d_obs/0.4e4);
    [xv, ~, xi] = unique(x_pred);
    [yv, ~, yi] = unique(y_pred);
    x2 = (repmat(x_obs, 1, numel(xv)) - repmat(xv', numel(x_obs), 1)).^2;
    y2 = (repmat(y_obs, 1, numel(yv)) - repmat(yv', numel(y_obs), 1)).^2;
    for i = 1 : numel(x_pred)
        d_pred = sqrt(x2(:,xi(i)) + y2(:,yi(i)));
        c_mat = q_fun_obs .* repmat(fun(d_pred)', size(q_fun_obs,1),1);
        %c_mat = bsxfun(@times, q_fun_obs, fun(d_pred)');
        c_mat = triu(c_mat) + triu(c_mat, 1)';

        trans = sum(c_mat);
        w = trans / sum(trans);
        data(:, i) = (w * data_obs)';
    end
end
