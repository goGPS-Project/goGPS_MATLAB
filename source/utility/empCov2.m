function [emp_cov, dist] = empCov2(x_obs, y_obs, data_obs, n_classes)
% SYNTAX:
%   [emp_cov, dist] = empCov2(x_obs, y_obs, data_obs, n_classes)
%
% DESCRIPTION:
%   return the array of the empirical covariance funtion
%
% EXAMPLE:
%   [emp_cov, dist] = empCov2(e_obs(id_td), n_obs(id_td), td_obs);
%   figure; plot(dist, emp_cov, 'b'); hold on;

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 Geomatics Research & Development srl (GReD)
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

    if nargin < 4
        n_classes = 50;
    end

    [x_mesh, y_mesh] = meshgrid(x_obs, y_obs);
    d_obs = sqrt(abs(x_mesh - x_mesh').^2 + abs(y_mesh - y_mesh').^2);
    classes = ceil(d_obs / max(d_obs(:)) * (n_classes-1))+1;

    emp_cov = zeros(max(classes(:)), size(data_obs,2));
    for i = 1 : size(data_obs,2)
        corr = (data_obs(:,i)-mean(data_obs(:,i))) * (data_obs(:,i)-mean(data_obs(:,i)))';
        for c = 1 : max(classes)
            emp_cov(c, i) = mean(corr(serialize(triu(classes)==c)));
        end
    end
    emp_cov = mean(emp_cov, 2);
    dist = (0 : (n_classes - 1))' * max(d_obs(:))/n_classes;
end
