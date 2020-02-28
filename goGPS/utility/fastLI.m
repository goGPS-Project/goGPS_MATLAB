function y_pred = fastLI(y_obs, x_pred)
% SYNTAX:
%   y_pred = fastLI(y_obs, x_pred);
%
% INPUT:
%   y_obs  = row-vector of [1 x p_deg + 1] data values
%   y_pred = row-vector of [1 x numel(x_pred)], where interpolation is to be found (could be a single value)
%
% OUTPUT:
%   y_pred = a row-vector of interpolated y-values
%
% DESCRIPTION:
%   Lagrange interpolation algorithm supposing y_obs regularly sampled
%   The degree of the polinomial (p_deg) is equivalent to the number of
%   element of y_obs - 1
%
% ORIGINAL CODE:
%   Author: Dmitry Pelinovsky
%   Available online at: http://dmpeli.mcmaster.ca/Matlab/Math4Q3/Lecture2-1/LagrangeInter.m
%
% SEE ALSO:
%   LagrangeInter

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b7
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Author of the previous implementation: Dmitry Pelinovsky
%  Available online at: http://dmpeli.mcmaster.ca/Matlab/Math4Q3/Lecture2-1/LagrangeInter.m
%  Written by:       Andrea Gatti
%  Contributors:     ...
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

    n_obs = size(y_obs,2);
    n_pred = numel(x_pred);

    %y_pred = zeros(size(x_pred));
    tmp = ones(n_obs, n_pred);
    for i = 1 : n_pred
        for k = 1 : n_obs
            %y_pred(i) = y_obs * prod((x_pred(i) - subtractor - x_pred(i) * eye(11)) ./ denum, 2);
            for kk = 1 : (k-1) % start the inner loop through the data values for x (if k = 0 this loop is not executed)
                tmp(kk, i) = tmp(kk, i) .* ((x_pred(i) - k) / (kk - k)); % see the Lagrange interpolating polynomials
                %L(kk+1,:) = L(kk+1,:).*(xi - x(k+1))/(x(kk+1)-x(k+1)); % see the Lagrange interpolating polynomials
            end % end of the inner loop

            for kk = k+1 : n_obs % start the inner loop through the data values (if k = n this loop is not executed)
                tmp(kk, i) = tmp(kk, i) .* ((x_pred(i) - k) / (kk - k));
            end % end of the inner loop
        end
    end
    y_pred = y_obs * tmp;
end
