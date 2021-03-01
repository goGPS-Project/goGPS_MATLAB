function [data] = simpleFill2D(data, flags, fun)
% SYNTAX:
%    [data_filled] = simpleFill1D(data, flags)
%
% DESCRIPTION:
%    fill flagged data with a simple interpolation using MATLAB
%    interp1 'pchip', 'extrap'
%
% NOTE: data can be a matrix, the operation is executed column by column
%

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 Geomatics Research & Development srl (GRed)
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti
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

if nargin < 3
        % Correlation function
        fun = @(dist) exp(-(dist).^2);
end

    x_g = 1 : size(data,2);
    y_g = 1 : size(data,1);
    
    [x_mesh, y_mesh] = meshgrid(x_g, y_g);
    
    x_p = x_mesh(:);
    x_o = x_mesh(~flags);
    
    y_p = y_mesh(:);
    y_o = y_mesh(~flags);
    
    data_o = data(~flags);
    
    for x = 1 : size(data,2)
        id_p = (x_p == x);
        if sum(id_p) > 0
            tmp_x_p = x_p(id_p);
            tmp_y_p = y_p(id_p);
            d = sqrt((repmat(tmp_x_p, 1, size(x_o,1)) - repmat(x_o', size(tmp_x_p, 1), 1)).^2 + ...
                (repmat(tmp_y_p, 1, size(y_o,1)) - repmat(y_o', size(tmp_y_p, 1), 1)).^2);
            %d = d / sqrt(sum(size(data).^2));
            data(id_p) = ((fun(d) * data_o) ./ sum(fun(d),2));
        end            
    end
end
