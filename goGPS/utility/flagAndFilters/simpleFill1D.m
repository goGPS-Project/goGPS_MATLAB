function [data] = simpleFill1D(data, flags, method)
% SYNTAX:
%    [data_filled] = simpleFill1D(data, flags)
%
% DESCRIPTION:
%    fill flagged data with a simple interpolation using MATLAB
%    interp1 'pchip', 'extrap'
%
% NOTE: data can be a matrix, the operation is executed column bby column
%

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.2 beta 1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
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

%     if nargin == 2
%         half_win = 6;
%     end
%
%     lim = getOutliers(flags);
%     for j = 1 : size(lim, 1)
%         dt_lhs = data(max(1, lim(j, 1) - half_win) : lim(j, 1) - 1);
%         dt_rhs = data(lim(j, 2) + 1 : min(numel(data), lim(j, 2) + half_win));
%         data(lim(j,1):lim(j,2)) = interp1([ lim(j, 1) - numel(dt_lhs) : lim(j, 1) - 1,  lim(j, 2) + 1 : lim(j, 2) + numel(dt_rhs)], [dt_lhs; dt_rhs], lim(j,1) : lim(j,2), 'pchip','extrap');
%     end
    narginchk(2,3);
    if nargin == 2
        method = 'pchip';
    end
    t = 1 : size(data, 1);
    for r = 1 : size(data, 2)
        if any(~isnan(data(:, r))) && any(~isnan(flags(:, r)))
            jmp = find(flags(:, r));
            flags(:, r) = flags(:, r) | isnan(data(:, r));
            data(jmp, r) = interp1(t(~flags(:, r)), data(~flags(:, r), r), jmp, method,'extrap');
        end
    end
end
