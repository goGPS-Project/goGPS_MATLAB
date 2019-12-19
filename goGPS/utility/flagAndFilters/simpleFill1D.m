function [data] = simpleFill1D(data, flags, method, time)
% Fill in flags position
%
% SYNTAX:
%    [data_filled] = simpleFill1D(data, flags, <method>)
%
% DESCRIPTION:
%    fill flagged data with a simple interpolation using MATLAB
%    interp1 'linear', 'pchip' (default), 'extrap', 'last'
%
%   'last' is not a standard interpolation method, it uses the last
%   valid epoch of an arc to fill all the consecutive nan values
%
% NOTE: data can be a matrix, the operation is executed column by column
%

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
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
    narginchk(2,4);
    if nargin == 2 || isempty(method)
        method = 'pchip';
    end

    if strcmp(method, 'last')
        for s = 1 : size(data, 2)
            lim = getOutliers(flags(:,s));
            % if the data starts with nan do nothing
            if ~isempty(lim) && lim(1) == 1
                lim(1, :) = [];
            end
            for l = 1 : size(lim, 1)
                data(lim(l,1) : lim(l,2), s) = data(lim(l,1) - 1, s);
            end
        end
    else
        if nargin == 4 && ~isempty(time)
            t = time;
        else
            t = 1 : size(data, 1);
        end
        for r = 1 : size(data, 2)
            if any(~isnan(data(:, r))) && any(~isnan(flags(:, r))) && any(~flags(:, r))
                jmp = find(flags(:, r));
                flags(:, r) = flags(:, r) | isnan(data(:, r));
                if sum(~flags(:, r)) > 1
                    data(jmp, r) = interp1(t(~flags(:, r)), data(~flags(:, r), r), t(jmp), method, 'extrap');
                end
            end
        end
    end
end
