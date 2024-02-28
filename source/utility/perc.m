function perc_val = perc(data, p)
% SYNTAX:
%   perc_val = perc(data, thr)
%
% INPUT:
%   data    array of values
%   p       percentile requested (0 1]
%
% OUTPUT:
%   perc_val
%
% DESCRIPTION:
%   returns percentile of the values in data

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
    data(isnan(data)) = [];
if size(data,1) == 1 || size(data,2) == 1
    data = sort(data(:));
else
        for c = 1:size(data,2)
            data(:,c) = sort(data(:,c));
        end
    end
    if numel(data) > 0
        for c = 1:size(data,2)
            for v = 1:numel(p)
                perc_val(c,v) = data(min(max(round(size(data,1) * p(v)),1), size(data,1)),c);
end
        end
    else
        perc_val = 0;
    end
end
