function detrended_data = strongDeTrend(data, thr_perc, thr_perc_global, n_sigma)
% Returns the mean per column with a certaint percentile of values in data
% The code uses only the data filtered by the requested percentile to estimate 
% an std column by column
% With this estimation only the data within n_sigma range are used for the 
% robust mean estimation
%
% SYNTAX:
%   smean = strongMean(data, thr_perc, thr_perc_global, n_sigma)
%
% INPUT:
%   data                matrix of values
%   thr_perc            percentile requested [0 1] per column
%   thr_perc_global     percentile requested [0 1] global
%   n_sigma             number of sigmas
%
% OUTPUT:
%   smean   strong mean per column

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b6
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

    if nargin < 2 || ~isempty(thr_perc) 
        thr_perc = 0.99;
    end
    
    if nargin < 3 || ~isempty(thr_perc_global) 
        thr_perc = 0.97;
    end
    
    if nargin < 4 || ~isempty(n_sigma) 
        n_sigma = 5;
    end
    detrended_data = data - cumsum(repmat(strongMean(diff(data), 1, 0.97, 5), size(data, 1), 1));
end
