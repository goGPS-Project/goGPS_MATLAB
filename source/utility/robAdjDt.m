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
%  Contributors:     Andrea Gatti, Giulio Tagliaferro ...
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

function [dts] = robAdjDt(dph)
    % wrapper of fallback in case no mex
    % it is not the same code, instead of huber it uses the simple median
    
    % Fast simple approach: dt = median(sensor, 2, 'omitnan');
    
    [n_ep, n_stream] = size(dph);
    thrs = 0.02;
    
    dts = zeros(n_ep, 1);

    % for all the epoch
    for k = 1:n_ep
        % get the phases of the row
        dph_tmp = dph(k, :);
        id_ok = isfinite(dph_tmp);
        dph_tmp2 = dph_tmp(id_ok);
        if ~isempty(dph_tmp2) % if we have phases
            w = ones(size(dph_tmp2));
            j = 0;
            dt = 1e9;
            dt_prev = -1e9;
            while (j < 20 && abs(dt - dt_prev) > 0.005) % limit the reweight to 30 or less than 0.005 improvement
                dt_prev = dt;
                tmp = dph_tmp2 .* w;
                dt = sum(tmp) / sum(w); % weighted mean
                
                ares_n = abs(dph_tmp2 - dt) / thrs; % absolute residuals
                w = ones(size(dph_tmp2));
                idx_rw = find(ares_n > 1); % residual to be reweighted
                if ~isempty(idx_rw)
                    w(idx_rw) = 1 ./ ares_n(idx_rw).^2; % compute the weight
                end
                j = j + 1;  
            end
            dts(k) = dt; % put in the results
        end
    end
end


