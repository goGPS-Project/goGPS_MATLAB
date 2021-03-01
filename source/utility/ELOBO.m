function [x_k, s02_k, v_k, Cxx_k, N_inv] = ELOBO(A, Q, y0, b, N_inv, v, x, s02, subset)
% SYNTAX:
%   [x_k, s2_k, v_k, Cxx_k, N_inv] = ELOBO(A, Q, y, b, N_inv, v, x, s02, subset)
%
% INPUT:
%   A     : design matrix of the full problem
%   Q     : cofactor matrix of the full problem
%   y     : observation vector of the full problem
%   Ninv  : inverse of normal matrix of the full problem
%   x     : estimated parameters of the full problem
%   s02   : a posteriori sigma of the full problem
%   subset: index of observations of y that compose the testing block
%
% OUTPUT:
%   x_k  : estimated parameters without the testing block
%   s2_k : a posteriori sigma without the testing block
%   v_k  : estimated residuals without the testing block
%   Cxx_k: parameters covariance without the testing block
%
% DESCRIPTION:
%   perform ELOBO on blocks of correlated observations
%   identify one (block) outlier
%   reject it
%   re-estimate unknowns
%   according to the theory in "L. Biagi and S. Caldera. An efficient leave one block out approach to identify outliers.Journal of Applied Geodesy, Volume 7, Issue 1, pages 11..19, 2013"
%
%   this version is optimized to manage blocks!

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 Geomatics Research & Development srl (GReD)
%  Written by:       Stefano Caldera
%  Contributors:     Stefano Caldera, Andrea Gatti
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

    if isempty(N_inv)
        N_inv = cholinv(A'/Q*A);
    end
    [n_obs, n_col] = size(A);
    Ak = A(subset,:);
    Qk = Q(subset,subset);
    vk = v(subset,:);
    Bk = N_inv * Ak';
    Ck = Ak * Bk;
    Kk = Qk - Ck;
    if n_obs-n_col-length(subset) < 1
        log = Core.getLogger();
        log.addWarning(sprintf('Cluster cannot be checked. Redudancy = %d', n_obs - n_col - length(subset)));
    else
        BKk_inv = Bk / Kk;
        x_k = x - BKk_inv * vk;
        s02_k = ((n_obs - n_col) * s02 - vk' / Kk * vk) / (n_obs - n_col - length(subset));
        N_inv_fin = N_inv + BKk_inv * Bk';
        Cxx_k = s02_k * N_inv_fin;
        id_ok = setdiff(1 : size(A, 1), subset);
        y_k = A(id_ok, :) * x_k + b(id_ok);
        v_k = y0(id_ok) - y_k;
    end
end
