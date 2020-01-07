function [up_bound, lo_bound] = success_rate(D, L, bias)

% SYNTAX:
%   [up_bound, lo_bound] = success_rate(D, L, bias);
%
% INPUT:
%   Qahat = vcm of ambiguity
%   D = diagonal matrix of decorrelated Qahat
%   L = lower triangular matrix of decorrelated Qahat
%   bias = detected bias, if available
%
% OUTPUT:
%   up_bound = upper bound of success rate based on ADOP
%   lo_bound = lower bound of success rate based on bootstrapping method
%
% DESCRIPTION:
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
%  Written by:       Hendy F. Suhandri
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

%%define some input
%------------------
n    = length(D);
beta = L^-1 * bias;

%adop = det(Qahat)^(1/m);
s_ad = (D).^1/(2*n);
adop = prod(s_ad);

%%compute upper bound
%--------------------
u =1/(2*adop);
q = normcdf(u,0,1);
up_bound = (2*q - 1)^n;

%%compute lower bound
%--------------------

for i = 1:n
    u1 = (1 - 2*beta(i))/(2*sqrt(D(i)));
    u2 = (1 + 2*beta(i))/(2*sqrt(D(i)));
    q1 = normcdf(u1,0,1);
    q2 = normcdf(u2,0,1);
    lo_bound = prod(q1 + q2 - 1);
    %lo_bound = prod(lo_bound);
end
