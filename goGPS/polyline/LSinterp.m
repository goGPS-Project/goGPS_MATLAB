function [m,q,s2m,s2q] = LSinterp(x0,y0,s2x,s2y,sxy)

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:
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


xm = mean(x0);
ym = mean(y0);

x0m = x0 - xm;
y0m = y0 - ym;

% approximate angular coefficient
if sum(x0m.^2) > 0
    m0 = sum(x0m.*y0m) / sum(x0m.^2);
else
    m0 = 0;
end

% weights
q = (s2y + m0^2 * s2x - 2*m0 * sxy).^(-1);

% normal matrix
N = [sum(q), sum(q.*x0m); sum(q.*x0m), sum(q.*(x0m.^2))];
invN = inv(N);

% known term
Y = [sum(q.*y0); sum(q.*y0.*x0m)];

% estimated parameters
par = invN * Y;

% output
m = par(2);
q = par(1) - m*xm;

s2m = invN(2,2);
s2q = invN(1,1);
