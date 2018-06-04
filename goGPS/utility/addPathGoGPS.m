% addPathGoGPS
% Script to add goGPS folders to path with black_list

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.6.0 alpha 2 - nightly
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2018 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Gatti Andrea
%  Contributors:     Gatti Andrea, ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

p = genpath(pwd);

% GACOS folder
lim = []; [lim(:,1), lim(:,2)] = regexp(p,'(?<=:)[^:]*20180416T085957[^:]*:');

for l = size(lim, 1) : -1 : 1
    p(lim(l, 1) : lim(l, 2)) = [];
end

% GACOS folder
lim = []; [lim(:,1), lim(:,2)] = regexp(p,'(?<=:)[^:]*GACOS[\/|\\]example[^:]*:');

for l = size(lim, 1) : -1 : 1
    p(lim(l, 1) : lim(l, 2)) = [];
end

% SINERGY folder
lim = []; [lim(:,1), lim(:,2)] = regexp(p,'(?<=:)[^:]*Sinergy[\/|\\]maps[^:]*:');

for l = size(lim, 1) : -1 : 1
    p(lim(l, 1) : lim(l, 2)) = [];
end

% GIT folders
lim = []; [lim(:,1), lim(:,2)] = regexp(p,'(?<=:)[^:]*git[^:]*:');

for l = size(lim, 1) : -1 : 1
    p(lim(l, 1) : lim(l, 2)) = [];
end

addpath(p);
