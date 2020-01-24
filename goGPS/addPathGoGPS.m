function addPathGoGPS()
% Script to add goGPS folders to path with black_list

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b6
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti, ...
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
    [l1, l2] = regexp(p,'(?<=:)[^:]*20180416T085957[^:]*:');

    for l = size(l1, 2) : -1 : 1
        p(l1(l) : l2(l)) = [];
    end

    % GACOS folder
    [l1, l2] = regexp(p,'(?<=:)[^:]*GACOS[\/|\\]example[^:]*:');

    for l = size(l1, 2) : -1 : 1
        p(l1(l) : l2(l)) = [];
    end

    % SINERGY folder
    [l1, l2] = regexp(p,'(?<=:)[^:]*Sinergy[\/|\\]maps[^:]*:');

    for l = size(l1, 2) : -1 : 1
        p(l1(l) : l2(l)) = [];
    end

    % GIT folders
    [l1, l2] = regexp(p,'(?<=:)[^:]*\.git[^:]*:');

    for l = size(l1, 2) : -1 : 1
        p(l1(l) : l2(l)) = [];
    end

    addpath(p);
end
