function addPathGoGPS(flag_rem)
% Script to add goGPS folders to path with black_list

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 Geomatics Research & Development srl (GRed)
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

    if ~isdeployed
        p = genpath(pwd);

        % bin folder
        [l1, l2] = regexp(p,'(?<=:)[^:]*bin[\/|\\][^:]*:');

        for l = size(l1, 2) : -1 : 1
            p(l1(l) : l2(l)) = [];
        end

        [l1, l2] = regexp(p,'(?<=:)[^:]*bin[^:]*:');

        for l = size(l1, 2) : -1 : 1
            p(l1(l) : l2(l)) = [];
        end

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

        % SVN folders
        [l1, l2] = regexp(p,'(?<=:)[^:]*\.svn[^:]*:');

        % LOG folders
        [l1, l2] = regexp(p,'(?<=:)[^:]*\log[^:]*:');
        
        for l = size(l1, 2) : -1 : 1
            p(l1(l) : l2(l)) = [];
        end
        
        if nargin == 1 && flag_rem
            rmpath(p);
        else
            addpath(p);
        end
    end
end
