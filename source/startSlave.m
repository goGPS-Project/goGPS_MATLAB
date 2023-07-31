
function startSlave(com_dir)
% SYNTAX:
%   startSlave(<com_dir>);
%
% INPUT:
%   com_dir       communication directory
%
% DESCRIPTION:
%   function launcher for goGPS slaves
%
% EXAMPLE:
%   startSlave('../com');
%
% COMPILATION STRING:
%    tic; mcc -v -d ../bin/ -m startSlave -a tai-utc.dat -a cls.csv -a icpl.csv -a nals.csv -a napl.csv -a remote_resource.ini -a credentials.txt -a app_settings.ini -a icons/*.png -a utility/thirdParty/guiLayoutToolbox/layout/+uix/Resources/*.png -R -singleCompThread; toc;
%

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0
%
%--------------------------------------------------------------------------
%  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
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

    if nargin == 1
        fprintf('Starting slave, my com dir is "%s"\n', com_dir);
        log = Logger.getInstance;
        log.setColorMode(0);
        %log.setVerbosityLev(0);
        gos = Go_Slave.getInstance(com_dir);
        gos.live;
        exit
    else
        fprintf('I need a parameter: com_dir\n');
    end
end
