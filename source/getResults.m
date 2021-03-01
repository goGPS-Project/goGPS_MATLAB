% DESCRIPTION:
% Run this script to get the objects used by goGPS 
% when you finished your execution.
% Parsing these objects you can retrive all the results.
%
% OUTPUT:
%   core      the core processor object containing all the goGPS structures
%   rec       the last session array of Receivers
%   rec_list  when enabled in SESSION settings it stores all the Receivers
%             for all the processed sessions
%

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 Geomatics Research & Development srl (GReD)
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

core = Core.getCurrentCore();
rec = core.rec;

log = Core.getLogger();
log.addMarkedMessage('Now you should be able to see 2 new variables:');
log.addMessage(log.indent(' - core      the core processor object containing all the goGPS structures'));
log.addMessage(log.indent(' - rec       the array of Receivers'));
log.newLine();



