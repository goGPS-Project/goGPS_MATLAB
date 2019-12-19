function [ntripstring] = NTRIP_string_generator(nmeastring)

% SYNTAX:
%   NTRIP_string_generator(nmeastring);
%
% INPUT:
%   nmeastring = NMEA sentence for Virtual Reference Station
%
% DESCRIPTION:
%   Builds the string used to connect to a positioning service by means of NTRIP.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
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

global master_ip ntrip_user ntrip_pw ntrip_mountpoint

%-------------------------------------------------------------------------------
% ENCODE USER AND PASSWORD (BASE 64 ENCODING)
%-------------------------------------------------------------------------------

encoded_auth = base64encode(sprintf('%s:%s',ntrip_user,ntrip_pw),'');

%-------------------------------------------------------------------------------
% NTRIP STRING GENERATION
%-------------------------------------------------------------------------------

if (nargin  == 0)
    nmeastring = [];
else
    nmeastring = sprintf('%s\r\n',nmeastring);
end

ntripstring = sprintf('GET /%s HTTP/1.1\r\nHost: %s\r\nUser-Agent: NTRIP goGPS\r\nConnection: close\r\nAuthorization: Basic %s\r\n\r\n%s\r\n', ntrip_mountpoint, master_ip, encoded_auth, nmeastring);
