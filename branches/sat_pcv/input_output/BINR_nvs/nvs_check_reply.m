function [out] = nvs_check_reply(serialObj, check)

% SYNTAX:
%   [out] = nvs_check_reply(serialObj, check);
%
% INPUT:
%   serialObj = serial Object identifier
%   check  = column vector containing the values to look for in the reply (dec)
%
% OUTPUT:
%   out = acknowledge outcome
%
% DESCRIPTION:
%   Check acknowledge reply after issuing a command to NVS receivers.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%----------------------------------------------------------------------------------------------
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
%--------------------------------------------------------------------------
%--------------------

out = 0;

%time acquisition
start_time = toc;

%maximum waiting time
dtMax = 1;
dt = 0.1;
t = 1;

reply_1 = 0;
reply_2 = 0;

while (reply_1 ~= reply_2) || (reply_1 == 0)
    
    %time acquisition
    current_time = toc;
    
    %check if maximum waiting time is expired
    if (current_time - start_time > dtMax)
        return
    end
    
    % serial port checking
    reply_1 = get(serialObj, 'BytesAvailable');
    pause(dt);
    reply_2 = get(serialObj, 'BytesAvailable');
    
    if (mod(dtMax,t) == 0)
        fprintf('.')
    end
    
    t = t + dt;
end

reply = fread(serialObj, reply_1, 'uint8');

if (~isempty(strfind(reply',check')))
    out = 1;
end
