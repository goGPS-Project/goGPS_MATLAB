function [out] = remove_double_10h(msg, wait_dlg)

% SYNTAX:
%   [out] = remove_double_10h(data_in, wait_dlg);
%
% INPUT:
%   msg = BINR binary message that may contain double 10h values
%   wait_dlg = optional handler to waitbar figure
%
% OUTPUT:
%   out = BINR binary message with double 10h values compressed into one
%
% DESCRIPTION:
%   NVS BINR message data repeats 10h values twice, to distinguish them
%   from header/footer bytes. This function compresses them back to a
%   single value.

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
%----------------------------------------------------------------------------------------------

if (nargin < 2)
    %fprintf('Removing double 10h bytes from BINR data...\n');
else
    waitbar(0,wait_dlg,'Removing double 10h bytes from BINR data...')
end

%find all <DLE> bytes
pos = find(msg == 16);

%check the presence of <DLE><DLE>; in case remove the first
j = 1;
while (j <= length(pos))
    
    if (nargin == 2)
        waitbar(j/length(pos),wait_dlg)
    end
    
    if (length(msg) >= pos(j)+1)
        if (msg(pos(j)+1) ~= 16)
            %if it's a single <DLE>, leave it be
            pos(j) = [];
            j = j - 1;
        else
            %count the number of consecutive <DLE>s
            k = 1;
            while(msg(pos(j)+k) == 16)
                k = k + 1;
                if (length(msg) < pos(j)+k)
                    break
                end
            end
            rmv = j:2:j+k-1;
            len_rmv = length(rmv);
            pos(rmv) = [];
            j = j + len_rmv - 1;
        end
        j = j + 1;
    else
        pos(j) = [];
        break
    end
end
msg(pos) = [];
out = msg;

if (nargin == 2)
    waitbar(1,wait_dlg);
end
