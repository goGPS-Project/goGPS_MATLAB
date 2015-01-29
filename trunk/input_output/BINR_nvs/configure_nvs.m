function [serialObj] = configure_nvs(serialObj, COMportR, prot_par, rate)

% SYNTAX:
%   [serialObj] = configure_nvs(serialObj, COMportR, prot_par, rate);
%
% INPUT:
%   serialObj = handle to the rover serial object
%   COMportR = serial port the receiver is connected to
%   prot_par = receiver-specific parameters
%   rate = measurement rate to be set (default = 1 Hz)
%
% OUTPUT:
%   serialObj = handle to the rover serial object (it may have been re-created)
%
% DESCRIPTION:
%   Configure NVS receivers to be used with goGPS.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
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

if (nargin < 4)
    rate = 1;
end

% enable binary data output
fprintf('Enabling NVS receiver binary data output (NMEA output is automatically disabled)');

reply_BIN = nvs_message_format(serialObj, 'BIN', prot_par{2,1});
tries = 0;

while (~reply_BIN)
    tries = tries + 1;
    if (tries > 3)
        break
    end
    % close and delete old serial object
    try
        fclose(serialObj);
        delete(serialObj);
    catch
        stopasync(serialObj);
        fclose(serialObj);
        delete(serialObj);
    end
    % create new serial object
    serialObj = serial (COMportR,'BaudRate',prot_par{2,1});
    set(serialObj,'InputBufferSize',prot_par{3,1});
    fopen(serialObj);
    reply_BIN = nvs_message_format(serialObj, 'BIN', prot_par{2,1});
end

if (reply_BIN)
    fprintf(' done\n');
else
    fprintf(2, ' failed\n');
end

% set output rate (and raw measurement output)
fprintf('Enabling raw data output at %dHz measurement rate', rate);

reply_RATE = nvs_raw_output_rate(serialObj, rate);
tries = 0;

while (~reply_RATE)
    tries = tries + 1;
    if (tries > 3)
        break
    end
    % close and delete old serial object
    try
        fclose(serialObj);
        delete(serialObj);
    catch
        stopasync(serialObj);
        fclose(serialObj);
        delete(serialObj);
    end
    % create new serial object
    serialObj = serial (COMportR,'BaudRate',prot_par{2,1});
    set(serialObj,'InputBufferSize',prot_par{3,1});
    set(serialObj,'Parity','odd');
    fopen(serialObj);
    reply_RATE = nvs_raw_output_rate(serialObj, rate);
end

if (reply_RATE)
    fprintf(' done\n');
else
    fprintf(2, ' failed\n');
end
