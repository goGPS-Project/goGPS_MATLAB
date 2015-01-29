function [serialObj, reply_SAVE] = configure_ublox(serialObj, COMportR, prot_par, rate)

% SYNTAX:
%   [serialObj, reply_SAVE] = configure_ublox(serialObj, COMportR, prot_par, rate);
%
% INPUT:
%   serialObj = handle to the rover serial object
%   COMportR = serial port the receiver is connected to
%   prot_par = receiver-specific parameters
%   rate = measurement rate to be set (default = 1 Hz)
%
% OUTPUT:
%   serialObj = handle to the rover serial object (it may have been re-created)
%   reply_SAVE = flag to verify that the receiver previous configuration was saved
%
% DESCRIPTION:
%   Configure u-blox receivers to be used with goGPS.

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

% save receiver configuration
fprintf('Saving receiver configuration... ');

reply_SAVE = ublox_CFG_CFG(serialObj, 'save');
tries = 0;

while (~reply_SAVE)
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
    serialObj = serial(COMportR,'BaudRate',prot_par{2,1});
    set(serialObj,'InputBufferSize',prot_par{3,1});
    set(serialObj,'FlowControl','hardware');
    set(serialObj,'RequestToSend','on');
    fopen(serialObj);
    reply_SAVE = ublox_CFG_CFG(serialObj, 'save');
end

if (reply_SAVE)
    fprintf('done\n');
else
    fprintf(2, 'failed\n');
end

% set output rate
if (nargin < 4)
    rate = 1;
end
fprintf('Setting measurement rate to %dHz... ', rate);

reply_RATE = ublox_CFG_RATE(serialObj, 1000/rate, 1, 1);
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
    set(serialObj,'FlowControl','hardware');
    set(serialObj,'RequestToSend','on');
    fopen(serialObj);
    reply_RATE = ublox_CFG_RATE(serialObj, 1000/rate, 1, 1);
end

if (reply_RATE)
    fprintf('done\n');
else
    fprintf(2, 'failed\n');
end

% enable raw measurements output
fprintf('Enabling raw data output... ');

reply_RAW = ublox_CFG_MSG(serialObj, 'RXM', 'RAW', 1);
tries = 0;

while (~reply_RAW)
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
    set(serialObj,'FlowControl','hardware');
    set(serialObj,'RequestToSend','on');
    fopen(serialObj);
    reply_RAW = ublox_CFG_MSG(serialObj, 'RXM', 'RAW', 1);
end

if (reply_RAW)
    fprintf('done\n');
else
    fprintf(2, 'failed\n');
end

% disable subframe buffer output
fprintf('Disabling u-blox receiver subframe buffer (SFRB) messages... ');

reply_SFRB = ublox_CFG_MSG(serialObj, 'RXM', 'SFRB', 0);
tries = 0;

while (~reply_SFRB)
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
    set(serialObj,'FlowControl','hardware');
    set(serialObj,'RequestToSend','on');
    fopen(serialObj);
    reply_SFRB = ublox_CFG_MSG(serialObj, 'RXM', 'SFRB', 0);
end

if (reply_SFRB)
    fprintf('done\n');
else
    fprintf(2, 'failed\n');
end

% enable GGA messages, disable all other NMEA messages
fprintf('Configuring u-blox receiver NMEA messages:\n');

ublox_CFG_MSG(serialObj, 'NMEA', 'GGA', 1); fprintf('Enabling GGA...\n');
ublox_CFG_MSG(serialObj, 'NMEA', 'GLL', 0); fprintf('Disabling GLL ');
ublox_CFG_MSG(serialObj, 'NMEA', 'GSA', 0); fprintf('GSA ');
ublox_CFG_MSG(serialObj, 'NMEA', 'GSV', 0); fprintf('GSV ');
ublox_CFG_MSG(serialObj, 'NMEA', 'RMC', 0); fprintf('RMC ');
ublox_CFG_MSG(serialObj, 'NMEA', 'VTG', 0); fprintf('VTG ');
ublox_CFG_MSG(serialObj, 'NMEA', 'GRS', 0); fprintf('GRS ');
ublox_CFG_MSG(serialObj, 'NMEA', 'GST', 0); fprintf('GST ');
ublox_CFG_MSG(serialObj, 'NMEA', 'ZDA', 0); fprintf('ZDA ');
ublox_CFG_MSG(serialObj, 'NMEA', 'GBS', 0); fprintf('GBS ');
ublox_CFG_MSG(serialObj, 'NMEA', 'DTM', 0); fprintf('DTM ');
ublox_CFG_MSG(serialObj, 'PUBX', '00', 0); fprintf('PUBX00 ');
ublox_CFG_MSG(serialObj, 'PUBX', '01', 0); fprintf('PUBX01 ');
ublox_CFG_MSG(serialObj, 'PUBX', '03', 0); fprintf('PUBX03 ');
ublox_CFG_MSG(serialObj, 'PUBX', '04', 0); fprintf('PUBX04\n');

