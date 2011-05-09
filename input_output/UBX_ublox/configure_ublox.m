function [serialObj, reply_save] = configure_ublox(serialObj, COMportR, prot_par, rate)

% SYNTAX:
%   [serialObj, reply_save] = configure_ublox(serialObj, COMportR, prot_par, rate);
%
% INPUT:
%   serialObj = handle to the rover serial object
%   COMportR = serial port the receiver is connected to
%   prot_par = receiver-specific parameters
%   rate = measurement rate to be set (default = 1 Hz)
%
% OUTPUT:
%   serialObj = handle to the rover serial object (it may have been re-created)
%   reply_save = flag to verify that the receiver previous configuration was saved
%
% DESCRIPTION:
%   Configure u-blox receivers to be used with goGPS.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.2.0 beta
%
% Copyright (C) 2009-2011 Mirko Reguzzoni, Eugenio Realini
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
fprintf('Saving u-blox receiver configuration...\n');

reply_save = ublox_CFG_CFG(serialObj, 'save');
tries = 0;

while (~reply_save)
    tries = tries + 1;
    if (tries > 3)
        disp('It was not possible to save the receiver configuration.');
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
    reply_save = ublox_CFG_CFG(serialObj, 'save');
end

% set output rate
if (nargin < 4)
    rate = 1;
end
fprintf('Setting measurement rate to %dHz...\n', rate);

reply_RATE = ublox_CFG_RATE(serialObj, 1000/rate, 1, 1);
tries = 0;

while (~reply_RATE)
    tries = tries + 1;
    if (tries > 3)
        disp('It was not possible to set the receiver output rate to %dHz.\n', rate);
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

% enable raw measurements output
fprintf('Enabling u-blox receiver RAW measurements...\n');

reply_RAW = ublox_CFG_MSG(serialObj, 'RXM', 'RAW', 1);
tries = 0;

while (~reply_RAW)
    tries = tries + 1;
    if (tries > 3)
        disp('It was not possible to configure the receiver to provide RAW data.');
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

% disable subframe buffer output
fprintf('Disabling u-blox receiver subframe buffer (SFRB) messages...\n');

reply_SFRB = ublox_CFG_MSG(serialObj, 'RXM', 'SFRB', 0);
tries = 0;

while (~reply_SFRB)
    tries = tries + 1;
    if (tries > 3)
        disp('It was not possible to disable SFRB messages.');
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

