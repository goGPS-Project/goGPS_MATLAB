function [classOut, idOut] = ublox_UBX_codes(classIn, idIn)

% SYNTAX:
%   [classOut, idOut] = ublox_UBX_codes(classIn, idIn)
%
% INPUT:
%   classIn = u-blox message class (label - e.g. 'RXM')
%   idIn = u-blox message ID (label - e.g. 'RAW')
%
% OUTPUT:
%   classOut = u-blox message class (hex value - e.g. '02')
%   idOut = u-blox message ID (hex value - e.g. '10')
%
% DESCRIPTION:
%   Associate hex value/label pairs for u-blox message classes/IDs.

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

switch (classIn)

    case 'NAV'
        classOut = '01';
        switch (idIn)
            case 'POSECEF'  , idOut = '01';
            case 'POSLLH'   , idOut = '02';
            case 'POSUTM'   , idOut = '08';
            case 'DOP'      , idOut = '04';
            case 'STATUS'   , idOut = '03';
            case 'SOL'      , idOut = '06';
            case 'VELECEF'  , idOut = '11';
            case 'VELNED'   , idOut = '12';
            case 'TIMEGPS'  , idOut = '20';
            case 'TIMEUTC'  , idOut = '21';
            case 'CLOCK'    , idOut = '22';
            case 'SVINFO'   , idOut = '30';
            case 'DGPS'     , idOut = '31';
            case 'SBAS'     , idOut = '32';
            case 'EKFSTATUS', idOut = '40';
        end

    case 'RXM'
        classOut = '02';
        switch (idIn)
            case 'RAW'      , idOut = '10';
            case 'SVSI'     , idOut = '20';
            case 'SFRB'     , idOut = '11';
            case 'ALM'      , idOut = '30';
            case 'EPH'      , idOut = '31';
            case 'POSREQ'   , idOut = '40';
        end

    case 'INF'
        classOut = '04';
        switch (idIn)
            case 'ERROR'    , idOut = '00';
            case 'WARNING'  , idOut = '01';
            case 'NOTICE'   , idOut = '02';
            case 'TEST'     , idOut = '03';
            case 'DEBUG'    , idOut = '04';
            case 'USER'     , idOut = '07';
        end

    case 'ACK'
        classOut = '05';
        switch (idIn)
            case 'ACK'      , idOut = '01';
            case 'NAK'      , idOut = '00';
        end

    case 'CFG'
        classOut = '06';
        switch (idIn)
            case 'PRT'      , idOut = '00';
            case 'USB'      , idOut = '1B';
            case 'MSG'      , idOut = '01';
            case 'NMEA'     , idOut = '17';
            case 'RATE'     , idOut = '08';
            case 'CFG'      , idOut = '09';
            case 'TP'       , idOut = '07';
            case 'NAV2'     , idOut = '1A';
            case 'DAT'      , idOut = '06';
            case 'INF'      , idOut = '02';
            case 'RST'      , idOut = '04';
            case 'RXM'      , idOut = '11';
            case 'ANT'      , idOut = '13';
            case 'FXN'      , idOut = '0E';
            case 'SBAS'     , idOut = '16';
            case 'LIC'      , idOut = '80';
            case 'TM'       , idOut = '10';
            case 'TM2'      , idOut = '19';
            case 'TMODE'    , idOut = '1D';
            case 'EKF'      , idOut = '12';
        end

    case 'UPD'
        classOut = '09';
        switch (idIn)
            case 'DOWNL'    , idOut = '01';
            case 'UPLOAD'   , idOut = '02';
            case 'EXEC'     , idOut = '03';
            case 'MEMCPY'   , idOut = '04';
        end

    case 'MON'
        classOut = '0A';
        switch (idIn)
            case 'SCHD'     , idOut = '01';
            case 'IO'       , idOut = '02';
            case 'MSGPP'    , idOut = '06';
            case 'RXBUF'    , idOut = '07';
            case 'TXBUF'    , idOut = '08';
            case 'HW'       , idOut = '09';
            case 'IPC'      , idOut = '03';
            case 'USB'      , idOut = '0A';
            case 'VER'      , idOut = '04';
            case 'EXCEPT'   , idOut = '05';
        end

    case 'AID'
        classOut = '0B';
        switch (idIn)
            case 'REQ'      , idOut = '00';
            case 'DATA'     , idOut = '10';
            case 'INI'      , idOut = '01';
            case 'HUI'      , idOut = '02';
            case 'ALM'      , idOut = '30';
            case 'EPH'      , idOut = '31';
        end

    case 'TIM'
        classOut = '0D';
        switch (idIn)
            case 'TM'       , idOut = '02';
            case 'TM2'      , idOut = '03';
            case 'TP'       , idOut = '01';
            case 'SVIN'     , idOut = '04';
        end

    % following two are not really class labels, but they can be convenient
    case 'NMEA'
        classOut = 'F0';
        switch (idIn)
            case 'GGA'      , idOut = '00';
            case 'GLL'      , idOut = '01';
            case 'GSA'      , idOut = '02';
            case 'GSV'      , idOut = '03';
            case 'RMC'      , idOut = '04';
            case 'VTG'      , idOut = '05';
            case 'GRS'      , idOut = '06';
            case 'GST'      , idOut = '07';
            case 'ZDA'      , idOut = '08';
            case 'GBS'      , idOut = '09';
            case 'DTM'      , idOut = '0A';
        end

    case 'PUBX'
        classOut = 'F1';
        idOut = idIn;
end