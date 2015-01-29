function [data] = decode_fastrax_it03(msg, constellations, wait_dlg)

% SYNTAX:
%   [data] = decode_fastrax_it03(msg, constellations, wait_dlg);
%
% INPUT:
%   msg = binary message received by the fastrax_it03 receiver
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%   wait_dlg = optional handler to waitbar figure
%
% OUTPUT:
%   data = cell-array that contains the decoded fastrax_it03 messages
%          (message class and id are in the first cell-array field)
%
% DESCRIPTION:
%   Fastrax_it03 messages decoding (also in sequence).

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Code contributed by Ivan Reguzzoni
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

if (nargin < 2 || isempty(constellations))
    [constellations] = goGNSS.initConstellation(1, 0, 0, 0, 0, 0);
end

%----------------------------------------------------------------------------------------------
% FTX MESSAGE HEADER
%----------------------------------------------------------------------------------------------

header1 = '3C';      % header (hexadecimal value) - 60 (decimal value) - ascii '<'
header2 = '21';      % header (hexadecimal value) - 33 (decimal value) - ascii '!'

codeHEX = [header1 header2];              % initial hexadecimal stream
codeBIN = dec2bin(hex2dec(codeHEX),16);   % initial binary stream

pos_FTX = strfind(msg, codeBIN);          % message initial index

%----------------------------------------------------------------------------------------------
% MESSAGE STARTING POINT
%----------------------------------------------------------------------------------------------

% output variable initialization
data = cell(0);
% find the index of the first message, if any
if (~isempty(pos_FTX))
    pos  = pos_FTX(1);
else
    return
end

%----------------------------------------------------------------------------------------------
% MESSAGE DECODING LOOP
%----------------------------------------------------------------------------------------------

% counter initialization
i = 0;

if (nargin == 3)
    waitbar(0,wait_dlg,'Decoding rover stream...')
end

while (pos + 15 <= length(msg))
    
    if (nargin == 3)
        waitbar(pos/length(msg),wait_dlg)
    end

    % check if there is an FTX header
    if (strcmp(msg(pos:pos+15),codeBIN))

        % counter increment
        i = i + 1;

        % skip the Fastrax_it03 header (16 bit)
        pos = pos + 16;

        if (pos + 39 <= length(msg))

            % ..:: Source ::..          (3 bit + 5 bits)
            % Source Node
            s_node = fbin2dec(msg(pos:pos+2));  pos = pos + 3;
            % Source Task
            s_task = fbin2dec(msg(pos:pos+4));  pos = pos + 5;

            % ..:: Destination ::..     (3 bit + 5 bits)
            % Destination Node
            d_node = fbin2dec(msg(pos:pos+2));  pos = pos + 3;
            % Destination Task
            d_task = fbin2dec(msg(pos:pos+4));  pos = pos + 5;

            % ..:: Msg Id ::.. (1 byte)
            msg_id = fbin2dec(msg(pos:pos+7));  pos = pos + 8;

            % ..:: Transaction Id ::..  (1 bit + 7 bits)
            % Response bit
            R = msg(pos);                         pos = pos + 1;
            % Transaction Msg
            trans_id = fbin2dec(msg(pos:pos+6));  pos = pos + 7;

            % ..:: Payload size ::..  (1 byte)
            LEN = fbin2dec(msg(pos:pos+7));  pos = pos + 8;

            %   pos + payload (0-255) + check + >
            if (pos +    (2*LEN)*8    + (8+8) + 8 - 1 <= length(msg))

                % checksum
                checksum = checksumFTX (msg(pos:pos+(2*LEN)*8-1), LEN);

                if strcmp(msg(pos+(2*LEN)*8:pos+(2*LEN)*8 + 15), checksum)

                    % message identification
                    switch msg_id

                        % TRACK (Track)                
                        case  3, 
                            [data(:,i)] = decode_FTX_TRACK(msg(pos:pos+((2*LEN)*8)-1), constellations);

                        % RAW (raw measurement) - PSEUDO
                        case  4, 
                            [data(:,i)] = decode_FTX_PSEUDO(msg(pos:pos+((2*LEN)*8)-1), constellations);

                        % EPH (ephemerides)
                        case 10, 
                            [data(:,i)] = decode_FTX_EPH(msg(pos:pos+((2*LEN)*8)-1), constellations);

                    end

                else
                    %fprintf('Checksum error!\n');
                end

                % skip the message body
                pos = pos + (2*LEN)*8;

                % skip the 2 checksum bytes
                pos = pos + 16;
                
                % skip > after check
                if strcmp(msg(pos:pos+7),dec2bin(62,8))
                    pos = pos + 8;
                else
                    %fprintf('Error: End not found!\n');
                end
                
            else
                break
            end

        else
            break
        end

    % check if there are other packages
    else

        % find the index of the first FTX message, if any
        pos_FTX = strfind(msg(pos:end),codeBIN);
        if ~isempty(pos_FTX)
            pos = pos + pos_FTX(1) - 1;
        else
            break
        end
    end

end

% // Acquisition related messages.
% #define ACQ_DATA_MSG_ID                        1
% #define PRN_STATUS_MSG_ID                      2
% 
% // Tracking related messages.
% #define TRACK_MSG_ID                           3
% #define PSEUDO_MSG_ID                          4
% #define TRACK_CONF_MSG_ID                      5
% 
% // AGC related messages.
% #define AGC_MSG_ID                             6
% 
% // Navigation related messages.
% #define NAV_FIX_MSG_ID                         7
% 
% // Time related messages.
% //#define GPS_TIME_MSG_ID                      8
% 
% // Frame decoding related messages.
% #define RAW_ALMANAC_MSG_ID                     9
% #define RAW_EPHEMERIS_MSG_ID                   10
% #define SV_HEALTH_MSG_ID                       11
% #define UTC_IONO_MODEL_MSG_ID                  12
% 
% // Predictions to search and OBS.
% #define PRN_PRED_MSG_ID                        13
% #define FREQ_PRED_MSG_ID                       14
% 
% // Bit decoding related messages.
% #define SUBFRAME_MSG_ID                        15
% #define RESERVED_MSG_17                        17   // Likely to be implemented - SUBFRAME_INFO_MSG_ID
% #define BIT_STREAM_MSG_ID                      18   
% 
% #define DBGTRACE_MSG_ID                        19   // Debugging trace message
% 
% // #define SENSOR_MEAS_MSG_ID                  20   // INS sensor measurements
% 
% #define ERROR_REPORT_MSG_ID                    21
% 
% #define TRACK_AID_MSG_ID                       22
% // Non-maskable messages.
% 
% // System messages start at 64 == 0x40.
% #define START_MSG_ID                           64
% #define STOP_MSG_ID                            65
% #define SLEEP_MSG_ID                           66
% #define STATUS_MSG_ID                          67
% #define ITALK_CONF_MSG_ID                      68
% #define SYSINFO_MSG_ID                         69
% #define ITALK_TASK_ROUTE_MSG_ID                70
% 
% // Parameter handling.
% #define PARAM_CTRL_MSG_ID                      71
% #define PARAMS_CHANGED_MSG_ID                  72
% 
% // Start & Stop completed messages
% #define START_COMPLETED_MSG_ID                 73
% #define STOP_COMPLETED_MSG_ID                  74
% 
% // Logging
% #define LOG_CMD_MSG_ID                         75
% 
% // System status
% #define SYSTEM_START_MSG_ID                    76
% 
% 
% // Non-maskable messages related to the GPS core.
% // Core messages start at 80 == 0x50.
% 
% // (acquisition)
% #define STOP_SEARCH_MSG_ID                     79
% #define SEARCH_MSG_ID                          80
% #define PRED_SEARCH_MSG_ID                     81
% #define SEARCH_DONE_MSG_ID                     82
% #define SE_DBG_MSG_ID                          84
% 
% // (tracking)
% #define DATA_INIT_MSG_ID                       85 
% #define TRACK_DEBUG_MSG_ID                     86
% #define TRACK_REINIT_MSG_ID                    87
% #define TRACK_DROP_MSG_ID                      88
% #define TRACK_FAST_MSG_ID                      89
% #define TRACK_STATUS_MSG_ID                    90
% 
% #define TRACK_DEBUG_DLL_MSG_ID                 91
% 
% // (search)
% #define HANDOVER_DATA_MSG_ID                   92
% #define CORE_SYNC_MSG_ID                       93 // 
% 
% // Navigation related non-maskable messages start at 96 == 0x60.
% #define WAAS_RAWDATA_MSG_ID                    96
% #define OP_DIFF_MSG_ID                         97
% #define ASSISTANCE_MSG_ID                      98
% #define PULL_FIX_MSG_ID                        99
% #define WAAS_FASTCORR_DBG_MSG_ID               100
% #define WAAS_SLOWCORR_MSG_ID                   101
% #define WAAS_IONO_DBG_MSG_ID                   102
% 
% // Debugging and memory access messages start at 112 = 0x70.
% // Memctrl messages.
% #define MEMCTRL_MSG_ID                          112
% 
% // Debug messages.
% #define DEBUG_CMD_MSG_ID                        113 // Debugging command message for common debug info querying
% 
% #define DEBUG_CPU_LOAD_MSG_ID                   118
% #define DEBUG_LOOPBACK_MSG_ID                   119
% #define DEBUG_SERVTERM_SWITCH_MSG_ID            120
% 
% // User messages start at 128 = 0x80
% 
% // The last message id is reserved for the STOP_TASK message.
% #define STOP_TASK_MSG_ID                        255

