function [obs_GPS, obs_GLO, obs_SBS] = RINEX_get_obs(file_RINEX, sat, sat_types, obs_types)

% SYNTAX:
%   [obs_GPS, obs_GLO, obs_SBS] = RINEX_get_obs(file_RINEX, sat, sat_types, obs_types);
%
% INPUT:
%   file_RINEX = observation RINEX file
%   sat  = list of all visible satellites
%   sat_types = ordered list of satellite types ('G' = GPS, 'R' = GLONASS, 'S' = SBAS)
%   obs_types = observations types (e.g. C1L1P1...)
%
% OUTPUT:
%   obs_GPS = GPS observations
%   obs_GLO = GLONASS observations
%   obs_SBS = SBAS observations
%
% DESCRIPTION:
%   Acquisition of RINEX observation data (code, phase and signal-to-noise ratio).

%----------------------------------------------------------------------------------------------
%                           goGPS v0.2.0 beta
%
% Copyright (C) 2009-2011 Mirko Reguzzoni, Eugenio Realini
%
% Partially based on GRABDATA.M (EASY suite) by Kai Borre
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

num_sat = length(sat);

%observation types
[col_L1, col_L2, col_C1, col_P1, col_P2, col_S1, col_S2, col_D1, col_D2] = obs_type_find(obs_types);
num_obs_types = size(obs_types,2)/2;

%observations structure initialization
obs_GPS.L1 = zeros(32,1);
obs_GLO.L1 = zeros(32,1);
obs_SBS.L1 = zeros(32,1);

obs_GPS.L2 = zeros(32,1);
obs_GLO.L2 = zeros(32,1);
obs_SBS.L2 = zeros(32,1);

obs_GPS.C1 = zeros(32,1);
obs_GLO.C1 = zeros(32,1);
obs_SBS.C1 = zeros(32,1);

obs_GPS.P1 = zeros(32,1);
obs_GLO.P1 = zeros(32,1);
obs_SBS.P1 = zeros(32,1);

obs_GPS.P2 = zeros(32,1);
obs_GLO.P2 = zeros(32,1);
obs_SBS.P2 = zeros(32,1);

obs_GPS.S1 = zeros(32,1);
obs_GLO.S1 = zeros(32,1);
obs_SBS.S1 = zeros(32,1);

obs_GPS.S2 = zeros(32,1);
obs_GLO.S2 = zeros(32,1);
obs_SBS.S2 = zeros(32,1);

obs_GPS.D1 = zeros(32,1);
obs_GLO.D1 = zeros(32,1);
obs_SBS.D1 = zeros(32,1);

obs_GPS.D2 = zeros(32,1);
obs_GLO.D2 = zeros(32,1);
obs_SBS.D2 = zeros(32,1);

obs_GPS.TMP1 = zeros(32,1);
obs_GLO.TMP1 = zeros(32,1);
obs_SBS.TMP1 = zeros(32,1);

obs_GPS.TMP2 = zeros(32,1);
obs_GLO.TMP2 = zeros(32,1);
obs_SBS.TMP2 = zeros(32,1);

%data read and assignment
for s = 1 : num_sat
    lin = fgetl(file_RINEX);
    %more than 5 observation types --> 2 lines to be read
    if num_obs_types > 5
        lin_add = fgetl(file_RINEX);
        %add padding if necessary
        if (length(lin) < 80)
            for i = 1 : 80 - length(lin)
                lin = [lin ' '];
            end
        end
        lin = [lin lin_add];
    end
    %more than 10 observation types --> 3 lines to be read
    if num_obs_types > 10
        lin_add = fgetl(file_RINEX);
        %add padding if necessary
        if (length(lin) < 160)
            for i = 1 : 160 - length(lin)
                lin = [lin ' '];
            end
        end
        lin = [lin lin_add];
    end

    for k = 1 : num_obs_types

        %data save
        if (length(lin) < 2+16*k-2) | (isempty(sscanf(lin(2+16*(k-1):16*k-2),'%f')))
            obs = 0;
        else
            obs = sscanf(lin(2+16*(k-1):16*k-2),'%f');
        end
        
        %check and assign the observation type
        switch k
            case col_L1
                switch sat_types(s)
                    case 'G'
                        obs_GPS.L1(sat(s)) = obs;
                    case 'R'
                        obs_GLO.L1(sat(s)) = obs;
                    case 'S'
                        obs_SBS.L1(sat(s)) = obs;
                end
                if (numel(lin)>=16*k)
                    snr = sscanf(lin(16*k),'%f');
                    %convert signal-to-noise ratio
                    snr = snr * 6;
                    
                    if (~isempty(snr))
                        switch sat_types(s)
                            case 'G'
                                obs_GPS.TMP1(sat(s)) = snr;
                            case 'R'
                                obs_GLO.TMP1(sat(s)) = snr;
                            case 'S'
                                obs_SBS.TMP1(sat(s)) = snr;
                        end
                    end
                end
            case col_L2
                switch sat_types(s)
                    case 'G'
                        obs_GPS.L2(sat(s)) = obs;
                    case 'R'
                        obs_GLO.L2(sat(s)) = obs;
                    case 'S'
                        obs_SBS.L2(sat(s)) = obs;
                end
                if (numel(lin)>=16*k)
                    snr = sscanf(lin(16*k),'%f');
                    %convert signal-to-noise ratio
                    snr = snr * 6;
                    
                    if (~isempty(snr))
                        switch sat_types(s)
                            case 'G'
                                obs_GPS.TMP2(sat(s)) = snr;
                            case 'R'
                                obs_GLO.TMP2(sat(s)) = snr;
                            case 'S'
                                obs_SBS.TMP2(sat(s)) = snr;
                        end
                    end
                end
            case col_C1
                switch sat_types(s)
                    case 'G'
                        obs_GPS.C1(sat(s)) = obs;
                    case 'R'
                        obs_GLO.C1(sat(s)) = obs;
                    case 'S'
                        obs_SBS.C1(sat(s)) = obs;
                end
            case col_P1
                switch sat_types(s)
                    case 'G'
                        obs_GPS.P1(sat(s)) = obs;
                    case 'R'
                        obs_GLO.P1(sat(s)) = obs;
                    case 'S'
                        obs_SBS.P1(sat(s)) = obs;
                end
            case col_P2
                switch sat_types(s)
                    case 'G'
                        obs_GPS.P2(sat(s)) = obs;
                    case 'R'
                        obs_GLO.P2(sat(s)) = obs;
                    case 'S'
                        obs_SBS.P2(sat(s)) = obs;
                end
            case col_S1
                switch sat_types(s)
                    case 'G'
                        obs_GPS.S1(sat(s)) = obs;
                    case 'R'
                        obs_GLO.S1(sat(s)) = obs;
                    case 'S'
                        obs_SBS.S1(sat(s)) = obs;
                end
            case col_S2
                switch sat_types(s)
                    case 'G'
                        obs_GPS.S2(sat(s)) = obs;
                    case 'R'
                        obs_GLO.S2(sat(s)) = obs;
                    case 'S'
                        obs_SBS.S2(sat(s)) = obs;
                end
            case col_D1
                switch sat_types(s)
                    case 'G'
                        obs_GPS.D1(sat(s)) = obs;
                    case 'R'
                        obs_GLO.D1(sat(s)) = obs;
                    case 'S'
                        obs_SBS.D1(sat(s)) = obs;
                end
            case col_D2
                switch sat_types(s)
                    case 'G'
                        obs_GPS.D2(sat(s)) = obs;
                    case 'R'
                        obs_GLO.D2(sat(s)) = obs;
                    case 'S'
                        obs_SBS.D2(sat(s)) = obs;
                end
        end
    end
    switch sat_types(s)
        case 'G'
            if (~obs_GPS.S1(sat(s)))
                obs_GPS.S1(sat(s)) = obs_GPS.TMP1(sat(s));
            end
            if (~obs_GPS.S2(sat(s)))
                obs_GPS.S2(sat(s)) = obs_GPS.TMP2(sat(s));
            end
        case 'R'
            if (~obs_GLO.S1(sat(s)))
                obs_GLO.S1(sat(s)) = obs_GLO.TMP1(sat(s));
            end
            if (~obs_GLO.S2(sat(s)))
                obs_GLO.S2(sat(s)) = obs_GLO.TMP2(sat(s));
            end
        case 'S'
            if (~obs_SBS.S1(sat(s)))
                obs_SBS.S1(sat(s)) = obs_SBS.TMP1(sat(s));
            end
            if (~obs_SBS.S2(sat(s)))
                obs_SBS.S2(sat(s)) = obs_SBS.TMP2(sat(s));
            end
    end
end

%remove unnecessary fields
obs_GPS = rmfield(obs_GPS,'TMP1');
obs_GPS = rmfield(obs_GPS,'TMP2');
obs_GLO = rmfield(obs_GLO,'TMP1');
obs_GLO = rmfield(obs_GLO,'TMP2');
obs_SBS = rmfield(obs_SBS,'TMP1');
obs_SBS = rmfield(obs_SBS,'TMP2');