function [cur_line, obs_struct] = RINEX_get_obs(buf, cur_line, nSat, sat, sat_types, obs_col, nObsTypes, cc, active_sys, first_id_sys)

% SYNTAX:
%   [obs_struct] = RINEX_get_obs(file_RINEX, nSat, sat, sat_types, obs_col, nObsTypes, constellations);
%
% INPUT:
%   file_RINEX = observation RINEX file
%   nSat = number of available satellites (NOTE: RINEX v3.xx does not provide 'sat' and 'sat_types'))
%   sat = list of all available satellites
%   sat_types = ordered list of satellite types ('G' = GPS, 'R' = GLONASS, 'S' = SBAS) in the observation
%   obs_col = structure defining in which columns each observation type is to be found
%   nObsTypes = number of available observations
%   cc = Constellation_Collector object, contains the satus of the satellite systems in use
%
% OUTPUT:
%   obs_struct = struct with observations of enabled constellations
%
% DESCRIPTION:
%   Acquisition of RINEX observation data (code, phase and signal-to-noise ratio).

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
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

%total number of satellites (according to enabled constellations)
nSatTot = cc.getNumSat();

%array to contain the starting index for each constellation in the total array (with length = nSatTot)
sat_types_id = zeros(size(sat_types));

%output observations structure initialization
obs_struct = struct('L1', zeros(nSatTot,1), 'L2', zeros(nSatTot,1), ...
    'C1', zeros(nSatTot,1), ...
    'P1', zeros(nSatTot,1), 'P2', zeros(nSatTot,1), ...
    'S1', zeros(nSatTot,1), 'S2', zeros(nSatTot,1), ...
    'D1', zeros(nSatTot,1), 'D2', zeros(nSatTot,1));

obs_tmp = zeros(nSatTot,2);

if (~isempty(sat_types)) %RINEX v2.xx

    %convert constellations letter to starting index in the total array
    sat_types_id(sat_types == 'G') = first_id_sys(1);
    sat_types_id(sat_types == 'R') = first_id_sys(2);
    sat_types_id(sat_types == 'E') = first_id_sys(3);
    sat_types_id(sat_types == 'J') = first_id_sys(4);
    sat_types_id(sat_types == 'C') = first_id_sys(5);
    sat_types_id(sat_types == 'I') = first_id_sys(6);
    sat_types_id(sat_types == 'S') = first_id_sys(7);

    %observation types
    nLinesToRead = ceil(nObsTypes/5);  % I read a maximum of 5 obs per line => this is the number of lines to read
    nObsToRead = nLinesToRead * 5;     % Each line contains 5 observations

    %data read and assignment
    lin = char(32*uint8(ones(16*nObsToRead,1))'); % preallocate the line to read

    % Mask to filter all the possible observations (max 15)
    mask = false(16,nObsToRead);
    mask(2:14,:) = true;
    % preallocate a matrix of 15 strings (of length 14 characters)
    % notice that each observation element has a max length of 13 char,
    % the first character is added as a padding to separate the strings for
    % the command sscanf that now can be launched once for each satellite
    strObs = char(ones(14,nObsToRead)*32);

    for s = 1 : nSat

        %DEBUG
        if (sat(s) > 32)
            sat(s) = 32; %this happens only with SBAS; it's already fixed in the multi-constellation version
        end

        lin = char(lin*0+32); % clear line -> fill with spaces

        if (sat_types_id(s) ~= 0)
            % read all the lines containing the observations needed
            for l = 1 : (nLinesToRead)
                cur_line = cur_line + 1; linTmp = buf{cur_line};
                linLengthTmp = length(linTmp);
                lin((80*(l-1))+(1:linLengthTmp)) = linTmp;  %each line has a maximum length of 80 characters
            end
            linLength = 80*(nLinesToRead-1)+linLengthTmp;

            % convert the lines read from the RINEX file to a single matrix
            % containing all the observations
            strObs(1:13,:) = (reshape(lin(mask(:)),13,nObsToRead));
            fltObs = sscanf(strObs, '%f'); % read all the observations in the string
            obsId = 0; % index of the current observation
            % start parsing the observation string
            for k = find(sum(strObs == ' ') < 14)
                obsId = obsId+1;
                %obs = sscanf(lin(mask(:,k)), '%f');
                obs = fltObs(obsId);

                %check and assign the observation type
                if (any(~(k-obs_col.C1)))
                    obs_struct.C1(sat_types_id(s)+sat(s)-1) = obs;
                elseif (any(~(k-obs_col.L1)))
                    obs_struct.L1(sat_types_id(s)+sat(s)-1) = obs;
                    if (linLength>=16*k)
                        %convert signal-to-noise ratio
                        % faster conversion of a single ASCII character into an int
                        snr = mod((lin(16*k)-48),16);
                        obs_tmp(sat_types_id(s)+sat(s)-1, 1) = 6 * snr;
                    end
                elseif (any(~(k-obs_col.D1)))
                    obs_struct.D1(sat_types_id(s)+sat(s)-1) = obs;
                elseif (any(~(k-obs_col.S1)))
                    obs_struct.S1(sat_types_id(s)+sat(s)-1) = obs;
                elseif (any(~(k-obs_col.P1)))
                    obs_struct.P1(sat_types_id(s)+sat(s)-1) = obs;
                elseif (any(~(k-obs_col.L2)))
                    obs_struct.L2(sat_types_id(s)+sat(s)-1) = obs;
                    if (linLength>=16*k)
                        %convert signal-to-noise ratio
                        % faster conversion of a single ASCII character into an int
                        snr = mod((lin(16*k)-48),16);
                        obs_tmp(sat_types_id(s)+sat(s)-1, 2) = 6 * snr;
                    end
                elseif (any(~(k-obs_col.D2)))
                    obs_struct.D2(sat_types_id(s)+sat(s)-1) = obs;
                elseif (any(~(k-obs_col.S2)))
                    obs_struct.S2(sat_types_id(s)+sat(s)-1) = obs;
                elseif (any(~(k-obs_col.P2)))
                    obs_struct.P2(sat_types_id(s)+sat(s)-1) = obs;
                end
            end
            if (~obs_struct.S1(sat_types_id(s)+sat(s)-1))
                obs_struct.S1(sat_types_id(s)+sat(s)-1) = obs_tmp(sat_types_id(s)+sat(s)-1, 1);
            end
            if (~obs_struct.S2(sat_types_id(s)+sat(s)-1))
                obs_struct.S2(sat_types_id(s)+sat(s)-1) = obs_tmp(sat_types_id(s)+sat(s)-1, 2);
            end
        else
            %skip all the observation lines for the unused satellite
            cur_line = cur_line + nLinesToRead;
        end
    end

else %RINEX v3.xx

    for s = 1 : nSat

        %read the line for satellite 's'
        cur_line = cur_line + 1; linTmp = buf{cur_line};

        %read the constellation ID
        sysId = linTmp(1);

        %read the satellite PRN/slot number
        %satId = str2num(linTmp(2:3));
        % faster conversion of a couple of ASCII character into an 2 digits int
        satId = uint8(linTmp(2)-48)*10+linTmp(3)-48;

        %number of observations to be read on this line
        nObsToRead = nObsTypes.(sysId);

        %data read and assignment
        lin = char(32*uint8(ones(16*nObsToRead,1))'); % preallocate the line to read

        %keep only the part of 'lin' containing the observations
        linTmp = linTmp(4:end);
        linLengthTmp = length(linTmp);

        %fill in the 'lin' variable with the actual line read (may be shorter than expected)
        lin(1:linLengthTmp) = linTmp;
        linLength = length(lin);

        %compute the index in the total array and check if the current
        %constellation is required (if not, skip the line)
        switch (sysId)
            case 'G'
                if active_sys(1)
                    index = first_id_sys(1) + satId - 1;
                else
                    continue
                end
            case 'R'
                if active_sys(2)
                    index = first_id_sys(2) + satId - 1;
                else
                    continue
                end
            case 'E'
                if active_sys(3)
                    index = first_id_sys(3) + satId - 1;
                else
                    continue
                end
            case 'J'
                if active_sys(5)
                    index = first_id_sys(4) + satId - 1;
                else
                    continue
                end
            case 'C'
                if active_sys(4)
                    index = first_id_sys(5) + satId - 1;
                else
                    continue
                end
            case 'I'
                if active_sys(6)
                    index = first_id_sys(6) + satId - 1;
                else
                    continue
                end
            case 'S'
                if active_sys(6)
                    index = first_id_sys(7) + satId - 1;
                else
                    continue
                end
        end

        % Mask to filter all the possible observations
        mask = false(16,nObsToRead);
        mask(2:14,:) = true;
        % preallocate a matrix of n strings (of length 14 characters)
        % notice that each observation element has a max length of 13 char,
        % the first character is added as a padding to separate the strings for
        % the command sscanf that now can be launched once for each satellite
        strObs = char(ones(14,nObsToRead)*32);

        % convert the lines read from the RINEX file to a single matrix
        % containing all the observations
        strObs(1:13,:) = (reshape(lin(mask(:)),13,nObsToRead));
        fltObs = sscanf(strObs, '%f'); % read all the observations in the string
        obsId = 0; % index of the current observation
        % start parsing the observation string
        for k = find(sum(strObs == ' ') < 14)
                obsId = obsId+1;
                %obs = sscanf(lin(mask(:,k)), '%f');
                obs = fltObs(obsId);

                %check and assign the observation type
                if (any(~(k-obs_col.(sysId).C1)))
                    obs_struct.C1(index) = obs;
                elseif (any(~(k-obs_col.(sysId).L1)))
                    obs_struct.L1(index) = obs;
                    if (linLength>=16*k)
                        %convert signal-to-noise ratio
                        % faster conversion of a single ASCII character into an int
                        snr = mod((lin(16*k)-48),16);
                        obs_tmp(index, 1) = 6 * snr;
                    end
                elseif (any(~(k-obs_col.(sysId).D1)))
                    obs_struct.D1(index) = obs;
                elseif (any(~(k-obs_col.(sysId).S1)))
                    obs_struct.S1(index) = obs;
                elseif (any(~(k-obs_col.(sysId).P1)))
                    obs_struct.P1(index) = obs;
                elseif (any(~(k-obs_col.(sysId).L2)))
                    obs_struct.L2(index) = obs;
                    if (linLength>=16*k)
                        %convert signal-to-noise ratio
                        % faster conversion of a single ASCII character into an int
                        snr = mod((lin(16*k)-48),16);
                        obs_tmp(index, 2) = 6 * snr;
                    end
                elseif (any(~(k-obs_col.(sysId).P2)))
                    obs_struct.P2(index) = obs;
                elseif (any(~(k-obs_col.(sysId).D2)))
                    obs_struct.D2(index) = obs;
                elseif (any(~(k-obs_col.(sysId).S2)))
                    obs_struct.S2(index) = obs;
                end
        end
        if (~obs_struct.S1(index))
            obs_struct.S1(index) = obs_tmp(index, 1);
        end
        if (~obs_struct.S2(index))
            obs_struct.S2(index) = obs_tmp(index, 2);
        end
    end
end
