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
%                           goGPS v0.3.1 beta
%
% Copyright (C) 2009-2012 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Damiano Triglione (2012)
% Portions of code contributed by Andrea Gatti (2013)
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

    nSat = length(sat);

    % Convert constellations letter to index
    nConstellations = 3; % Number of constellations available 1: GPS, 2: GLONASS, 3: SBAS
    % nConstellations = 4; % Number of constellations available 1: GPS, 2: GLONASS, 3: SBAS 4: QZSS
    idGPS = 1;      sat_types_id = uint8(sat_types == 'G');
    idGLONASS = 2;  sat_types_id(sat_types == 'R') = idGLONASS;
    idSBAS = 3;     sat_types_id(sat_types == 'S') = idSBAS;
    %idQZSS = 4;     sat_types_id(sat_types == 'J') = idQZSS;
    
    %observation types
    [col_L1, col_L2, col_C1, col_P1, col_P2, col_S1, col_S2, col_D1, col_D2] = obs_type_find(obs_types);
    nObsTypes = size(obs_types,2)/2;
    nLinesToRead = ceil(nObsTypes/5);  % I read a maximum of 5 obs per line => this is the number of lines to read
    nObsToRead = nLinesToRead * 5;     % Each line contains 5 observations
    
    %observations structure initialization
    obs_struct = struct('L1', zeros(32,nConstellations), 'L2', zeros(32,nConstellations), ...
                        'C1', zeros(32,nConstellations), ...
                        'P1', zeros(32,nConstellations), 'P2', zeros(32,nConstellations), ...
                        'S1', zeros(32,nConstellations), 'S2', zeros(32,nConstellations), ...
                        'D1', zeros(32,nConstellations), 'D2', zeros(32,nConstellations));
                    
    obs_tmp = struct('TMP1', zeros(32,nConstellations), 'TMP2', zeros(32,nConstellations));
    
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
        lin = char(lin*0+32); % clear line -> fill with spaces
        
        % read all the lines containing the observations needed
        for l = 1 : (nLinesToRead)
            linTmp = fgetl(file_RINEX);
            linLengthTmp = length(linTmp);
            lin((80*(l-1))+(1:linLengthTmp)) = linTmp;  %each line has a maximum lenght of 80 characters
        end
        linLength = 80*(nLinesToRead-1)+linLengthTmp;
        
        % convert the lines read from the RINEX file to a single matrix
        % containing all the observations
        strObs(1:13,:) = (reshape(lin(mask(:)),13,nObsToRead)); 
        fltObs = sscanf(strObs, '%f'); % read all the numbers in the stringe
        obsId = 0; % index of the actual observation
        % start parsing the string of the observation
        for k = 1 : min(nObsTypes, floor(linLength/16))
            % check if the element is empty
            if(~strcmp(strObs(:,k)','              ')) % if the current val is missing (full of spaces)            
                obsId = obsId+1;
                %obs = sscanf(lin(mask(:,k)), '%f');
                obs = fltObs(obsId);
                
                %check and assign the observation type
                switch k
                    case col_L1
                        obs_struct.L1(sat(s),sat_types_id(s)) = obs;
                        if (linLength>=16*k)
                            %convert signal-to-noise ratio
                            % faster conversion of a single ASCII character into an int
                            snr = mod((lin(16*k)-48),16);
                            obs_tmp.TMP1(sat(s),sat_types_id(s)) = 6 * snr;
                        end
                    case col_L2
                        obs_struct.L2(sat(s),sat_types_id(s)) = obs;
                        if (linLength>=16*k)
                            %convert signal-to-noise ratio
                            % faster conversion of a single ASCII character into an int
                            snr = mod((lin(16*k)-48),16);
                            obs_tmp.TMP2(sat(s),sat_types_id(s)) = 6 * snr;
                        end
                    case col_C1
                        obs_struct.C1(sat(s),sat_types_id(s)) = obs;
                    case col_P1
                        obs_struct.P1(sat(s),sat_types_id(s)) = obs;
                    case col_P2
                        obs_struct.P2(sat(s),sat_types_id(s)) = obs;
                    case col_D1
                        obs_struct.D1(sat(s),sat_types_id(s)) = obs;
                    case col_D2
                        obs_struct.D2(sat(s),sat_types_id(s)) = obs;
                    case col_S1
                        obs_struct.S1(sat(s),sat_types_id(s)) = obs;
                    case col_S2
                        obs_struct.S2(sat(s),sat_types_id(s)) = obs;
                end
            end
        end
        if (~obs_struct.S1(sat(s)))
            obs_struct.S1(sat(s),sat_types_id(s)) = obs_tmp.TMP1(sat(s),sat_types_id(s));
        end
        if (~obs_struct.S2(sat(s)))
            obs_struct.S2(sat(s),sat_types_id(s)) = obs_tmp.TMP2(sat(s),sat_types_id(s));
        end
    end
    
    clear obs_tmp;
    
    % Distribute the observations in constellation separate structures
    obs_GPS = copyObsStruct(obs_struct, idGPS);
    obs_GLO = copyObsStruct(obs_struct, idGLONASS);
    obs_SBS = copyObsStruct(obs_struct, idSBAS);    
    % obs_QZSS = copyObsStruct(obs_struct, idQZSS); 
    
% -------------------------------------------------------------------------  
% End of function - start nested function declaration
% -------------------------------------------------------------------------    
    
    % Divide the multi constellation structure into single structures for each
    % constellation type
    function obsStruct = copyObsStruct(multiObsStruct, satTypeId)
        obsStruct.L1 = multiObsStruct.L1(:,satTypeId);
        obsStruct.L2 = multiObsStruct.L2(:,satTypeId);
        obsStruct.C1 = multiObsStruct.C1(:,satTypeId);
        obsStruct.P1 = multiObsStruct.P1(:,satTypeId);
        obsStruct.P2 = multiObsStruct.P2(:,satTypeId);
        obsStruct.D1 = multiObsStruct.D1(:,satTypeId);
        obsStruct.D2 = multiObsStruct.D2(:,satTypeId);
        obsStruct.S1 = multiObsStruct.S1(:,satTypeId);
        obsStruct.S2 = multiObsStruct.S2(:,satTypeId);
    end
end