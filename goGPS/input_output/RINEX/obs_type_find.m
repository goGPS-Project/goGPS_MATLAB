function [Obs_columns, nObs_types] = obs_type_find(Obs_types, sysId)

% SYNTAX:
%   [Obs_columns, nObs_types] = obs_type_find(Obs_types, sysId);
%
% INPUT:
%   Obs_types = cell of strings containing observation types
%               RINEX v2.xx --> e.g. L1C1P1...
%               RINEX v3.xx --> e.g. C1CL1CD1C...
%   sysId = cell-array containing one-letter identifiers for constellations
%
% OUTPUT:
%   Obs_columns = structure containing the column number of each observation type
%                 in the following fields:
%                   .L1 = L1 column
%                   .L2 = L2 column
%                   .C1 = C1 column
%                   .P1 = P1 column
%                   .P2 = P2 column
%                   .S1 = S1 column
%                   .S2 = S2 column
%                   .D1 = D1 column
%                   .D2 = D2 column
%                 In the case of RINEX v3.xx, an additional field is added
%                 for specifying the constellation, e.g.:
%                   .G.L1 (GPS)
%                   .R.L1 (GLONASS)
%   nObs_types = number of available observation types
%
% DESCRIPTION:
%   Detection of the column index for phase observations (L1, L2), for
%   code observations (C1, P1, P2), SNR ratios (S1, S2) and Doppler
%   measurements (D1, D2).

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 5 Merry Christmas
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Stefano Caldera
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

if (isempty(sysId)) %RINEX v2.xx

    nObs_types = size(Obs_types{1},2)/2;

    %search L1 column
    s1 = strfind(Obs_types{1}, 'L1');
    s2 = strfind(Obs_types{1}, 'LA');
    s = [s1 s2];
    col_L1 = (s+1)/2;

    %search L2 column
    s1 = strfind(Obs_types{1}, 'L2');
    s2 = strfind(Obs_types{1}, 'LC');
    s = [s1 s2];
    col_L2 = (s+1)/2;

    %search C1 column
    s1 = strfind(Obs_types{1}, 'C1');
    s2 = strfind(Obs_types{1}, 'CA');
    s = [s1 s2];
    col_C1 = (s+1)/2;

    %search P1 column
    s1 = strfind(Obs_types{1}, 'P1');
    s2 = strfind(Obs_types{1}, 'CA'); %QZSS does not use P1
    s = [s1 s2];
    col_P1 = (s+1)/2;

    %if RINEX v2.12 and GPS/GLONASS P1 observations are not available
    if (length(col_P1) ~= 2 && ~isempty(s2))
        %keep QZSS CA observations as C1
        col_P1 = [];
    end

    %search P2 column
    s1 = strfind(Obs_types{1}, 'P2');
    s2 = strfind(Obs_types{1}, 'CC');
    s = [s1 s2];
    col_P2 = (s+1)/2;

    %search S1 column
    s1 = strfind(Obs_types{1}, 'S1');
    s2 = strfind(Obs_types{1}, 'SA');
    s = [s1 s2];
    col_S1 = (s+1)/2;

    %search S2 column
    s1 = strfind(Obs_types{1}, 'S2');
    s2 = strfind(Obs_types{1}, 'SC');
    s = [s1 s2];
    col_S2 = (s+1)/2;

    %search D1 column
    s1 = strfind(Obs_types{1}, 'D1');
    s2 = strfind(Obs_types{1}, 'DA');
    s = [s1 s2];
    col_D1 = (s+1)/2;

    %search D2 column
    s1 = strfind(Obs_types{1}, 'D2');
    s2 = strfind(Obs_types{1}, 'DC');
    s = [s1 s2];
    col_D2 = (s+1)/2;

    Obs_columns.L1 = col_L1;
    Obs_columns.L2 = col_L2;
    Obs_columns.C1 = col_C1;
    Obs_columns.P1 = col_P1;
    Obs_columns.P2 = col_P2;
    Obs_columns.S1 = col_S1;
    Obs_columns.S2 = col_S2;
    Obs_columns.D1 = col_D1;
    Obs_columns.D2 = col_D2;

else %RINEX v3.xx
    for c = 1 : length(sysId)

        nObs_types.(sysId{c}) = size(Obs_types.(sysId{c}),2)/3;

        switch sysId{c}
            case 'G' %GPS
                idC1 = {'C1C'};               %L1
                idP1 = {'C1W';'C1P'};         %L1
                idL1 = {'L1W';'L1X';'L1C'};   %L1
                idS1 = {'S1W';'S1P';'S1C'};   %L1
                idD1 = {'D1W';'D1P';'D1C'};   %L1
                %--------------------------------
                %idC2 = {'C2W';'C2P';'C2C};    %L2,L2C (precedence to L2)
                idP2 = {'C2W';'C2P';'C2C';'C2S';'C2L';'C2X'};    %L2,L2C
                idL2 = {'L2W';'L2P';'L2C';'L2S';'L2L';'L2X'};    %L2,L2C
                idS2 = {'S2W';'S2P';'S2C';'S2S';'S2L';'S2X'};    %L2,L2C
                idD2 = {'D2W';'D2P';'D2C';'D2S';'D2L';'D2X'};    %L2,L2C
                %--------------------------------
                %idP2 = {'C5X'};               %L5
                %idL2 = {'L5X'};               %L5
                %idS2 = {'S5X'};               %L5
                %idD2 = {'D5X'};               %L5
                %--------------------------------
            case 'R' %GLONASS
                idC1 = {'C1C'};               %R1
                idP1 = {'C1P'};               %R1
                idL1 = {'L1C'};               %R1
                idS1 = {'S1C'};               %R1
                idD1 = {'D1C'};               %R1
                %--------------
                idP2 = {'C2P'};               %R2
                idL2 = {'L2P'};               %R2
                idS2 = {'S2P'};               %R2
                idD2 = {'D2P'};               %R2
            case 'E' %Galileo
                idC1 = {'C1X';'C1C'};         %E1
                idP1 = {'...'};               %E1
                idL1 = {'L1X';'L1C'};         %E1
                idS1 = {'S1X';'S1C'};         %E1
                idD1 = {'D1X';'D1C'};         %E1
                %--------------------------------
                idP2 = {'C5X';'C5Q'};        %E5a
                idL2 = {'L5X';'L5Q'};        %E5a
                idS2 = {'S5X';'S5Q'};        %E5a
                idD2 = {'D5X';'D5Q'};        %E5a
                %--------------------------------
                %idP2 = {'C7X';'C7Q'};        %E5b
                %idL2 = {'L7X';'L7Q'};        %E5b
                %idS2 = {'S7X';'S7Q'};        %E5b
                %idD2 = {'D7X';'D7Q'};        %E5b
                %--------------------------------
                %idP2 = {'C8X';'C8Q'};         %E5
                %idL2 = {'L8X';'L8Q'};         %E5
                %idS2 = {'S8X';'S8Q'};         %E5
                %idD2 = {'D8X';'D8Q'};         %E5
            case 'C' %Compass/Beidou
                idC1 = {'C1I';'C1Q';'C2I'};   %B1
                idP1 = {'...'};               %B1
                idL1 = {'L1I';'L1Q';'C2I'};   %B1
                idS1 = {'S1I';'S1Q';'S2I'};   %B1
                idD1 = {'D1I';'D1Q';'D2I'};   %B1
                %--------------------------------
                idL2 = {'L7I';'L7Q'};         %B2
                idP2 = {'C7I';'C7Q'};         %B2
                idS2 = {'S7I';'S7Q'};         %B2
                idD2 = {'D7I';'D7Q'};         %B2
            case 'J' %QZSS
                idC1 = {'C1X';'C1C'};         %J1
                idP1 = {'...'};               %J1
                idL1 = {'L1X';'L1C'};         %J1
                idS1 = {'S1X';'S1C'};         %J1
                idD1 = {'D1X';'D1C'};         %J1
                %--------------------------------
                idP2 = {'C2X';'C2C'};         %J2
                idL2 = {'L2X';'L2C'};         %J2
                idS2 = {'S2X';'S2C'};         %J2
                idD2 = {'D2X';'D2C'};         %J2
        end

        %search L1 column
        for i = 1 : length(idL1)
            s = strfind(Obs_types.(sysId{c}), idL1{i}); if (~isempty(s)), break, end;
        end
        col_L1 = (s+2)/3;

        %search L2 column
        for i = 1 : length(idL2)
            s = strfind(Obs_types.(sysId{c}), idL2{i}); if (~isempty(s)), break, end;
        end
        col_L2 = (s+2)/3;

        %search C1 column
        for i = 1 : length(idC1)
            s = strfind(Obs_types.(sysId{c}), idC1{i}); if (~isempty(s)), break, end;
        end
        col_C1 = (s+2)/3;

        %search P1 column
        for i = 1 : length(idP1)
            s = strfind(Obs_types.(sysId{c}), idP1{i}); if (~isempty(s)), break, end;
        end
        col_P1 = (s+2)/3;

        %search P2 column
        for i = 1 : length(idP2)
            s = strfind(Obs_types.(sysId{c}), idP2{i}); if (~isempty(s)), break, end;
        end
        col_P2 = (s+2)/3;

        %search S1 column
        for i = 1 : length(idS1)
            s = strfind(Obs_types.(sysId{c}), idS1{i}); if (~isempty(s)), break, end;
        end
        col_S1 = (s+2)/3;

        %search S2 column
        for i = 1 : length(idS2)
            s = strfind(Obs_types.(sysId{c}), idS2{i}); if (~isempty(s)), break, end;
        end
        col_S2 = (s+2)/3;

        %search D1 column
        for i = 1 : length(idD1)
            s = strfind(Obs_types.(sysId{c}), idD1{i}); if (~isempty(s)), break, end;
        end
        col_D1 = (s+2)/3;

        %search D2 column
        for i = 1 : length(idD2)
            s = strfind(Obs_types.(sysId{c}), idD2{i}); if (~isempty(s)), break, end;
        end
        col_D2 = (s+2)/3;

        Obs_columns.(sysId{c}).L1 = col_L1;
        Obs_columns.(sysId{c}).L2 = col_L2;
        Obs_columns.(sysId{c}).C1 = col_C1;
        Obs_columns.(sysId{c}).P1 = col_P1;
        Obs_columns.(sysId{c}).P2 = col_P2;
        Obs_columns.(sysId{c}).S1 = col_S1;
        Obs_columns.(sysId{c}).S2 = col_S2;
        Obs_columns.(sysId{c}).D1 = col_D1;
        Obs_columns.(sysId{c}).D2 = col_D2;
    end
end
