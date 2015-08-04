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

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni,Eugenio Realini
% Portions of code contributed by Damiano Triglione (2012)
%
% Partially based on FOBS_TYP.M (EASY suite) by Kai Borre
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
                idL1 = {'L1C'};
                idL2 = {'L2W'};
                idC1 = {'C1C'};
                idP1 = {'C1P'};
                idP2 = {'C2W'};
                idS1 = {'S1C'};
                idS2 = {'S2W'};
                idD1 = {'D1C'};
                idD2 = {'D2W'};
            case 'R' %GLONASS
                idL1 = {'L1C'};
                idL2 = {'L2P'};
                idC1 = {'C1C'};
                idP1 = {'C1P'};
                idP2 = {'C2P'};
                idS1 = {'S1C'};
                idS2 = {'S2P'};
                idD1 = {'D1C'};
                idD2 = {'D2P'};
            case 'E' %Galileo
                idL1 = {'L1A';'L1B';'L1C';'L1X';'L1Z'};
                idL2 = {'L5X'};
                idC1 = {'C1A';'C1B';'C1C';'C1X';'C1Z'};
                idP1 = {'...'}; % <-- ?
                idP2 = {'C5X'};
                idS1 = {'S1A';'S1B';'S1C';'S1X';'S1Z'};
                idS2 = {'S5X'};
                idD1 = {'D1A';'D1B';'D1C';'D1X';'D1Z'};
                idD2 = {'D5X'};
            case 'C' %Compass/Beidou
                idL1 = {'L1I';'L2I'};
                idL2 = {'L7I'};
                idC1 = {'C1I';'C2I'};
                idP1 = {'...'}; % <-- ?
                idP2 = {'C7I'};
                idS1 = {'S1I';'S2I'};
                idS2 = {'S7I'};
                idD1 = {'D1I';'D2I'};
                idD2 = {'D7I'};
            case 'J' %QZSS
                idL1 = {'L1C'};
                idL2 = {'L2C'};
                idC1 = {'C1C'};
                idP1 = {'...'}; %QZSS does not use P1
                idP2 = {'C2C'};
                idS1 = {'S1C'};
                idS2 = {'S2C'};
                idD1 = {'D1C'};
                idD2 = {'D2C'};
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
