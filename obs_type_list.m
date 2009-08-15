function [Obs_types, ifound_types, eof] = obs_type_list(file)

% SYNTAX:
%   [Obs_types, ifound_types, eof] = obs_type_list(file);
%
% INPUT:
%   file = RINEX observation file
%
% OUTPUT:
%   Obs_types = string containing observation types (e.g. L1C1P1...)
%   ifound_types = boolean variable to check the correct acquisition of basic information
%   eof = boolean variable to check the end of file
%
% DESCRIPTION:
%   RINEX observation file header analysis.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 pre-alpha
%
% Copyright (C) Kai Borre 
% Kai Borre 09-23-97
%
% Adapted by Mirko Reguzzoni, Eugenio Realini, 2009
%
%----------------------------------------------------------------------------------------------

fid=fopen(file,'rt');
eof=0;
ifound_types=0;
Obs_types=[];
while ((ifound_types ==0) & (eof==0))			 
    line = fgetl(fid);
    answer = findstr(line,'END OF HEADER');
    if  ~isempty(answer)
        eof=1;
    end;
    if (line == -1)
        eof = 1; 
    end;
    answer = findstr(line,'# / TYPES OF OBSERV');
    if ~isempty(answer)
        [NObs, line] = strtok(line);
        NoObs = str2num(NObs);
        for k = 1:NoObs
           [ot, line] = strtok(line);
           Obs_types = [Obs_types ot];
        end;
        ifound_types = 1;
    end;
end;
