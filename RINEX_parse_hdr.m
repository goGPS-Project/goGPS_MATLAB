function [Obs_types, pos_M, ifound_types] = RINEX_parse_hdr(file)

% SYNTAX:
%   [Obs_types, pos_M, ifound_types] = RINEX_parse_hdr(file);
%
% INPUT:
%   file = pointer to RINEX observation file
%
% OUTPUT:
%   Obs_types = string containing observation types (e.g. L1C1P1...)
%   pos_M = master station approximate position
%   ifound_types = boolean variable to check the correct acquisition of basic information
%
% DESCRIPTION:
%   RINEX observation file header analysis.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.1 alpha
%
% Copyright (C) Kai Borre
% Kai Borre 09-23-97
%
% Adapted by Mirko Reguzzoni, Eugenio Realini, 2009
%
%----------------------------------------------------------------------------------------------

ifound_types = 0;
Obs_types = [];
pos_M = [];

%parse first line
line = fgetl(file);

%check if the end of the header or the end of the file has been reached
while isempty(findstr(line,'END OF HEADER')) & (line ~= -1)
    answer = findstr(line,'# / TYPES OF OBSERV');
    if ~isempty(answer)
        NoObs = sscanf(line(1:6),'%d');
        for k = 1 : NoObs
           ot = sscanf(line(k*6+1:k*6+6),'%s');
           Obs_types = [Obs_types ot];
        end
        ifound_types = 1;
    end
    answer = findstr(line,'APPROX POSITION XYZ');
    if ~isempty(answer)
        X = sscanf(line(1:14),'%f');
        Y = sscanf(line(15:28),'%f');
        Z = sscanf(line(29:42),'%f');
        pos_M = [X; Y; Z];
    end
    
    %parse next line
    line = fgetl(file);
end
