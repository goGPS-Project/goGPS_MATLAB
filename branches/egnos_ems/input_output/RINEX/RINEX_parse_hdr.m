function [Obs_types, pos_M, ifound_types, interval] = RINEX_parse_hdr(file)

% SYNTAX:
%   [Obs_types, pos_M, ifound_types, interval] = RINEX_parse_hdr(file);
%
% INPUT:
%   file = pointer to RINEX observation file
%
% OUTPUT:
%   Obs_types = string containing observation types (e.g. L1C1P1...)
%   pos_M = master station approximate position
%   ifound_types = boolean variable to check the correct acquisition of basic information
%   interval = Observation interval in seconds
%
% DESCRIPTION:
%   RINEX observation file header analysis.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.1 beta
%
% Copyright (C) Kai Borre
% Kai Borre 09-23-97
%
% Adapted by Mirko Reguzzoni, Eugenio Realini, 2009
% Portions of code contributed by Damiano Triglione, 2012
%----------------------------------------------------------------------------------------------

ifound_types = 0;
Obs_types = [];
pos_M = [];
interval = 1; %default to 1 second (1 Hz observations)

%parse first line
line = fgetl(file);

%check if the end of the header or the end of the file has been reached
while isempty(strfind(line,'END OF HEADER')) && ischar(line)
    %NOTE1: findstr is obsolete, so strfind is used
    %NOTE2: ischar is better than checking if line is the number -1.
    answer = strfind(line,'# / TYPES OF OBSERV');
    if ~isempty(answer)
        NoObs  = sscanf(line(1:6),'%d');
        NoObs2 = 0;
        if (NoObs > 9)
            NoObs2 = NoObs - 9;
            NoObs  = 9;
            line2  = fgetl(file);
        end
        for k = 1 : NoObs
           ot = sscanf(line(k*6+1:k*6+6),'%s');
           Obs_types = [Obs_types ot];
        end
        for k = 1 : NoObs2
           ot = sscanf(line2(k*6+1:k*6+6),'%s');
           Obs_types = [Obs_types ot];
        end
        
        ifound_types = 1;
    end
    answer = strfind(line,'APPROX POSITION XYZ');
    if ~isempty(answer)
        X = sscanf(line(1:14),'%f');
        Y = sscanf(line(15:28),'%f');
        Z = sscanf(line(29:42),'%f');
        pos_M = [X; Y; Z];
    end
    answer = strfind(line,'INTERVAL');
    if ~isempty(answer)
        interval = sscanf(line(1:10),'%f');
    end
    
    %parse next line
    line = fgetl(file);
end
