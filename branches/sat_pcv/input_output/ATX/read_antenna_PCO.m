function [antPCO] = read_antenna_PCO(filename, antmod, date)

% SYNTAX:
%   [antPCO] = read_antenna_PCO(filename, antmod, date);
%
% INPUT:
%   filename = antenna phase center offset/variation file
%   antmod   = cell-array containing antenna model strings
%   date     = observation dates
%
% OUTPUT:
%   antPCO = antenna phase center offset [m]
%
% DESCRIPTION:
%   Extracts antenna phase center offset values from a PCO/PCV file in ATX format.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.2 beta
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

antPCO = zeros(3,1,length(antmod));

if (~isempty(filename))
    
    fid = fopen(filename,'r');
    
    if (fid ~= -1)
        
        found = 0;
        
        while (~feof(fid) && found < length(antmod))
            %parse the next line
            line = fgetl(fid);
            
            %check if any of the requested antenna models is in this line
            for m = 1 : length(antmod)
                validity_start = [];
                validity_end   = [];
                answer1 = strfind(line,[antmod{m} ' ']);
                answer2 = strfind(line,'TYPE / SERIAL NO');
                if (isempty(answer2))
                    break
                end
                if (~isempty(answer1) && ~isempty(answer2) && ~isempty(find(antmod{m} ~= ' ', 1)))

                    %get only the first offset, i.e. for L1
                    while (isempty(strfind(line,'NORTH / EAST / UP')))
                        line = fgetl(fid);
                        
                        %look for "VALID FROM" and "VALID UNTIL" lines (if satellite antenna)
                        if (strfind(line, 'VALID FROM'))
                            validity_start = [str2num(line(3:6)) str2num(line(11:12)) str2num(line(17:18)) str2num(line(23:24)) str2num(line(29:30)) str2num(line(34:43))];
                            line = fgetl(fid);
                            if (strfind(line, 'VALID UNTIL'))
                                validity_end = [str2num(line(3:6)) str2num(line(11:12)) str2num(line(17:18)) str2num(line(23:24)) str2num(line(29:30)) str2num(line(34:43))];
                            else
                                validity_end = Inf;
                            end
                        end
                    end
                    N = sscanf(line(1:10),'%f');
                    E = sscanf(line(11:20),'%f');
                    U = sscanf(line(21:30),'%f');

                    if (~isempty(validity_start)) %satellite antenna
                        if (datenum(date(1,:)) > datenum(validity_start) && datenum(date(end,:)) < datenum(validity_end))
                            antPCO(:,1,m) = [N; E; U]*1e-3;
                            
                            break
                        end

                    else  %receiver antenna
                        antPCO(:,1,m) = [E; N; U]*1e-3;
                        
                        found = found + 1;
                        
                        break
                    end
                end
            end
        end
        
        fclose(fid);
    end
else
    fprintf('Warning: PCO/PCV file not loaded.\n');
end
