function [antPCO] = read_antenna_PCO(filename, antmod)

% SYNTAX:
%   [antPCO] = read_antenna_PCO(filename, antmod);
%
% INPUT:
%   filename = antenna phase center offset/variation file
%   antmod = cell-array containing antenna model strings
%
% OUTPUT:
%   antPCO = antenna phase center offset [m]
%
% DESCRIPTION:
%   Extracts antenna phase center offset values from a PCO/PCV file in ATX format.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.2 beta
%
% Copyright (C) 2009-2013 Mirko Reguzzoni,Eugenio Realini
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
                answer1 = strfind(line,antmod{m});
                answer2 = strfind(line,'TYPE / SERIAL NO');
                if (~isempty(answer1) && ~isempty(answer2) && ~isempty(find(antmod{m} ~= ' ', 1)))
                    %get only the first offset, i.e. for L1
                    while (isempty(strfind(line,'NORTH / EAST / UP')))
                        line = fgetl(fid);
                    end
                    N = sscanf(line(1:10),'%f');
                    E = sscanf(line(11:20),'%f');
                    U = sscanf(line(21:30),'%f');
                    antPCO(:,1,m) = [E; N; U]*1e-3;
                    
                    found = found + 1;
                end
            end
        end
        
        fclose(fid);
    end
else
    fprintf('Warning: PCO/PCV file not loaded.\n');
end
