function undersamplingRINEX(filenameIN, filenameOUT, base, step, delta)

% SYNTAX:
%   undersamplingRINEX(filenameIN, filenameOUT, base, step, delta);
%
% INPUT:
%   filenameIN  = input RINEX observation file
%   filenameOUT = output RINEX observation file
%   base  = base timing (e.g. if first epoch should be 13  4  5 11 35  1, then base = 1) [sec]
%   step  = new sampling rate [sec]
%   delta = original sampling rate [sec]
%
% OUTPUT:
%
%
% DESCRIPTION:
%   Undersamples a RINEX observation file.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni,Eugenio Realini
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

%default values
if (nargin < 3)
    base  = 0;
end

if (nargin < 4)
    step  = 30;
end

if (nargin < 5)
    delta = 1;
end

% ------------------------------------------
% header copy and paste
% ------------------------------------------

fidIN  = fopen(filenameIN,'rt');
fidOUT = fopen(filenameOUT,'wt');

line = fgets(fidIN);
while (isempty(strfind(line,'END OF HEADER')))
   fprintf(fidOUT,'%s',line);
   line = fgets(fidIN);
end
fprintf(fidOUT,line);

% fclose(fidIN);
% fclose(fidOUT);

% ------------------------------------------
% undersampling
% ------------------------------------------

while (feof(fidIN) == 0)

    % time (seconds) extraction
	sec = roundmod(str2num(line(17:27)),delta);
	%fprintf('%d\n',sec);           % debugging

	if (mod(sec-base,step) == 0)    % if it is a header message
		fprintf(fidOUT,'%s',line);
		line = fgets(fidIN);
		%fprintf('%s',line);        % debugging
		while (feof(fidIN) == 0) & ...
		      ((line(1) ~= ' ') | (line(2) == ' ') | ...
		       (line(3) == ' ') | (line(4) ~= ' '))
			fprintf(fidOUT,'%s',line);
			line = fgets(fidIN);
			%fprintf('%s',line);    % debugging
		end
	else                            % if it is a data message
		line = fgets(fidIN);
		%fprintf('%s',line);        % debugging
		while (feof(fidIN) == 0) & ...
		      ((line(1) ~= ' ') | (line(2) == ' ') | ...
		       (line(3) == ' ') | (line(4) ~= ' '))
			line = fgets(fidIN);
			%fprintf('%s',line);    % debugging
		end
	end
end

fclose(fidIN);      % input file closure
fclose(fidOUT);     % output file closure 
