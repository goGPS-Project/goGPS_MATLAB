function undersamplingRINEX(filenameIN, filenameOUT, base, step, delta, wait_dlg)

% SYNTAX:
%   undersamplingRINEX(filenameIN, filenameOUT, base, step, delta, wait_dlg);
%
% INPUT:
%   filenameIN  = input RINEX observation file
%   filenameOUT = output RINEX observation file
%   base  = base timing (e.g. if first epoch should be 13  4  5 11 35  1, then base = 1) [sec]
%   step  = new sampling rate [sec]
%   delta = original sampling rate [sec]
%   wait_dlg = optional handler to waitbar figure (optional)
%
% OUTPUT:
%
%
% DESCRIPTION:
%   Undersamples a RINEX observation file.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b6
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
   if (~isempty(strfind(line,'INTERVAL')))
       line = sprintf('%10.3f                                                  INTERVAL            \n', step);
   end
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
        if (feof(fidIN) == 1) && ~isempty(line)
            fprintf(fidOUT,'%s',line);
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
