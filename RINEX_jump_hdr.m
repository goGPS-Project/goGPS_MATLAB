function RINEX_jump_hdr(fid)

% SYNTAX:
%   RINEX_jump_hdr(fid);
%
% INPUT:
%   fid = pointer to the RINEX observation file
%
% DESCRIPTION:
%   Parse a RINEX observation file until the tag 'END OF HEADER' is found.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 alpha
%
% Copyright (C) 2009-2010 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Graduate School for Creative Cities, Osaka City University, Japan
%
% Partially based on ANHEADER.M (EASY suite) by Kai Borre
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

flag = 0;

while (flag == 0)
    lin = fgets(fid);
    answer = findstr(lin,'END OF HEADER');

    if ~isempty(answer)
        flag = 1;
    end
end
