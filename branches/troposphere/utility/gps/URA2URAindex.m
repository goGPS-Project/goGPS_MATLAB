function [URAindex] = URA2URAindex(URA)

% SYNTAX:
%   [URAindex] = URA2URAindex(URA);
%
% INPUT:
%   URA = User Range Accuracy value in meters [vector]
%
% OUTPUT:
%   URAindex = User Range Accuracy index [vector]
%
% DESCRIPTION:
%   Converts URA values to URA indexes. From the Interface Specification 
%   document revision E (IS-GPS-200E), page 146.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
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

URAtable = [0.01 0.02 0.03 0.04 0.06 0.08 0.11 0.15 0.21 0.30 0.43 0.60 ...
    0.85 1.20 1.70 2.40 3.40 4.85 6.85 9.65 13.65 24.00 48.00 96.00 ...
    192.00 384.00 768.00 1536.00 3072.00 6144.00];

URAindextable = -15 : 1 : 15;

URAindex = zeros(size(URA));

for i = 1 : length(URA)
    if (isempty(URA(i)) || isnan(URA(i)))
        URAindex(i) = -16;
        continue
    end
    
    d = URAtable - URA(i);
    [~, j] = min(abs(d));
    if (d(j) >= 0)
        URAindex(i) = URAindextable(j);
    else
        URAindex(i) = URAindextable(j+1);
    end
end
