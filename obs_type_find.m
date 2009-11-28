function [col_L1, col_L2, col_P1_C1, col_P2] = obs_type_find(Obs_types)

% SYNTAX:
%   [col_L1, col_L2, col_P1_C1, col_P2] = obs_type_find(Obs_types);
%
% INPUT:
%   Obs_types = string containing observation types
%
% OUTPUT:
%   col_L1 = L1 column
%   col_L2 = L2 column
%   col_P1_C1 = P1/C1 column
%   col_P2 = P2 column
%
% DESCRIPTION:
%   Selection of the column index for phase observations (L1, L2) and for
%   code observations (C1/P1, P2).

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Media Center, Osaka City University, Japan
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

%search L1 column
s = findstr(Obs_types, 'L1');
col_L1 = (s+1)/2;

%search L2 column
s = findstr(Obs_types, 'L2');
col_L2 = (s+1)/2;

%search C1/P1 column
s = findstr(Obs_types, 'C1');
if (isempty(s))
    s = findstr(Obs_types, 'P1');
end;
col_P1_C1 = (s+1)/2;

%search P2 column
s = findstr(Obs_types, 'P2');
col_P2 = (s+1)/2;