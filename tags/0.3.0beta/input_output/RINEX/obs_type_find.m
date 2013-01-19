function [col_L1, col_L2, col_C1, col_P1, col_P2, col_S1, col_S2, col_D1, col_D2] = obs_type_find(Obs_types)

% SYNTAX:
%   [col_L1, col_L2, col_C1, col_P1, col_P2, col_S1, col_S2, col_D1, col_D2] = obs_type_find(Obs_types);
%
% INPUT:
%   Obs_types = string containing observation types
%
% OUTPUT:
%   col_L1 = L1 column
%   col_L2 = L2 column
%   col_C1 = C1 column
%   col_P1 = P1 column
%   col_P2 = P2 column
%   col_S1 = S1 column
%   col_S2 = S2 column
%   col_D1 = D1 column
%   col_D2 = D2 column
%
% DESCRIPTION:
%   Selection of the column index for phase observations (L1, L2), for
%   code observations (C1, P1, P2), SNR ratios (S1, S2) and Doppler
%   measurements (D1, D2).

%----------------------------------------------------------------------------------------------
%                           goGPS v0.3.0 beta
%
% Copyright (C) 2009-2012 Mirko Reguzzoni,Eugenio Realini
% Portions of code contributed by Damiano Triglione (2012)
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
s = strfind(Obs_types, 'L1'); %findstr is obsolete, so strfind is used
col_L1 = (s+1)/2;

%search L2 column
s = strfind(Obs_types, 'L2');
col_L2 = (s+1)/2;

%search C1 column
s = strfind(Obs_types, 'C1');
col_C1 = (s+1)/2;

%search P1 column
s = strfind(Obs_types, 'P1');
col_P1 = (s+1)/2;

%search P2 column
s = strfind(Obs_types, 'P2');
col_P2 = (s+1)/2;

%search S1 column
s = strfind(Obs_types, 'S1');
col_S1 = (s+1)/2;

%search S2 column
s = strfind(Obs_types, 'S2');
col_S2 = (s+1)/2;

%search D1 column
s = strfind(Obs_types, 'D1');
col_D1 = (s+1)/2;

%search D2 column
s = strfind(Obs_types, 'D2');
col_D2 = (s+1)/2;