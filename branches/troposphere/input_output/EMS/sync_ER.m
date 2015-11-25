function [sync_data] = sync_ER(data_E, time_E, week_R, time_R)

% SYNTAX:
%   [sync_data] = sync_ER(data_E, time_E, week_R, time_R);
%
% INPUT:
%   data_E = EGNOS data to be synchronized with GPS data
%   time_E = time tags related to matrix data_E
%   week_R = GPS week for rover data
%   time_R = GPS seconds-of-week for rover data
%
% OUTPUT:
%   sync_data = matrix woth EGNOS data synchronized with GPS rover data
%
% DESCRIPTION:
%   Tool for synchronizing EGNOS data with the rover observations.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Code contributed by Giuliano Sironi, 2011
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
sync_data = zeros(length(time_R), length(data_E(1,:)));
 
for i = 1 : length(time_R)
    r = find(time_E(:,1)*604800 + time_E(:,2) <= week_R(i)*604800 + round(time_R(i)), 1, 'last');
    if ~isempty(r)
        sync_data(i,:) = data_E(r,:);
    end
end
