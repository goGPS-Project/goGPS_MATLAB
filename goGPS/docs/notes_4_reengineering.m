%----------------------------------------------------------------------------------------------
%                           goGPS v0.5.9
% Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
% Written by:       Gatti Andrea
% Contributors:     Gatti Andrea, ...
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

% First step: reorganization of the parameters:
% Let's create a singleton class with the usefull parameters

% ambiguity status, it's in obs_rover
[] = LS_SA_CP(obs_rover, sat_info, proc_settings)

% Lagrange estimation
% suppose input x regularly sampled
% suppose odd polynomial
function y_pred = fastLI (first_x, step_x, y_obs, x_pred, p_deg)

    halfWinMin = ((p_deg - 1) / 2);
    halfWinMax = n_obs - ((p_deg + 1) / 2);
    
    n_pred = numel(x_pred);
    n_obs  = numel(y_obs);
    
    y_pred = zeros(size(x_pred));
    
    % find ids of the closer node to the estimation
    closer_node = round((x_pred - first_x) / step_x) + 1;

    % border checking - I should never predict close to the border
    i = 1;
    while (i < n_obs) && (closer_node(i) < halfWinMin)
        closer_node(i) = halfWinMin;
        i = i+1;
    end        
    i = n_obs;
    while (i > 0) && (closer_node(i) > halfWinMax)
        closer_node(i) = halfWinMin;
        i = i-1;
    end
    
    p_win    = (-halfWinMin : halfWinMin); % processing window relative ids
    p_around = [-halfWinMin : -1, 1 : halfWinMin]; % processing window relative ids without the center
    
    
    % for each element to be predicted
    for k = 1 : n_pred
        m = mean(y_obs(closer_node + p_win));
        
        L = 0;
        for j = p_win + 1 
            Lj = 1;
            for i = p_around + 1 
                Lj = Lj * (x_pred(k) - (first_x + (closer_node - 1) * step_x) ) / ((j-i) * step_x);
            end
            L = L + (y_obs(closer_node + j) - m) * Lj;
        end
        
        y_pred(k) = L + m;
    end
end
