function [Q] = cofactor_matrix_SA(elR, snr_R)

% SYNTAX:
%   [Q] = cofactor_matrix_SA(elR, snr_R)
%
% INPUT:
%   elR = satellite elevations
%   snr_R = signal-to-noise ratio
%
% OUTPUT:
%   Q = code-code or phase-phase co-factor matrix
%
% DESCRIPTION:
%   Co-factor matrix construction on the basis of selected weighting
%   strategy (determined by "weights" global variable).

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
% Andrea Nardo, 22-Apr-2013: added exponential weighting function (weights == 4)
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

global weights snr_a snr_0 snr_1 snr_A elea

%total number of visible satellites
n = length(elR);

if (weights == 0 || (~any(elR) || ~any(snr_R)))
    
    %code-code or phase-phase co-factor matrix Q construction
    Q = eye(n);
    
else
    if (weights == 1)
        
        %weight vectors (elevation)
        q_R = 1 ./ (sin(elR * pi/180).^2);
        
    elseif (weights == 2)
        
        %weight vectors (signal-to-noise ratio)
        q_R = 10.^(-(snr_R-snr_1)/snr_a) .* ((snr_A/10.^(-(snr_0-snr_1)/snr_a)-1)./(snr_0-snr_1).*(snr_R-snr_1)+1);
        q_R(snr_R >= snr_1) = 1;
        
    elseif (weights == 3)
        %weight vectors (elevation and signal-to-noise ratio)
        q_R = 1 ./ (sin(elR * pi/180).^2) .* (10.^(-(snr_R-snr_1)/snr_a) .* ((snr_A/10.^(-(snr_0-snr_1)/snr_a)-1)./(snr_0-snr_1).*(snr_R-snr_1)+1));
        q_R(snr_R >= snr_1) = 1;
        
    elseif (weights == 4)
        %weight vectors (elevation, exponential function)
        eleref = min(elR)* pi/180; % this is the value for the elevation cut-off angle
        q_R = (1 + elea*exp(-(elR * pi/180)/eleref)).^2;
    end
    
    %code-code or phase-phase co-factor matrix Q construction
    Q = diag(q_R);
end
