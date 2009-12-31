function [Q] = cofactor_matrix(elR, elM, snr_R, snr_M, sat, pivot)

% SYNTAX:
%   [Q] = cofactor_matrix(elR, elM, snr_R, snr_M, sat, pivot)
%
% INPUT:
%   elR = ROVER position (X,Y,Z)
%   elM = ROVER-SATELLITE code pseudorange
%   snr_R = MASTER position (X,Y,Z)
%   snr_M = MASTER-SATELLITE code pseudorange
%   sat = visible satellite configuration
%   pivot = pivot satellite
%
% OUTPUT:
%   Q = code-code or phase-phase co-factor matrix
%
% DESCRIPTION:
%   Co-factor matrix construction on the basis of selected weighting
%   strategy (determined by "weights" global variable).

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Media Center, Osaka City University, Japan
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

global weights snr_a snr_0 snr_1 snr_A

%total number of visible satellites
nsat = size(sat,1);
n = nsat - 1;

%PIVOT satellite index
i = find(pivot == sat);

if (weights == 0)
    
    %code-code or phase-phase co-factor matrix Q construction
    Q = 2*ones(n) + 2*eye(n);
    
else
    if (weights == 1)
        
        %weight vectors (elevation)
        q_R = 1 ./ (sin(elR * pi/180).^2);
        q_M = 1 ./ (sin(elM * pi/180).^2);
        
    elseif (weights == 2)
        
        %weight vectors (signal-to-noise ratio)
        q_R = 10.^(-(snr_R-snr_1)/snr_a) .* ((snr_A/10.^(-(snr_0-snr_1)/snr_a)-1)./(snr_0-snr_1).*(snr_R-snr_1)+1);
        q_R(snr_R >= snr_1) = 1;
        q_M = 10.^(-(snr_M-snr_1)/snr_a) .* ((snr_A/10.^(-(snr_0-snr_1)/snr_a)-1)./(snr_0-snr_1).*(snr_M-snr_1)+1);
        q_M(snr_M >= snr_1) = 1;
        
    elseif (weights == 3)
        %weight vectors (elevation and signal-to-noise ratio)
        q_R = 1 ./ (sin(elR * pi/180).^2) .* (10.^(-(snr_R-snr_1)/snr_a) .* ((snr_A/10.^(-(snr_0-snr_1)/snr_a)-1)./(snr_0-snr_1).*(snr_R-snr_1)+1));
        q_R(snr_R >= snr_1) = 1;
        q_M = 1 ./ (sin(elM * pi/180).^2) .* (10.^(-(snr_M-snr_1)/snr_a) .* ((snr_A/10.^(-(snr_0-snr_1)/snr_a)-1)./(snr_0-snr_1).*(snr_M-snr_1)+1));
        q_M(snr_M >= snr_1) = 1;
        
    end
    
    q_RP = q_R(i,1);             % ROVER-PIVOT
    q_MP = q_M(i,1);             % MASTER-PIVOT
    q_RS = q_R(sat~=pivot);      % ROVER-generic satellite (without pivot)
    q_MS = q_M(sat~=pivot);      % MASTER-generic satellite (without pivot)
    %q_RS(sat==pivot) = [];       % ROVER-generic satellite (without pivot)
    %q_MS(sat==pivot) = [];       % MASTER-generic satellite (without pivot)
    
    %code-code or phase-phase co-factor matrix Q construction
    Q = (q_RP + q_MP) * ones(n) + diag(q_RS + q_MS);
end