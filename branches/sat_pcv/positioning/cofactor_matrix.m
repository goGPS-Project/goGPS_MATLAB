function [Q] = cofactor_matrix(elR, elM, snr_R, snr_M, pivot_index)

% SYNTAX:
%   [Q] = cofactor_matrix(elR, elM, snr_R, snr_M, pivot_index);
%
% INPUT:
%   elR = satellite elevations (ROVER)
%   elM = satellite elevations (MASTER)
%   snr_R = signal-to-noise ratio (ROVER)
%   snr_M = signal-to-noise ratio (MASTER)
%   pivot = pivot satellite index
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
%-------------------------------------------------------s---------------------------------------

global weights snr_a snr_0 snr_1 snr_A elea

%total number of visible satellites
n = length(elR) - 1;

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

    elseif (weights == 4)
        %weight vectors (elevation, exponential function)
        eleref = min(elR)* pi/180; % this is the value for the elevation cut-off angle
        q_R = (1 + elea*exp(-(elR * pi/180)/eleref)).^2;
        eleref = min(elM)* pi/180; % this is the value for the elevation cut-off angle
        q_M = (1 + elea*exp(-(elM * pi/180)/eleref)).^2;
    end

    q_RP = q_R(pivot_index,1); % ROVER-PIVOT
    q_MP = q_M(pivot_index,1); % MASTER-PIVOT
    q_R(pivot_index) = [];
    q_M(pivot_index) = [];
    q_RS = q_R;                % ROVER-generic satellite (without pivot)
    q_MS = q_M;                % MASTER-generic satellite (without pivot)

    %code-code or phase-phase co-factor matrix Q construction
    Q = (q_RP + q_MP) * ones(n) + diag(q_RS + q_MS);
end
