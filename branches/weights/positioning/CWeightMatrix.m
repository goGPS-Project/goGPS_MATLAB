%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.2 beta
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

classdef CWeightMatrix < handle %#codegen
    %CWEIGHTMATRIX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant = true)
        %weight function parameters
        snr_a = 30;
        snr_0 = 10;
        snr_1 = 50;
        snr_A = 30;
        m = 3; %number of unknown parameters
        elea  = 10; % default value for the exponential elevation weight function
    end
    
    properties
        pivot_index; % pivot satellite index
        weights;
    end
    
    methods
        function obj = CWeightMatrix()
            obj.pivot_index = 0;
            obj.weights = 0;
        end
        
        % SYNTAX:
        %   [Q] = cofactor_matrix(elR, elM, snr_R, snr_M, pivot_index);
        %
        % OUTPUT:
        %   Q = code-code or phase-phase co-factor matrix
        %
        % DESCRIPTION:
        %   Co-factor matrix construction on the basis of selected weighting
        %   strategy (determined by "weights" global variable).
        
        function weightMatrix = getCofactorMatrixDD(obj, numSat, elR, elM, snr_R, snr_M)
            
            if (obj.weights == 0)
                %code-code or phase-phase co-factor matrix Q construction
                weightMatrix = 2*ones(numSat) + 2*eye(numSat);
            else
                
                q_R = getNonTrivialWeight( elR , snr_R );
                q_M = getNonTrivialWeight( elM , snr_M );
                
                q_RP = q_R( obj.pivot_index , 1 ); % ROVER-PIVOT
                q_MP = q_M( obj.pivot_index , 1 ); % MASTER-PIVOT
                q_R(obj.pivot_index) = [];
                q_M(obj.pivot_index) = [];
                q_RS = q_R;                % ROVER-generic satellite (without pivot)
                q_MS = q_M;                % MASTER-generic satellite (without pivot)
                
                %code-code or phase-phase co-factor matrix Q construction
                weightMatrix = (q_RP + q_MP) * ones(obj.numSat) + diag(q_RS + q_MS);
            end
        end
        
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
        function weightMatrix = cofactor_matrix_SA(obj, elR, snr_R)
            %total number of visible satellites
            n = length(elR);
            
            if (obj.weights == 0 || (~any(elR) || ~any(snr_R)))
                
                %code-code or phase-phase co-factor matrix Q construction
                weightMatrix = eye(n);
                
            else
                q_R = getNonTrivialWeight( elR , snr_R );
                
                %code-code or phase-phase co-factor matrix Q construction
                weightMatrix = diag(q_R);
            end
        end
    end
    
    methods (Access = private)
        function q = getNonTrivialWeight(obj, el, snr)
            if ( obj.weights == 1 )                
                %weight vectors (elevation)
                q = getSinElevationWeight(el);
                
            elseif ( obj.weights == 2 )
                
                %weight vectors (signal-to-noise ratio)
                q = getSnrRatioWeight( snr );
                
            elseif ( obj.weights == 3 )
                %weight vectors (elevation and signal-to-noise ratio)
                q = getSinElevationWeight( el ) .* getSnrRatioWeight( snr );
                
            elseif ( obj.weights == 4 )
                %weight vectors (elevation, exponential function)
                q = getExpElevationWeight( el );
            end            
        end
        
        %weight vectors (elevation, exponential function)
        function q = getExpElevationWeight( el )
            eleref = degtorad(min(el)); % this is the value for the elevation cut-off angle
            q = power( 1 + obj.elea*exp(-(degtorad( el ))/eleref) , 2 );
        end
        
        %weight vectors (elevation)
        function q = getSinElevationWeight( el )
            q = 1 ./ power(sin(degtorad( el )),2);
        end
        
        %weight vectors (signal-to-noise ratio)
        function q = getSnrRatioWeight( snr )
            q = 10.^( -( snr - obj.snr_1 ) / obj.snr_a ) .* ...
                (  (obj.snr_A / 10.^( -(obj.snr_0 - obj.snr_1 ) / obj.snr_a ) - 1 ) ./ ( obj.snr_0 - obj.snr_1 ) .* ( snr - obj.snr_1 ) + 1 );
            q(snr >= obj.snr_1) = 1;
        end
    end
    
end

