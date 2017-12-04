%   CLASS Least_Squares
% =========================================================================
%
% DESCRIPTION
%   Class to store element and general solver of LS adjustment
%
% EXAMPLE
%   ls = Least_Squares();
%


%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.2 beta 1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Giulio Tagliaferro
%  Contributors:      ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
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
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------
classdef Least_Squares < handle
    properties
        % Input
        A       % design matrix                          [ n_obs x n_par ]
        b       % array of known term                    [ n_obs x 1]
        y0      % observation array                      [ n_obs x 1]
        Q       % Cofactor matrix array                  [ n_obs x n_obs]
        
        % Output
        x       % estimated parameters                   [n_par x 1]
        res     % array of the residuals                 [n_obs x 1]
        s02     % (sigma0)^2                             [1 x 1]
        Cxx     % Covariance matrix of the parameters    [n_par x n_par]
        % Cyy   % Covariance matrix of the observations  [n_obs x n_obs] rarely estimated / not stored;
        
        % pre computed values
        iQ      % inverse of Cofactor matrix array       [ n_obs x n_obs]
        P       % A' / iQ                                [ n_par x n_obs]
        N       % Normal matrix                          [ n_par x n_par]
        iN      % inverse of the Normal Matrix           [ n_par x n_par]
       
    end
    
    properties (GetAccess = private)
         updated % flag to sync when the output parameters have been computed / the input change
        
        % [x res s02 Cxx iQ P/N iN]
    end
    
    methods
        function this = Least_Squares()
            % Constructor method
            % SYNTAX: this = Least_Squares()
            this.init()
        end
        
        function import(this, A, b, y0, Q, x, Cxx, s02, res)
            % Clear the object content
            % SYNTAX:
            %    this.import(A, b, y0, Q)
            %    this.import(A, b, y0, Q, x, Cxx, s02, res)
            this.A  = A;
            this.b  = b;
            this.y0 = y0;
            this.Q  = Q;
            
            if nargin == 8
                this.x = x;
                this.res = res;
                this.s02 = s02;
                this.Cxx = Cxx;
                this.updated = true(8,1);
            else
                this.updated = false(8,1);
            end
        end
        
        function [x, res, s02, Cxx] = solve(this, column)
            % Solve a LS problem
            % SYNTAX: [x, res, s02, Cxx, Cyy] = this.solve
            y0 = this.y0;
            b = this.b;
            if nargin > 1
                A = this.A(:,column);
            else
                A = this.A;
            end
            Q = this.Q;
            
            
            % least-squares solution
            if ~this.updated(6)
                this.P = A' / Q;
                this.N = this.P * A;
                this.updated(6) = true;
            end
            
            if nargout < 4
                if ~this.updated(1)
                    this.x = this.N \ (this.P * (y0 - b));
                    this.updated(1) = true;
                    x = this.x;
                else
                    x = this.x;
                end
            else
                % If I need to compute Cxx it is bbetter to pre compute the inverse of the normal matrix
                if ~this.updated(7)
                    try
                        this.iN = cholinv(full(this.N));
                    catch
                        this.iN = this.N^-1;
                    end
                    this.updated(7) = true;
                end
                if ~this.updated(1)
                    x = this.iN * this.P * (y0 - b);
                    this.x = x;
                    this.updated(1) = true;
                else
                    x = this.x;
                end
            end
            if nargout > 1
            % estimation of the variance of the observation error
            if ~this.updated(2)
                y_hat = A * x + b;
                res = y0 - y_hat;
                this.res = res;
                this.updated(2) = true;
            else
                res = this.res;
            end
            end
            
            if nargout > 2
                if ~this.updated(3)
                    [n, m] = size(A);
                    s02 = (res' * (Q \ res)) / (n - m);
                    this.s02 = s02;
                    this.updated(3) = true;
                else
                    s02 = this.s02;
                end
                if nargout > 3
                    if ~this.updated(4)
                        % covariance matrix
                        Cxx = s02 * this.iN;
                        this.Cxx = Cxx;
                        this.updated(4) = true;
                    else
                        Cxx = this.Cxx;
                    end
                end
            end
            
            function getCyy(this)
                % Get correlation matrix of the observations
                Cyy = s02 * A * iN * A';
            end
        end
        function clearUpdated(this)
            % DESCRIPTION: set upadted to flase
            this.updated = false(size(this.updated));
        end
        
    end
    
    methods (Access = 'private')
        function init(this)
            % Clear the object content
            % SYNTAX: this.init();
            
            this.A  = [];
            this.y0 = [];
            this.b  = [];
            this.Q  = [];
            
            this.iQ  = [];
            this.x   = [];
            this.Cxx = [];
            this.s02 = [];
            this.res = [];
            this.updated = false(8,1    );
        end
    end
    
    methods (Static)
        function [x, res, s02, Cxx, Cyy] = solver(y0, b, A, Q)
            % Solve a LS problem
            % SYNTAX: [x, res, s02, Cxx, Cyy] = solver(y0, b, A, Q)
            
            % least-squares solution
            P = A' / Q;
            N = P * A;
            
            if nargout < 4
                x = N \ P * (y0 - b);
            else
                % If I need to compute Cxx it is bbetter to pre compute the inverse of the normal matrix
                try
                    iN = cholinv(full(N));
                catch
                    iN = N^-1;
                end
                
                x = iN * P * (y0 - b);
            end
            % estimation of the variance of the observation error
            y_hat = A * x + b;
            res = y0 - y_hat;
            
            if nargout > 2
                [n, m] = size(A);
                s02 = (res' * (Q \ res)) / (n - m);
                if nargout > 3
                    % covariance matrix
                    Cxx = s02 * iN;
                    
                    if nargout > 4
                        A = A';
                        Cyy = s02 * A' * iN * A;
                    end
                end
            end
        end
        
        % Optimized versions of solver
        function [x] = solverOut1(y0, b, A, Q)
            % Solve a LS problem Optimized for one output
            % SYNTAX: [x] = solverOut1(y0, b, A, Q)
            
            % least-squares solution
            P = A' / Q;
            N = P * A;
            
            x = (P * A) \ P * (y0 - b);
        end
        
        function [x, res] = solverOut2(y0, b, A, Q)
            % Solve a LS problem Optimized for 2 outputs
            % SYNTAX: [x, res] = solverOut2(y0, b, A, Q)
            
            % least-squares solution
            P = A' / Q;
            N = P * A;
            
            x = N \ P * (y0 - b);
            
            % estimation of the variance of the observation error
            y_hat = A * x + b;
            res = y0 - y_hat;
        end
        
        function [x, res, s02] = solverOut3(y0, b, A, Q)
            % Solve a LS problem Optimized for 3 outputs
            % SYNTAX: [x, res, s02] = solverOut2(y0, b, A, Q)
            
            [n, m] = size(A);
            
            % least-squares solution
            P = A' / Q;
            N = P * A;
            
            x = N \ P * (y0 - b);
            
            % estimation of the variance of the observation error
            y_hat = A * x + b;
            res = y0 - y_hat;
            
            s02 = (res' /  Q * res) / (n - m);
        end
    end
end
