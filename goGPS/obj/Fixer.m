%   CLASS Fixer
% =========================================================================
%
% DESCRIPTION
%   Class to manage fixing
%
% EXAMPLE
%   Net = Network();
%
% SEE ALSO
% FOR A LIST OF CONSTANTs and METHODS use doc Network

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 2
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2018 Mirko Reguzzoni, Eugenio Realini
%  Written by:                    Andrea Gatti
%  On the basis of the work of:   Andrea Nardo, 
%  Contributors:                  Eugenio Realini and Hendy F. Suhandri
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
classdef Fixer < handle
    properties (GetAccess = private, SetAccess = private)
        p0 = 0.001;
        mu = [];
    end
    
    methods (Access = private)
        function this = Fixer()
            % Creator
            % nothing necessary, this is mainly a Static class
        end
        
        function [amb_fixed, is_fixed, l_fixed] = fixAmbiguities(this, amb_float, C_amb_amb, approach)
            % Fix ambiguities using the selected approach 
            %
            % INPUT
            %   amb_float   array of float ambiguities
            %   C_amb_amb   covariance matrix of the ambiguities
            %
            %   approach    chose an approach amongo these:
            %                - 'lambda_ILS'     
            %                - 'lambda_bootstrapping'
            %                - 'lambda_partial'
            %                - 'bayesian'
            %                - 'best_integer_equivariant'
            %
            %
            % OUTPUTS:
            %
            %   amb_fixed: Array of size (n x ncands) with the estimated integer
            %              candidates, sorted according to the corresponding squared norms,
            %              best candidate first.
            %              For integer rounding and bootstrapping: ncands = 1
            %   is_fixed:  flag: true when fixing is computed successfully
            %   l_fixed:   Logical array of fixed ambiguities
            %
            % SYNTAX
            %   [amb_fixed, is_fixed, l_fixed] = this.fixAmbiguities(amb_float, C_amb_amb, approach)
            
            amb_fixed = amb_float;
            is_fixed = 0;
            l_fixed = false(size(amb_float));
            
            amb_ok = (abs(diag(C_amb_amb)) < 1); %& abs(fracFNI(amb_float)) < 0.3; % fix only valid ambiguities
            switch approach
                case {'lambda_ILS'}
                    try
                        [tmp_amb_fixed, sq_norm, success_rate] = LAMBDA(amb_float(amb_ok), full(10 * C_amb_amb(amb_ok, amb_ok)), 1, 'P0', this.p0, 'mu', this.mu);
                        %[tmp_amb_fixed,sqnorm,success_rate]=LAMBDA(amb_float(amb_ok), full(10 * C_amb_amb(amb_ok, amb_ok)),4,'P0',this.p0,'mu',mu);
                        
                        mu = ratioinv(this.p0, 1 - success_rate, length(tmp_amb_fixed));
                        ratio = sq_norm(1) / sq_norm(2);
                        
                        amb_fixed = repmat(amb_fixed, 1, size(tmp_amb_fixed, 2));
                        amb_fixed(amb_ok, :) = tmp_amb_fixed;
                        is_fixed  = ratio <= mu;
                        if is_fixed
                            is_fixed = is_fixed + ~all(amb_ok);
                        end
                        l_fixed   = abs(rem(amb_fixed,1)) < 1e-5;
                    catch
                        log = Logger.getInstance();
                        log.addError('LAMBDA fixing crashed, keeping float solution');
                        % Keep it float
                        amb_fixed = amb_float;
                        is_fixed = 0;
                        l_fixed = false(size(amb_float));
                    end
                case {'lambda_bootstrapping'}
                    try
                        [tmp_amb_fixed,sqnorm,success_rate]=LAMBDA(amb_float(amb_ok), full(10 * C_amb_amb(amb_ok, amb_ok)),4,'P0',this.p0,'mu',mu);
                        
                        mu = ratioinv(this.p0, 1 - success_rate, length(tmp_amb_fixed));
                        ratio = sq_norm(1) / sq_norm(2);
                        
                        amb_fixed = repmat(amb_fixed, 1, size(tmp_amb_fixed, 2));
                        amb_fixed(amb_ok, :) = tmp_amb_fixed;
                        is_fixed  = ratio <= mu;
                        if is_fixed
                            is_fixed = is_fixed + ~all(amb_ok);
                        end
                        l_fixed   = abs(rem(amb_fixed,1)) < 1e-5;
                    catch
                        log = Logger.getInstance();
                        log.addError('LAMBDA fixing crashed, keeping float solution');
                        % Keep it float
                        amb_fixed = amb_float;
                        is_fixed = 0;
                        l_fixed = false(size(amb_float));
                    end
                case {'lambda_partial'}
                    try
                        [tmp_amb_fixed, sq_norm, success_rate,~,~,nfx,mu] = LAMBDA(amb_float(amb_ok), full(10 * C_amb_amb(amb_ok, amb_ok)), 5, 'P0', 0.995, 'mu', this.mu);
                        
                        mu = ratioinv(this.p0, 1 - success_rate, length(tmp_amb_fixed));
                        ratio = sq_norm(1) / sq_norm(2);
                        
                        amb_fixed = repmat(amb_fixed, 1, size(tmp_amb_fixed, 2));
                        amb_fixed(amb_ok, :) = tmp_amb_fixed;
                        is_fixed  = ratio <= mu;
                        if is_fixed
                            is_fixed = is_fixed + ~all(amb_ok);
                        end
                        l_fixed   = abs(rem(amb_fixed,1)) < 1e-5;
                    catch
                        log = Logger.getInstance();
                        log.addError('LAMBDA fixing crashed, keeping float solution');
                        % Keep it float
                        amb_fixed = amb_float;
                        is_fixed = 0;
                        l_fixed = false(size(amb_float));
                    end
                case {'bayesian_with_monte_carlo'}
                    [tmp_amb_fixed] = this.bayesianAmbFixing(amb_float(amb_ok), full(10 * C_amb_amb(amb_ok, amb_ok)));
                    amb_fixed(amb_ok, :) = tmp_amb_fixed;
                    is_fixed = true;
                    l_fixed = amb_ok;
                case {'best_integer_equivariant'}
                    [tmp_amb_fixed] = this.BIE(amb_float(amb_ok), full(10 * C_amb_amb(amb_ok, amb_ok)));
                    amb_fixed(amb_ok, :) = tmp_amb_fixed;
                    is_fixed = true;
                    l_fixed = amb_ok;
                    
            end
        end
    end
    
    methods (Access = public, Static)
        
        function [amb_fixed, is_fixed, l_fixed] = fix(amb_float, C_amb_amb, approach)
            % Fix ambiguities using the selected approach 
            %
            % INPUT
            %   amb_float   array of float ambiguities
            %   C_amb_amb   covariance matrix of the ambiguities
            %
            %   approach    chose an approach amongo these:
            %                - 'lambda'     (yes) it's the only one available at the moment'
            %
            % OUTPUTS:
            %
            %   amb_fixed: Array of size (n x ncands) with the estimated integer
            %              candidates, sorted according to the corresponding squared norms,
            %              best candidate first.
            %              For integer rounding and bootstrapping: ncands = 1
            %   is_fixed:  flag: true when fixing is computed successfully
            %   l_fixed:   Logical array of fixed ambiguities
            %
            % SYNTAX
            %   [amb_fixed, is_fixed, l_fixed] = Fixer.fix(amb_float, C_amb_amb, approach)
            
            this = Fixer();
            if nargin < 3 || isempty(approach)
                approach = 'lambda';
            end
            [amb_fixed, is_fixed, l_fixed] = this.fixAmbiguities(amb_float, C_amb_amb, approach);
        end
        
        function [bie, sols, p_sols] = BIE(a,Q,n_cand)
            % compute the best integer equivariant 
            %
            % SYNTAX:
            %    [bie, sols, p_sols] = Fixer.BIE(a, Q, <n_cand>)
            if nargin < 3
                n_cand = 8;
            end
            % compute best integer equivariant using fir 8 terms
            if exist('decorrel') > 0
                [dQ,Z,L,D,da,iZt] = decorrel(Q,a);
            else
                dQ = Q;
                da = a;
                iZt = eye(size(Q));
            end
            [zfixed,sqnorm] = ssearch(da,L,D,n_cand);
            p_sols = exp(-1/2*sqnorm);
            p_sols = p_sols./sum(p_sols);
            bie = sum(zfixed.*repmat(p_sols,size(zfixed,1),1),2);
            bie = iZt*bie;
            if nargout > 2
                sols = iZt*zfixed;
            end
        end
        
        function [mass_pt] = bayesianAmbFixing(a,Q)
            % Same as best integer equivarint but  using a montecarlo simulation
            %
            % SYNTAX:
            %    [mass_pt] = bayesianAmbFixing(a,Q)
            if exist('decorrel') > 0
                [dQ,Z,L,D,da,iZt] = decorrel(Q,a);
            else
                dQ = Q;
                da = a;
                iZt = eye(size(Q));
            end
            R = chol(dQ);
            n_samples = 500000;
            %nc = det(2*pi*Q)^(0.5);
            %nc*exp(-1/2*u'*iQ*u);
            ra = repmat(da,1,n_samples)+ R*randn(numel(a),n_samples);
            ra(abs(fracFNI(ra)) > 0.1) = nan;
            ra = round(ra);
            n_valid_sml = sum(~isnan(ra),2);
            mass_pt = sum(nan2zero(ra),2)./n_valid_sml;
            %mass_pt = mode(ra,2);
            mass_pt = iZt*(mass_pt);
        end

        
    end
end
