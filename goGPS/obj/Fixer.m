
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
%    |___/                    v 1.0 beta 4 ION
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:                    Andrea Gatti, Giulio Tagliaferro
%  On the basis of the work of:   Andrea Nardo (on the LAMBDA approach), 
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
            
            
            if ~strcmp(approach, 'sequential_best_integer_equivariant')
                % Regularize the ambiguity covariance matrix
                % => improves LAMBDA approach
                 amb_ok = (abs(diag(C_amb_amb)) < 1e3); %& abs(fracFNI(amb_float)) < 0.3; % fix only valid ambiguities
%                  C_amb_amb(amb_ok, amb_ok) = eigRegularizer(full(C_amb_amb(amb_ok, amb_ok)), 1e6);
%                C_amb_amb = (C_amb_amb + C_amb_amb') ./ 2; % Force it to be symmetric
            else
                amb_ok = (abs(diag(C_amb_amb)) < 1e10); %& abs(fracFNI(amb_float)) < 0.3; % fix only valid ambiguities
                C_amb_amb = (C_amb_amb + C_amb_amb') ./ 2; % Force it to be symmetric
                               
            end
            
            switch approach
                case {'lambda_ILS'}
                    [tmp_amb_fixed, sq_norm, success_rate] = LAMBDA(amb_float(amb_ok), full(C_amb_amb(amb_ok, amb_ok)), 1, 'P0', this.p0, 'mu', this.mu);
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
                case {'round_then_ILS'}
                    % round amboguities very nera to zero then apply
                    % lambbda
                    amb_ok = find(amb_ok);
                    sigma_float = sqrt(diag(C_amb_amb(amb_ok, amb_ok)));
                    pmf = Fixer.oneDimPMF(amb_float(amb_ok), sigma_float);
                    idx_round = abs(fracFNI(amb_float(amb_ok))) < 0.1 & pmf > 0.999; % we can round right away no need for LAMBDA
                    amb_fixed = repmat(amb_fixed, 1, 2);
                    l_fixed   = false(size(amb_fixed));
                    l_fixed(amb_ok(idx_round),:) = true;
                    amb_fixed(amb_ok(idx_round),:) = repmat(round(amb_float(amb_ok(idx_round))),1,2);
                    amb_float(amb_ok(~idx_round)) = amb_float(amb_ok(~idx_round)) - C_amb_amb(amb_ok(~idx_round), amb_ok(idx_round))*  inv(C_amb_amb(amb_ok(idx_round), amb_ok(idx_round))) * (amb_float(amb_ok(idx_round)) - amb_fixed(amb_ok(idx_round),1));
                    C_amb_amb(amb_ok(~idx_round), amb_ok(~idx_round))   =   C_amb_amb(amb_ok(~idx_round), amb_ok(~idx_round))   -   C_amb_amb(amb_ok(~idx_round), amb_ok(idx_round))*  inv(C_amb_amb(amb_ok(idx_round), amb_ok(idx_round))) * C_amb_amb(amb_ok(idx_round), amb_ok(~idx_round));
                    [amb_fixed(amb_ok(~idx_round),:), sq_norm, success_rate] = LAMBDA(amb_float(amb_ok(~idx_round)), full(C_amb_amb(amb_ok(~idx_round), amb_ok(~idx_round))), 1, 'P0', this.p0, 'mu', this.mu);
                    %[tmp_amb_fixed,sqnorm,success_rate]=LAMBDA(amb_float(amb_ok), full(10 * C_amb_amb(amb_ok, amb_ok)),4,'P0',this.p0,'mu',mu);
                    mu = ratioinv(this.p0, 1 - success_rate, length(amb_ok(~idx_round)));
                    ratio = sq_norm(1) / sq_norm(2);
                    is_fixed  = ratio <= mu;
                    if is_fixed
                        l_fixed(amb_ok(~idx_round),:) = true;

                    end
                    is_fixed = true;
                case {'lambda_bootstrapping'}
                    [tmp_amb_fixed,sqnorm,success_rate]=LAMBDA(amb_float(amb_ok), full(C_amb_amb(amb_ok, amb_ok)),4,'P0',this.p0,'mu',this.mu);
                    
                    mu = ratioinv(this.p0, 1 - success_rate, length(tmp_amb_fixed));
                    
                    amb_fixed = repmat(amb_fixed, 1, size(tmp_amb_fixed, 2));
                    amb_fixed(amb_ok, :) = tmp_amb_fixed;
                    is_fixed = true;
                    l_fixed   = amb_ok;
                case {'lambda_partial'}
                    [tmp_amb_fixed, sq_norm, success_rate,~,~,nfx,mu] = LAMBDA(amb_float(amb_ok), full(C_amb_amb(amb_ok, amb_ok)), 5, 'P0', 0.995, 'mu', this.mu);
                    is_fixed = true;
                    l_fixed   = amb_ok;
                    l_fixed(l_fixed) = abs(fracFNI(tmp_amb_fixed(:,1))) < 1e-9;
                    amb_fixed(amb_ok, 1) = tmp_amb_fixed(:,1);

                case {'bayesian_with_monte_carlo'}
                    [tmp_amb_fixed] = this.bayesianAmbFixing(amb_float(amb_ok), full( C_amb_amb(amb_ok, amb_ok)));
                    amb_fixed(amb_ok, :) = tmp_amb_fixed;
                    is_fixed = true;
                    l_fixed = amb_ok;
                case {'best_integer_equivariant'}
                    [tmp_amb_fixed] = this.BIE(amb_float(amb_ok), full(C_amb_amb(amb_ok, amb_ok)));
                    amb_fixed(amb_ok, :) = tmp_amb_fixed;
                    is_fixed = true;
                    l_fixed = amb_ok;
                case {'sequential_best_integer_equivariant'}
                    % boostrap solution starting from the most probable and
                    % not the one with lower formal errror
                    l_fixed = amb_ok;
                    [amb_fixed(amb_ok), l_fixed(amb_ok),  VCV_not_fixed] = this.mp_bootstrap(amb_float(amb_ok),(C_amb_amb(amb_ok, amb_ok)));
                    is_fixed = true;
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
                approach = 'lambda_ILS';
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
            if ~any(p_sols)
                p_sols(:) = 1;
                Core.getLogger.addWarning('Best integer equivariant found integer set probabilities too low')
            end
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
        
        function amb_fixed = bootstrap(amb_float,L)
            % compute the integer bootstrap estimator
            %
            % SYNTAX
            %  amb_fixed=bootstrap(amb_float,L)
            %
            % INPUT
            %  amb_float : float value of the desidered integer
            %  L : lower vtriangular decomposition (LDLt) of the Variance Covarince Matrix
            
            amb_fixed = zeros(size(amb_float));
            
            for i = 1:n-1
                amb_fixed(i) = amb_float(i);
                amb_float(i+1:end) = amb_float(i+1:end) - L(i+1:end,i)*(amb_float(i) - amb_fixed(i))
            end
            
            
            
            
        end
        function [amb_fixed ,l_fixed, VCV_not_fixed] = mp_bootstrap(amb_float,VCV)
            % compute the bootstrap solution not ordering ambiguites by
            % formal accuracy but by most probabl integer
            %
            % SYNTAX
            %   amb_fixed = mp_bootstrap(amb_float,L)
            n_a = length(amb_float);
            rounded = round(amb_float);
            sigma_float = sqrt(diag(VCV));
            pmf = Fixer.oneDimPMF(amb_float, sigma_float);
            [max_pmf,idx] = max(pmf);
            if max_pmf == 1 
                idx_one = find(pmf == 1);
                [~,idx_tmp] = min(abs(amb_float(idx_one) - round(amb_float(idx_one))));
                idx = idx_one(idx_tmp);
                
            end
            amb_fixed = amb_float;
            fixed_lid = false(n_a,1);
            i = 1;
            while   i <= n_a && pmf(idx) > 0.95
                    amb_fixed(idx) = round(amb_float(idx));
                    other_amb_lid = 1:n_a ~= idx;
                    amb_float(other_amb_lid) = amb_float(other_amb_lid) - VCV(other_amb_lid,idx)* 1/sigma_float(idx)^2 * (amb_float(idx) - amb_fixed(idx) ); % reduce the other amb for the fixed ones
                    VCV(other_amb_lid,other_amb_lid) = VCV(other_amb_lid,other_amb_lid) - VCV(other_amb_lid,idx) .* 1/sigma_float(idx)^2 * VCV(idx,other_amb_lid);  % reduce the vcv for the fixed ones
                    sigma_float = max(diag(VCV),1e-6);
                    sigma_float(fixed_lid) = 100; % this might happen beacuse we are reducing for ambiugities already fixed
                    sigma_float = sqrt(sigma_float);
                    fixed_lid(idx) = true;
                    nfixed_id = find(~fixed_lid);
                    pmf = Fixer.oneDimPMF(amb_float, sigma_float); % recompute the probability mass function
                    [max_pmf, mx_id] = max(pmf(nfixed_id)); % find the new most probable
                    if max_pmf == 1
                        idx_one = find(pmf(nfixed_id) == 1);
                        [~,idx_tmp] = min(abs(amb_float(nfixed_id(idx_one)) - round(amb_float(nfixed_id(idx_one)))));
                        mx_id = idx_one(idx_tmp);
                        
                    end
                    idx = nfixed_id(mx_id);
                    i = i +1;
            end
            VCV_not_fixed = VCV(~fixed_lid,~fixed_lid);
            l_fixed = fixed_lid;
        end
        
        function pmf = oneDimPMF(amb_float, sigma_float)
            % probability mass founction for the nearest integer
            % eq: 23.19 handbook of GNSS
            rounded = round(amb_float);
            pmf = normcdf((1 - 2*(amb_float -rounded))./(2*sigma_float)) + normcdf((1 + 2*(amb_float -rounded))./(2*sigma_float)) - 1;
        end
        
    end
end
