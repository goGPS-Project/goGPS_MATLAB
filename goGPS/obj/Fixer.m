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
        iar_method = 1;
        p0 = 0.001;
        mu = [];
    end
    
    methods (Access = private)
        function this = Fixer()
            % Creator
            % nothing necessary, this is mainly a Static class
        end
        
        function [amb_fixed, is_fixed, l_fixed] = fixAmbiguities(this, amb_float, C_amb_amb, approach, iar_method)
            % Fix ambiguities using the selected approach 
            %
            % INPUT
            %   amb_float   array of float ambiguities
            %   C_amb_amb   covariance matrix of the ambiguities
            %
            %   approach    chose an approach amongo these:
            %                - 'lambda'     (yes) it's the only one available at the moment'
            %
            %   iar_method  Integer Ambiguities Resolution method
            %               Lambda approach:
            %                   1: ILS method based on search-and-shrink [DEFAULT]
            %                   2: ILS method based enumeration in search
            %                   3: integer rounding method
            %                   4: integer bootstrapping method
            %                   5: PAR with the input P0 of user-defined success rate
            %                   6: ILS method with Ratio Test (uses search-and shrink)
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
            %   [amb_fixed, is_fixed, l_fixed] = this.fixAmbiguities(amb_float, C_amb_amb, approach, iar_method)
            
            amb_fixed = amb_float;
            is_fixed = 0;
            l_fixed = false(size(amb_float));
            
            amb_ok = (abs(diag(C_amb_amb)) < 1); %& abs(fracFNI(amb_float)) < 0.3; % fix only valid ambiguities
            switch approach
                case {'lambda'}
                    try
                    [tmp_amb_fixed, sq_norm, success_rate] = LAMBDA(amb_float(amb_ok), full(10 * C_amb_amb(amb_ok, amb_ok)), iar_method, 'P0', this.p0, 'mu', this.mu);
                    %[tmp_amb_fixed, sq_norm, success_rate,~,~,nfx,mu] = LAMBDA(amb_float(amb_ok), full(10 * C_amb_amb(amb_ok, amb_ok)), 5, 'P0', 0.995, 'mu', this.mu);
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
                case {'bayesian'}
                    [tmp_amb_fixed] = bayesianAmbFixing(amb_float(amb_ok), full(10 * C_amb_amb(amb_ok, amb_ok)));
                    amb_fixed(amb_ok, :) = tmp_amb_fixed;;
                    is_fixed = true;
                    l_fixed = amb_ok;
                    
            end
        end
    end
    
    methods (Access = public, Static)
        
        function [amb_fixed, is_fixed, l_fixed] = fix(amb_float, C_amb_amb, approach, iar_method)
            % Fix ambiguities using the selected approach 
            %
            % INPUT
            %   amb_float   array of float ambiguities
            %   C_amb_amb   covariance matrix of the ambiguities
            %
            %   approach    chose an approach amongo these:
            %                - 'lambda'     (yes) it's the only one available at the moment'
            %
            %   iar_method  Integer Ambiguities Resolution method
            %               Lambda approach:
            %                   1: ILS method based on search-and-shrink [DEFAULT]
            %                   2: ILS method based enumeration in search
            %                   3: integer rounding method
            %                   4: integer bootstrapping method
            %                   5: PAR with the input P0 of user-defined success rate
            %                   6: ILS method with Ratio Test (uses search-and shrink)
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
            %   [amb_fixed, is_fixed, l_fixed] = Fixer.fix(amb_float, C_amb_amb, approach, iar_method)
            
            this = Fixer();
            if nargin < 3 || isempty(approach)
                approach = 'lambda';
                iar_method = this.iar_method;
            end
            [amb_fixed, is_fixed, l_fixed] = this.fixAmbiguities(amb_float, C_amb_amb, approach, iar_method);
        end
    end
end
