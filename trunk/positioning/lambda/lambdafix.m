function [bcheck, acheck, Qzhat, Qbcheck] = lambdafix(bhat, ahat, Qbb, Qahat, Qba)

% SYNTAX:
%   [bcheck, acheck, Qzhat] = lambdafix(bhat, ahat, Qbb, Qahat, Qba);
%
% INPUT:
%   bhat  = position coordinates (float solution)
%   ahat  = ambiguities (float solution)
%   Qbb   = VCV-matrix (position block)
%   Qahat = VCV-matrix (ambiguity block)
%   Qba   = VCV-matrix (position-ambiguity covariance block)
%
% OUTPUT:
%   bcheck = output baseline (fixed or float depending on method and ratio test)
%   acheck = output ambiguities (fixed or float depending on method and ratio test)
%   Qzhat  = variance-covariance matrix of decorrelated ambiguities
%
% DESCRIPTION:
%   A wrapper for LAMBDA function to be used in goGPS.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Code contributed by Andrea Nardo
% Portions of code contributed by Eugenio Realini and Hendy F. Suhandri
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

global ratiotest mutest succ_rate fixed_solution IAR_method P0 mu flag_auto_mu flag_default_P0

Qbcheck=[];

if (flag_auto_mu)
    mu = [];
end

if (flag_default_P0)
    if (IAR_method == 5)
        P0 = 0.995;
    else
        P0 = 0.001;
    end
end

try
    % perform ambiguity resolution
    if (IAR_method == 0)
        %ILS enumeration (LAMBDA2)
        [U] = chol(Qahat); %compute cholesky decomposition
        Qahat = U'*U; %find back the vcm, now the off diag. comp. are identical
        [afixed,sqnorm,Qzhat,Z,D,L] = lambda_routine2(ahat,Qahat);
        % compute the fixed solution
        bcheck = bhat - Qba*cholinv(Qahat)*(ahat-afixed(:,1));
        Qbcheck = Qbb  - Qba*cholinv(Qahat)*Qba';
        acheck = afixed(:,1);
        % success rate
        Ps = prod(2*normcdf(0.5./sqrt(D))-1);
        %[up_bound, lo_bound] = success_rate(D,L,zeros(length(D)));
        
    elseif (IAR_method == 1 || IAR_method == 2)
        % ILS shrinking, method 1
        % ILS enumeration, method 2
        [afixed,sqnorm,Ps,Qzhat,Z]=LAMBDA(ahat,Qahat,IAR_method,'P0',P0,'mu',mu);
        % compute the fixed solution
        bcheck = bhat - Qba*cholinv(Qahat)*(ahat-afixed(:,1));
        acheck = afixed(:,1);
        Qbcheck = Qbb  - Qba*cholinv(Qahat)*Qba';
        
    elseif (IAR_method == 3 || IAR_method == 4)
        % Integer rounding, method 3
        % Integer bootstrapping, method 4
        [afixed,sqnorm,Ps,Qzhat,Z]=LAMBDA(ahat,Qahat,IAR_method,'P0',P0,'mu',mu);
        % compute the fixed solution
        bcheck = bhat - Qba*cholinv(Qahat)*(ahat-afixed(:,1));
        acheck = afixed(:,1);
        
    elseif (IAR_method == 5)
        % Partial Ambiguity Resolution, method 5
        %[afixed,sqnorm,Ps,Qahat,Z,nfx]=LAMBDA(ahat,Qahat,IAR_method,'P0',P0,'mu',mu);
        [afixed,sqnorm,Ps,Qzhat,Z,nfx]=LAMBDA(ahat,Qahat,IAR_method,'P0',P0,'mu',mu);
        %nfx = size(afixed, 1);
        %Z   = Z(:, 1:nfx);
        % in case of PAR afixed contains the decorrelated ambiguities
        if (nfx > 0)
            Qbz = Qba*Z;
            
            try
                %bcheck = bhat - Qbz *cholinv(Z'*Qahat*Z) * (Z'*ahat-afixed(:,1));
                bcheck = bhat - Qba *cholinv(Qahat) * (ahat-afixed(:,1));
            catch ME
                disp('Problems in PAR (lambdafix.m)');
                %keyboard;
            end
            
            % anyway we store the float ambiguities and their vcv-matrix... (to be improved)
            acheck = ahat;
            Qzhat = Qahat;
        else
            % keep float solution
            bcheck = bhat;
            acheck = ahat;
            Qzhat = Qahat;
        end
    end
catch
    % keep float solution
    bcheck = bhat;
    acheck = ahat;
    Qzhat = Qahat;
    
    fixed_solution = [fixed_solution 0];
    ratiotest = [ratiotest NaN];
    mutest    = [mutest NaN];
    succ_rate = [succ_rate NaN];
    
    return
end

% If IAR_method = 0 or IAR_method = 1 or IAR_method = 2 perform ambiguity validation through ratio test
if (IAR_method == 0 || IAR_method == 1 || IAR_method == 2)
    
    if (flag_auto_mu)
        if (1-Ps > P0)
            mu = ratioinv(P0,1-Ps,length(acheck));
        else
            mu = 1;
        end
    end
        
    ratio = sqnorm(1)/sqnorm(2);
    
    if ratio > mu
        % rejection; keep float baseline solution
        bcheck = bhat;
        acheck = ahat;
        
        fixed_solution = [fixed_solution 0];
    else
        fixed_solution = [fixed_solution 1];
    end
    
    ratiotest = [ratiotest ratio];
    mutest    = [mutest mu];
    
elseif (IAR_method == 5)

    if (nfx > 0)
        fixed_solution = [fixed_solution 1];
    else
        fixed_solution = [fixed_solution 0];
    end
else
    ratiotest = [ratiotest NaN];
    mutest    = [mutest NaN];
    
    fixed_solution = [fixed_solution 0];
end

succ_rate = [succ_rate Ps];
