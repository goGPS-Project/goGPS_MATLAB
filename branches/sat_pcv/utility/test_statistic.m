function [identify, BNR_xhat, BNR_y] = test_statistic(Qy_hat, Qy, A, alpha, e_hat, Qe_hat)

% SYNTAX:
%   [identify, BNR_xhat, BNR_y] = test_statistic(Qy_hat, Qy, A, alpha, e_hat, Qe_hat);
%
% INPUT:
%   Qy_hat = estimated vcm of observation
%   Qy     = vcm of observation
%   A      = design matrix
%   alpha  = significant level
%   e_hat  = noise vector
%   Qe_hat = vcm of noise
%
% OUTPUT:
%   T        = T-test, detection
%   w        = w-test, identification
%   BNR_xhat = bias-to-noise ratio of estimate x
%   BNR_y    = bias-to_noise ratio of y 
%
% DESCRIPTION:
%

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Code contributed by Hendy F. Suhandri
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

%%defining some parameters
%-------------------------
m  = size(A,1);
n  = size(A,2);

q1 = 1;
qmn= m - n; 

%%Local Redundancy (r)
%----------------------
%r = 1 - (diag(Qyy_head)./diag(Qyy));

R = eye(m) - Qy_hat * (Qy^-1);
r = diag(R);

%%Mean local redundancy
%-----------------------
r_mean = mean(r);

%%Minimum detectable bias
%-------------------------
nabla = 3*(sqrt(diag(Qy)));
nabla_mean = 3*mean(sqrt(diag(Qy)));

%%Non-centrality parameter
%--------------------------
lambda_o= r.*nabla.^2./diag(Qy);
lambda_o_mean= mean(lambda_o);

%%Critical value (q = m-n)
%-------------------------
ka  = chi2inv(1 - alpha,qmn);
beta=ncx2cdf(ka,qmn,lambda_o_mean);
kb  = ncx2inv(beta,qmn,lambda_o_mean);

%%Critical value (q = 1)
%-------------------------
kA  = norminv(1 - alpha/2,q1);

%%T-test (for q = m-n)
T = e_hat' * Qy^(-1) * e_hat;

%%w-test (for q = 1)
for i= 1:m
    c    = zeros(m,1);
    c(i) = 1;
    w(i) = (c' * Qy^-1 * e_hat)^2./(c' * Qy^-1 * Qe_hat * Qy^-1 * c);
end
w = w';

identify =0;
if T > ka  
    %Detection
     disp('Ho rejected')

    %Identification
    identify = find (abs(w)>kA);
    fprintf('there is blunder in observation number: %d', identify);
end

%Bias to Noise ratio
BNR_xhat = (1 - r_mean./r_mean)*lambda_o_mean;
BNR_y = BNR_xhat  + lambda_o_mean;
