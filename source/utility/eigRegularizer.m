function [ Nreg ] = eigRegularizer( N, k )

% -----------------------------------------------------------------
% matrizer regularization by checking eigenvalues
% -----------------------------------------------------------------
%
% N:   input square matrix
% k:   required condition number

% -----------------------------------------------------------------

[V, D] = eig(N); 					 % eigenvectors and eigenvalues
e = diag(D);                         % eigenvalues (ascending order)

% k = 1e6;                           % condition number
% pos = find( e < (e(end)/k) );
% e(pos) = e(end)/k;                 % regularized eigenvalues
% Nreg = V*diag(e)*V';               % regularized matrix

pos = e < (max(e) / k);
e(pos) = max(e) / k;                 % regularized eigenvalues
Nreg = V * diag(e) * V';             % regularized matrix
%Nreg = real(V*diag(e)*V');          % regularized matrix

% -----------------------------------------------------------------

% figure; imagesc(N); colorbar
% figure; imagesc(Nreg); colorbar
% figure; imagesc(N-Nreg); colorbar
