function [Qzhat,Z,L,D,zhat,iZt] = decorrel_v2 (Qahat,ahat)
%DECORREL: Decorrelate a (co)va0.5.6
%
%     0.5.6
%
% This routine creates a decorrelated Q-matrix, by finding the
% Z-matrix and performing the corresponding transformation.
%
% The method is described in:
% The routine is based on Fortran routines written by Paul de Jonge (TUD)
% and on Matlab-routines written by Kai Borre.
% The resulting Z-matrix can be used as follows:
% zhat = Zt * ahat; \hat(z) = Z' * \hat(a);
% Q_\hat(z) = Z' * Q_\hat(a) * Z
%
% Input arguments:
%   Qahat: Variance-covariance matrix of ambiguities (original)
%   ahat:  Original ambiguities (optional)
%
% Output arguments:
%   Qzhat: Variance-covariance matrix of decorrelated ambiguitie
%   Z:     Z-transformation matrix
%   L:     L matrix (from LtDL-decomposition of Qzhat)
%   D:     D matrix (from LtDL-decomposition of Qzhat)
%   zhat:  Transformed ambiguities (optional)
%   iZt:   inv(Z')-transformation matrix
%
% ----------------------------------------------------------------------
% Function.: decorrel
% Date.....: 19-MAY-1999 / modified 12-APRIL-2012
% Author...: Peter Joosten / Sandra Verhagen
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
% ----------------------------------------------------------------------
%Tests on Inputs ahat and Qahat
%Is the Q-matrix symmetric

if ~isequal(Qahat-Qahat'<1E-8,ones(size(Qahat)))
  error ('Variance-covariance matrix is not symmetric!');
end

%Is the Q-matrix positive-definite?
if sum(eig(Qahat)>0) ~= size(Qahat,1)
  error ('Variance-covariance matrix is not positive definite!');
end
% -----------------------
% --- Initialisations ---
% -----------------------
n    = size(Qahat,1);
iZt  = eye(n);
c_end   = n - 1;
sw   = true;
% --------------------------
% --- LtDL decomposition ---
% --------------------------
[L,D] = ldldecom (Qahat);

% ------------------------------------------
% --- The actual decorrelation procedure ---
% ------------------------------------------
while sw
    c  = n;   % loop for column from n to 1
    sw = false;
    while ( ~sw ) && (c > 1)
        c = c - 1;  % the ith column
        if (c <= c_end)
            for r = c+1:n
                mu = round(L(r,c));
                if mu % if mu not equal to 0
                    L(r:n,c) = L(r:n,c) - mu * L(r:n,r);
                    iZt(:,r) = iZt(:,r) + mu * iZt(:,c);  % iZt is inv(Zt) matrix
                end
            end
        end
        delta = D(c) + D(c+1) * L(c+1, c) * L(c+1, c);
        if (delta < D(c+1))
            lambda       = D(c+1) * L(c+1,c) / delta;
            eta          = D(c) / delta;
            D(c)         = eta * D(c+1);
            D(c+1)       = delta;
            L(c:c+1,1:c-1) = [ -L(c+1,c) 1 ; eta lambda ] * L(c:c+1,1:c-1);
            L(c+1,c)     = lambda;
            % swap rows i and i+1
            L(c+2:n,c:c+1) = L(c+2:n,c+1:-1:c);
            iZt(:,c:c+1) = iZt(:,c+1:-1:c);
            c_end        = c;
            sw           = true;
        end
    end
end
% ---------------------------------------------------------------------
% --- Return the transformed Q-matrix and the transformation-matrix ---
% --- Return the decorrelated ambiguities, if they were supplied    ---
% ---------------------------------------------------------------------
Z = round(inv(iZt)');
Qzhat = Z' * Qahat * Z;
if nargin == 2 && nargout >= 5
  zhat = Z' * ahat;
end
return
% ----------------------------------------------------------------------
% End of routine: decorrel
% ----------------------------------------------------------------------
